#!/usr/bin/env python
# -*- coding: utf-8 -*-
# written by Edward Troccoli & Christoph Federrath, 2025

import h5py
import os, time
import glob
import argparse
import numpy as np
from scipy.stats import binned_statistic_2d
import dill
import gc, copy
from tqdm import tqdm
from cfpack import print, stop
import cfpack as cfp

# === define simulations to work on ===
sim_paths = ["../N256M0p2HDRe2500/", "../N256M5HDRe2500/","../N512M0p2HDRe2500/", "../N512M5HDRe2500/","../N1024M0p2HDRe2500/", "../N1024M5HDRe2500/","../N2048M0p2HDRe2500HP/", "../N2048M5HDRe2500HP/"]
# =====================================

# create figure output path
fig_path = "../Figures/"
if not os.path.isdir(fig_path):
    cfp.run_shell_command('mkdir '+fig_path)

# for getting global simulation parameters
def params(model_name):
    class ret:
        if 'M0p2' in model_name: Mach = 0.2
        if 'M5'   in model_name: Mach = 5
        if '2048' in model_name: N = 2048
        if '1024' in model_name: N = 1024
        if '512'  in model_name: N =  512
        if '256'  in model_name: N =  256
        t_turb = 0.5 / Mach
    return ret

# global 2D-PDF settings
vars_2Dpdf = [["dens", "ekdr"], ["vort", "ekdr"]]

# host settings
hostname = cfp.get_hostname()
if hostname.find("gadi") != -1:
    mpicmd = 'mpirun -np '
if hostname.find("setonix") != -1 or hostname.find("nid") != -1:
    mpicmd = 'srun -n '
if hostname.find("sng.lrz.de") != -1:
    mpicmd = 'mpiexec -n '

def compute_ekdr_size_fractal_dim_file(out_path, filename, overwrite=False, fix_neg=True):
    from cfpack.mpi import MPI, comm, myPE
    from cfflash import flashlib as fl
    norm_ekdr = params(out_path).t_turb / params(out_path).Mach**2 # EKDR norm
    fname_pkl = out_path+os.path.basename(filename)+"_ekdr_vs_size.pkl"
    if not os.path.isfile(fname_pkl) or overwrite:
        gg = fl.FlashGG(filename)
        dr = gg.D[0][0]
        half_diag = 0.5*np.sqrt(3)
        nbins_r = int(np.ceil(half_diag/dr))
        rad_bins = np.concatenate([[0.0], np.linspace(0.5*dr, (nbins_r-0.5)*dr, num=nbins_r)])
        ekdr_bins = cfp.symlogspace(-1e20, 1e20, lin_thresh=1e-30, num_lin=3, num=2*1000+3)
        bins = [rad_bins, ekdr_bins]
        minmax_obj = gg.GetMinMax("ekdr")
        centre = minmax_obj.max_loc
        print("min, max, max_loc = ", minmax_obj.min, minmax_obj.max, centre)
        bs = gg.binned_statistic(x="radius", y="ekdr", centre=centre, statistic='sum', bins=bins, use_hist = False)
        # get cumulative distribution and normalise by number of cells and by the EKDR unit
        if fix_neg:
            bs.y = np.nancumsum( bs.y * (bs.y > 0) )
        else:
            bs.y = np.nancumsum(bs.y)
        bs.y = bs.y * norm_ekdr / np.prod(gg.NMax) # normalise
        if myPE == 0: # only the master PE writes
            # save the data to file
            with open(fname_pkl, "wb") as fobj:
                print("Writing '"+fname_pkl+"'", color="magenta")
                dill.dump(bs, fobj)
    # read from file
    if MPI: comm.Barrier()
    print("Read '"+fname_pkl+"'", color="green")
    ret = dill.load(open(fname_pkl, "rb"))
    return ret

def get_ekdr_size_fractal_dim(path, overwrite=False):
    from cfpack.mpi import MPI, comm, myPE
    # compute ekdr vs. size
    print(f'Computing kinetic energy dissipation rate vs. size for fractal dimension analysis for: {path}', color='cyan')
    # create file output dir
    out_path = path + "FracDim/"
    if myPE == 0:
        if not os.path.isdir(out_path):
            cfp.run_shell_command('mkdir '+out_path)
    # output file
    fname_pkl = out_path+"aver_ekdr_vs_size.pkl"
    if not os.path.isfile(fname_pkl) or overwrite:
        bs_y = []
        dump_range = [20, 100]
        for d in range(dump_range[0], dump_range[1]+1, 1):
            filename = "Turb_hdf5_plt_cnt_{:04d}".format(d)
            bs = compute_ekdr_size_fractal_dim_file(out_path, path+filename, overwrite=overwrite)
            bs_y.append(bs.y)
        # setup a class to store edges and the averaged pdf data.
        class bsdat:
            def __init__(self, x_, y_, y_std_):
                self.x = x_
                self.y = y_
                self.y_std = y_std_ 
        if myPE == 0: # only the master rank writes to disk
            with open(fname_pkl, "wb") as fobj:
                print("Writing '"+fname_pkl+"' averged over dump range:", dump_range, color="magenta")
                y = np.mean(np.log10(np.stack(bs_y, axis=0)), axis=0)
                y_std = np.std(np.log10(np.stack(bs_y, axis=0)), axis=0)
                bsobj = bsdat(bs.xc, y, y_std)
                dill.dump(bsobj, fobj)
    if MPI: comm.Barrier()
    print("Read '"+fname_pkl+"'", color="green")
    bsdat = dill.load(open(fname_pkl, "rb"))
    # return averaged 2D PDF
    return bsdat

@cfp.timer_decorator
def compute_2d_pdf_file(out_path, filename, vars, bins, norms=[1.0,1.0], overwrite=False):
    import flashlib as fl
    from cfpack.mpi import MPI, comm, nPE, myPE
    print("Total number of MPI ranks = "+str(nPE))
    if MPI: comm.Barrier()
    fname_pkl = out_path+os.path.basename(filename)+"_2Dpdf_"+vars[0]+"_"+vars[1]+".pkl"
    if not os.path.isfile(fname_pkl) or overwrite:
        # read data
        gg = fl.FlashGG(filename, verbose=1)
        print("reading x and y data and computing binned_statistic_2d for vars", vars, "...", color="red")
        dsets = copy.deepcopy(vars)
        if dsets[0] == 'vort': dsets[0] = 'vorticity'
        MyBlocks = gg.GetMyBlocks(myPE, nPE) # domain decomposition
        counts_loc = []
        for b in tqdm(MyBlocks, disable=(myPE!=0), desc=f"[{myPE}]"): # loop over local list of block for the MPI rank
            x, y = gg.ReadBlockVar(b, dsets=dsets)
            x = x.flatten()*norms[0]
            y = y.flatten()*norms[1]
            counts_, x_edges_, y_edges_, binnum_ = binned_statistic_2d(x, y, np.ones_like(x, dtype=np.float32), statistic='count', bins=bins)
            counts_loc.append(counts_)
            del x; del y; del binnum_; del counts_
            gc.collect()
        # MPI reduction operation
        counts_loc = np.sum(np.array(counts_loc), axis=0)
        if MPI:
            counts_all = np.zeros(np.shape(counts_loc))
            comm.Reduce(counts_loc, counts_all, op=MPI.SUM) # master gets total sum
        else:
            counts_all = counts_loc
        if myPE == 0: # only the master PE writes
            class ret:
                # compute 2D counts
                counts, x_edges, y_edges, = counts_all, x_edges_, y_edges_
                # compute bin areas for normalization
                dx = np.diff(x_edges)[0]  # bin width in x
                dy = np.diff(y_edges)[0]  # bin width in y
                bin_area = dx * dy
                # normalize to get PDF (probability density)
                pdf = counts / (np.sum(counts) * bin_area)  # ensures sum(pdf * bin_area) = 1
            # save the data to file
            with open(fname_pkl, "wb") as fobj:
                print("Writing '"+fname_pkl+"'", color="magenta")
                dill.dump(ret, fobj)
    # read from file
    if MPI: comm.Barrier()
    print("Read '"+fname_pkl+"'", color="green")
    ret = dill.load(open(fname_pkl, "rb"))
    return ret

# computes or reads (if file already exists) 2D PDF for variables 'vars' in simulation 'path'
def get_2d_pdf(path, vars, overwrite=False):
    from cfpack.mpi import myPE
    # helper function to get normalisations
    def get_2Dpdf_norms(model_name, vars):
        norm_ekdr = params(model_name).t_turb / params(model_name).Mach**2
        norm_dens = 1.0
        norm_vort = 1.0 / (params(model_name).Mach * params(model_name).N)
        if vars == vars_2Dpdf[0]:
            norms = [norm_dens, norm_ekdr]
        if vars == vars_2Dpdf[1]:
            norms = [norm_vort, norm_ekdr]
        return norms
    # get 2D PDF
    print(f'\nComputing 2D-PDF for: {path}', color='cyan')
    # create file output dir
    out_path = path + "PDFs/"
    if not os.path.isdir(out_path):
        if myPE == 0: cfp.run_shell_command('mkdir '+out_path)
    # set binning
    if vars[1] == "ekdr": bins_y = np.logspace(-14, 6, 400)
    if vars[0] == "dens": bins_x = np.logspace(-6, 4, 200)
    if vars[0] == "vort": bins_x = np.logspace(-9, 1, 200)
    norms = get_2Dpdf_norms(path, vars)
    fname_pkl = out_path+"aver_2Dpdf_"+vars[0]+"_"+vars[1]+".pkl"
    if not os.path.isfile(fname_pkl) or overwrite:
        pdf_data = []
        dump_range = [20, 100]
        for d in range(dump_range[0], dump_range[1]+1, 1):
            filename = "Turb_hdf5_plt_cnt_{:04d}".format(d)
            po = compute_2d_pdf_file(out_path, path+filename, vars, bins=[bins_x,bins_y], norms=norms, overwrite=overwrite)
            pdf_data.append(po.pdf)
        # setup a class to store edges and the averaged pdf data.
        class pdat:
            pdf = np.mean(np.stack(pdf_data, axis=0), axis=0)
            x_edges = po.x_edges
            y_edges = po.y_edges
            variables = vars
        if myPE == 0: # only the master rank writes to disk
            with open(fname_pkl, "wb") as fobj:
                print("Writing '"+fname_pkl+"' averged over dump range:", dump_range, color="magenta")
                dill.dump(pdat, fobj)
    else:
        print("Read '"+fname_pkl+"'", color="green")
        pdat = dill.load(open(fname_pkl, "rb"))
    # return averaged 2D PDF
    return pdat

# function to compute vort, also checks if they already exist and runs only if they don't
def compute_vort_file(filename, ncpu=8, pixel=1024, overwrite=False):
    # first check if we need to generate anything
    run_derivative_var = False
    with h5py.File(filename, "r") as f:
        if not all(x in f for x in ['vorticity_x', 'vorticity_y', 'vorticity_z']):
            run_derivative_var = True
    if run_derivative_var or overwrite:
        # if plot file doesn't already contain vorticty components, generate them
        cfp.run_shell_command(mpicmd+str(ncpu)+' derivative_var '+filename+' -vort')
    else:
        cfp.print("Vorticity components already in '"+filename+"' -- skipping creation.", color='green')
    # taking slices
    for dim in ['x', 'y', 'z']:
        slice_filename = filename+"_vorticity_"+dim+"_slice_z.h5"
        outfile = os.path.dirname(slice_filename)+'/movie_files/'+os.path.basename(slice_filename)
        if not os.path.exists(outfile) or overwrite:
            cfp.run_shell_command(mpicmd+str(ncpu)+' projection '+filename+' -dset vorticity_'+dim+' -slice -pixel '+str(pixel)+' '+str(pixel))
            cfp.run_shell_command('mv '+slice_filename+' '+outfile)
        else:
            cfp.print("Slice file '"+outfile+"' already exists -- skipping.", color='green')

def compute_vort(overwrite=False):
    for sim_path in sim_paths:
        if params(sim_path).N == 2048:
            ncpu = 4096
        elif params(sim_path).N == 1024:
            ncpu = 512
        else:
            ncpu = 8
        dump_range = [20, 100]
        for d in range(dump_range[0], dump_range[1]+1, 1):
            plot_file = "Turb_hdf5_plt_cnt_{:04d}".format(d)
            compute_vort_file(sim_path+plot_file, ncpu=ncpu, pixel=params(sim_path).N, overwrite=overwrite)

# computes spectra using C++ pdfs function
def compute_spectra_file(filename, out_path='./', ncpu=8, overwrite=False):
    # Define expected output files
    extensions = ["_spect_vels.dat", "_spect_sqrtrho.dat", "_spect_dset_ekdr.dat"]
    # Check if file(s) exist
    call_spectra = False
    for ext in extensions:
        if not (os.path.exists(out_path+filename.split('/')[-1]+ext)) or overwrite:
            call_spectra = True
    if call_spectra:
        cfp.run_shell_command(mpicmd+str(ncpu)+' spectra '+filename+' -types 0 1 7 -dsets ekdr')
        time.sleep(0.1)
        for ext in extensions:
            cfp.run_shell_command("mv "+filename+ext+" "+out_path)
    else:
        print('Spectra files for '+filename+'already present - skipping.', color='green')

def compute_spectra(overwrite=False):
    for sim_path in sim_paths:
        if params(sim_path).N == 2048:
            ncpu = 512
        elif params(sim_path).N == 1024:
            ncpu = 64
        else:
            ncpu = 8
        plot_files = sorted(glob.glob(sim_path+"Turb_hdf5_plt_cnt_0???"))
        out_path = sim_path+'spectra/'
        if not os.path.isdir(out_path):
            cfp.run_shell_command('mkdir '+out_path)
        for plot_file in plot_files:
            compute_spectra_file(plot_file, out_path=out_path, ncpu=ncpu, overwrite=overwrite)

def clean_datfile(path, overwrite=False):
    print(f'Working on: {path}', color='green')
    # clean the .dat file
    datfile = path+'Turb.dat'
    if not os.path.isfile(datfile+'_cleaned') or overwrite:
        cfp.run_shell_command('flashlib.py datfile -clean '+datfile)
    out_path = path + "TimeEvol/"
    if not os.path.isdir(out_path):
        cfp.run_shell_command('mkdir '+out_path)

if __name__ == "__main__":

    # Argument parser setup
    parser = argparse.ArgumentParser(description="Main data analysis.")
    parser.add_argument("-ov", "--overwrite", action='store_true', default=False, help="Overwrite files")
    parser.add_argument("-clean_datfile", "--clean_datfile", action='store_true', default=False, help="Clean the .dat file(s).")
    parser.add_argument("-compute_vort", "--compute_vort", action='store_true', default=False, help="Compute vorticity.")
    parser.add_argument("-compute_2dpdfs", "--compute_2dpdfs", action='store_true', default=False, help="Compute 2D-PDFs.")
    parser.add_argument("-compute_spectra", "--compute_spectra", action='store_true', default=False, help="Compute spectra.")
    parser.add_argument("-compute_ekdr_size", "--compute_ekdr_size", action='store_true', default=False, help="Compute EKDR vs. size for fractal dimension.")
    args = parser.parse_args()

    if args.clean_datfile:
        for path in sim_paths:
            clean_datfile(path, overwrite=args.overwrite)

    if args.compute_vort:
        compute_vort(overwrite=args.overwrite)

    if args.compute_2dpdfs:
        for path in sim_paths:
            for vars in vars_2Dpdf: # loop over 2D PDF variables
                get_2d_pdf(path, vars, overwrite=args.overwrite)

    if args.compute_spectra:
        compute_spectra(overwrite=args.overwrite)

    if args.compute_ekdr_size:
        for path in sim_paths:
            bs = get_ekdr_size_fractal_dim(path, overwrite=args.overwrite)

