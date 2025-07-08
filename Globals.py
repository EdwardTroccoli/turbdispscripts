#!/usr/bin/env python
# -*- coding: utf-8 -*-
# written by Edward Troccoli & Christoph Federrath, 2025

import cfpack as cfp
import h5py
import os, time
import glob
import argparse
import numpy as np
from scipy.stats import binned_statistic_2d
import dill
import flashlib as fl
import gc, copy
from tqdm import tqdm
from cfpack.defaults import *
import cfpack as cfp

# === define simulations to work on ===
sim_paths = ["../N2048M0p2HDRe2500HP/", "../N2048M5HDRe2500HP/"]
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

@cfp.timer_decorator
def compute_2d_pdf_file(out_path, filename, vars, bins, norms=[1.0,1.0], overwrite=False):
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
    if vars[1] == "ekdr": bins_y = np.logspace(-8, 6, 250)
    if vars[0] == "dens": bins_x = np.logspace(-4, 3, 250)
    if vars[0] == "vort": bins_x = np.logspace(-6, 1, 250)
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
        if myPE == 0:
            with open(fname_pkl, "wb") as fobj:
                print("Writing '"+fname_pkl+"' averged over dump range:", dump_range, color="magenta")
                dill.dump(pdat, fobj) # only the master rank writes to disk
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
            ncpu = 2048
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
