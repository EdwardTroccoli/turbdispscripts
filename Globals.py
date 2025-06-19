#!/usr/bin/env python
# -*- coding: utf-8 -*-
# written by Edward Troccoli & Christoph Federrath, 2025

from cfpack.defaults import *
import cfpack as cfp
import os, time
import h5py
import glob
import argparse

# create figure output path
fig_path = "../Figures/"
if not os.path.isdir(fig_path):
    cfp.run_shell_command('mkdir '+fig_path)

# sim_paths = ["N256M0p2HDRe2500", "N512M0p2HDRe2500", "N1024M0p2HDRe2500", "N256M5HDRe2500", "N512M5HDRe2500", "N1024M5HDRe2500"]
#sim_paths = ["../N2048M0p2HDRe2500HP/", "../N2048M5HDRe2500HP/"]
sim_paths = ["../N1024M0p2HDRe2500/", "../N1024M5HDRe2500/"]

def params(model_name):
    class ret:
        if 'M0p2' in model_name: Mach = 0.2
        if 'M5'   in model_name: Mach = 5
        if '2048' in model_name: N = 2048
        if '1024' in model_name: N = 1024
        if '512'  in model_name: N =  512
        if '256'  in model_name: N =  256
        t_turb = 0.5 / Mach
        # do the color selection here as well, based on the model name / path; needs implementation...
        color = 'black'
    return ret

# function to compute vort, also checks if they already exist and runs only if they don't
def compute_vort_file(filename, ncpu=8, pixel=1024, overwrite=False):
    # first check if we need to generate anything
    run_derivative_var = False
    with h5py.File(filename, "r") as f:
        if not all(x in f for x in ['vorticity_x', 'vorticity_y', 'vorticity_z']):
            run_derivative_var = True
    if run_derivative_var or overwrite:
        # if plot file doesn't already contain vorticty components, generate them
        cfp.run_shell_command('mpirun -np '+str(ncpu)+' derivative_var '+filename+' -vort')
    else:
        cfp.print("Vorticity components already in '"+filename+"' -- skipping creation.", color='green')
    # taking slices
    for dim in ['x', 'y', 'z']:
        slice_filename = filename+"_vorticity_"+dim+"_slice_z.h5"
        outfile = os.path.dirname(slice_filename)+'/movie_files/'+os.path.basename(slice_filename)
        if not os.path.exists(outfile) or overwrite:
            cfp.run_shell_command('mpirun -np '+str(ncpu)+' projection '+filename+' -dset vorticity_'+dim+' -slice -pixel '+str(pixel)+' '+str(pixel))
            cfp.run_shell_command('mv '+slice_filename+' '+outfile)
        else:
            cfp.print("Slice file '"+outfile+"' already exists -- skipping.", color='green')

def compute_vort(overwrite=False):
    for sim_path in sim_paths:
        if params(sim_path).N == 2048:
            ncpu = 512
        else:
            ncpu = 8
        plot_files = sorted(glob.glob(sim_path+"Turb_hdf5_plt_cnt_0???"))
        for plot_file in plot_files:
            compute_vort_file(plot_file, ncpu=ncpu, pixel=params(sim_path).N, overwrite=overwrite)

# computes spectra using C++ pdfs function
def compute_spectra_file(filename, out_path='./', ncpu=8, overwrite=False):
    # Define expected output files
    extensions = ["_spect_vels.dat", "_spect_dset_ekdr.dat"]
    output_files = [filename.split('/')[-1]+"_spect_vels.dat", filename.split('/')[-1]+"_spect_dset_ekdr.dat"]
    # Check if file exists
    if not (os.path.exists(out_path+filename.split('/')[-1]+extensions[0]) and os.path.exists(out_path+filename.split('/')[-1]+extensions[1])) or overwrite:
        # run the spectra command
        cfp.run_shell_command(f'mpirun -np 64 spectra {filename} -types 0 1 -dsets ekdr')
        time.sleep(0.1)
        for ext in extensions:
            cfp.run_shell_command("mv "+filename+ext+" "+out_path)

def compute_spectra(overwrite=False):
    for sim_path in sim_paths:
        if params(sim_path).N == 2048:
            ncpu = 512
        else:
            ncpu = 8
        plot_files = sorted(glob.glob(sim_path+"Turb_hdf5_plt_cnt_0???"))
        for plot_file in plot_files:
            compute_spectra_file(plot_file, ncpu=ncpu, overwrite=overwrite)

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
    parser.add_argument("-compute_spectra", "--compute_spectra", action='store_true', default=False, help="Compute spectra files.")
    args = parser.parse_args()

    if args.clean_datfile:
        for path in sim_paths:
            clean_datfile(path, overwrite=args.overwrite)

    if args.compute_vort:
        compute_vort(overwrite=args.overwrite)

    if args.compute_spectra:
        compute_spectra(overwrite=args.overwrite)
