#!/usr/bin/env python
# -*- coding: utf-8 -*-
# written by Edward Troccoli and Christoph Federrath, 2025

import argparse
import timeit
import os, time
import cfpack as cfp
from cfpack.defaults import *
from turblib import aver_spect, write_spect, read_spect
from Globals import *
import glob

# computes spectra using C++ pdfs function
def compute_spectra(filename, out_path='./'):
    # Define expected output files
    extensions = ["_spect_vels.dat", "_spect_dset_ekdr.dat"]
    output_files = [filename.split('/')[-1]+"_spect_vels.dat", filename.split('/')[-1]+"_spect_dset_ekdr.dat"]
    # Check if file exists
    if not (os.path.exists(out_path+filename.split('/')[-1]+extensions[0]) and os.path.exists(out_path+filename.split('/')[-1]+extensions[1])):
        # run the spectra command
        cfp.run_shell_command(f'mpirun -np 8 spectra {filename} -types 0 1 -dsets ekdr')
        time.sleep(0.1)
        for ext in extensions:
            cfp.run_shell_command("mv "+filename+ext+" "+out_path)


# plotting function for spectras
def plot_spectra(dat, var):
    if var == "vels": ylabel=r'Power spectrum of $E_\mathrm{kin}$'
    if var == "ekdr": ylabel=r'Power spectrum of $\varepsilon_\mathrm{kin}$'
    y = 10**dat['col6']
    sigylo = y - 10**(dat['col6']-dat['col7'])
    sigyup = 10**(dat['col6']+dat['col7']) - y
    cfp.plot(x=dat['col1'], y=y, yerr=[sigylo,sigyup], shaded_err=True)
    cfp.plot(xlabel=r'Wavenumber $\mathrm{k}$', ylabel=ylabel, xlog=True, ylog=True,
            save=out_path+'aver_spectra'+ "_" + var + "_" + "M" +MachNumber[i] +'.pdf')


if __name__ == "__main__":

    # Argument parser setup
    parser = argparse.ArgumentParser(description="Create and plot PDFs of different variables from FLASH simulation data.")
    parser.add_argument("-ov", "--overwrite", action='store_true', default=False, help="Overwrite files")

    args = parser.parse_args()

    # Start timing the process
    start_time = timeit.default_timer()

    # loop over simulations
    for i, path in enumerate(sim_paths):

        print(f'Working on: {path}', color='green')

        # creates the figure output dir
        if not os.path.isdir(fig_path):
            cfp.run_shell_command('mkdir '+fig_path)

        # create output dir for spectra
        out_path = path + "spectra/"
        if not os.path.isdir(out_path):
            cfp.run_shell_command('mkdir '+out_path)

        # loop over simulation variables
        vars = ['ekdr', 'vels']
        for var in vars:
            spectra_aver_file = out_path+"aver_spectra_"+var+"_M"+MachNumber[i]+".dat"
            if not os.path.isfile(spectra_aver_file) or args.overwrite:
                for d in range(20, 101, 1):
                    filename = "Turb_hdf5_plt_cnt_{:04d}".format(d)
                    compute_spectra(path+filename, out_path=out_path) # compute the spectra by calling C++ 'spectra'
                if var == 'vels':
                    spectra_files = sorted(glob.glob(out_path+"*_spect_vels.dat"))
                elif var == 'ekdr':
                    spectra_files = sorted(glob.glob(out_path+"*_spect_dset_ekdr.dat"))
                aver_dat, header_aver = aver_spect(spectra_files) # average the spectras
                write_spect(spectra_aver_file, aver_dat, header_aver) # write averaged spectra

            # plot the spectra
            spectra_dat, spectra_header = read_spect(spectra_aver_file) # read the PDF data
            plot_spectra(spectra_dat, var)

    # End timing and output the total processing time
    stop_time = timeit.default_timer()
    print("Processing time: {:.2f}s".format(stop_time - start_time))
