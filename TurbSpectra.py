#!/usr/bin/env python
# -*- coding: utf-8 -*-
# written by Edward Troccoli and Christoph Federrath, 2025

import argparse
import timeit
import os
import cfpack as cfp
from cfpack.defaults import *
from turblib import aver_spect, write_spect, read_spect
from Globals import *
import glob
import flashlib as fl

# computes spectra using C++ pdfs function
def compute_spectra(filename):
    # Define expected output files
    base = os.path.splitext(filename)[0]
    output_file_1 = f"{base}_spect_vels.dat"
    output_file_2 = f"{base}_spect_dset_ekdr.dat"

    # Check if both files exist
    if not (os.path.exists(output_file_1) and os.path.exists(output_file_2)):
        # If either file is missing, run the spectra command
        cfp.run_shell_command(f'mpirun -np 8 spectra {filename} -types 0 1 -dsets ekdr')

# plotting function for spectras
def plot_spectra(pdf_dat):
    if var == "ekdr":
        cfp.plot(x=pdf_dat['col2'], y=pdf_dat['col8'])
        cfp.plot(xlabel='Kinetic Energy Dissipation Rate', ylabel='PDF of Kinetic Energy Dissipation Rate', 
                 ylog=True, save=out_path+'aver_1DPDF_'+var+'.pdf')
    elif var == "vels":
        cfp.plot(x=pdf_dat['col1'], y=pdf_dat['col5'])
        cfp.plot(xlabel=r'Wavenumber $\mathrm{k}$', ylabel='Velocity',
                xlog=True, show=True)


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
        
        # loop over simulation variables
        vars = ['vels']
        for var in vars:
            spectra_aver_file = "aver_spectra_"+var+".dat"
            if args.overwrite:
                for d in range(20, 23):
                    filename = "Turb_hdf5_plt_cnt_{:04d}".format(d)
                    compute_spectra(path+filename) # compute the spectra by calling C++ 'spectra'
                if var == 'vels':
                    spectra_files = sorted(glob.glob(path + "*_spect_vels.dat"))
                elif var == 'ekdr':
                    spectra_files = sorted(glob.glob(path + "*_spect_dset_ekdr.dat"))
                aver_dat, header_aver = aver_spect(spectra_files) # average the spectras
                write_spect(spectra_aver_file, aver_dat, header_aver) # write averaged spectra

        # plot the spectra
        out_path = path + "spectras/"
        if not os.path.isdir(out_path):
            cfp.run_shell_command('mkdir '+out_path)

        spectra_dat, spectra_header = read_spect(spectra_aver_file) # read the PDF data
        plot_spectra(spectra_dat)
    # End timing and output the total processing time
    stop_time = timeit.default_timer()
    print("Processing time: {:.2f}s".format(stop_time - start_time))
