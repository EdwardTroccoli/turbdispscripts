#!/usr/bin/env python
# -*- coding: utf-8 -*-
# written by Edward Troccoli, 2025

import argparse
import numpy as np
import timeit
import os
import cfpack as cfp
from cfpack.defaults import *
from turblib import aver_pdf, write_pdf, read_pdf
from Globals import *
import glob


def compute_pdf(path, variable):
    if variable == "ekdr":
        vmin = -2e5
        vmax = +2e5
        bw = 100
    else:
        print("Variable not implemented.", error=True)
    # loop over FLASH dump files
    for d in range(20, 101):
        filename = "Turb_hdf5_plt_cnt_{:04d}".format(d)
        cfp.run_shell_command(f'mpirun -np 8 pdfs {path+filename} -dset {variable} -vmin {vmin} -vmax {vmax} -bw {bw}')


def plot_pdf(pdf_dat):
    y_label = 'PDF'
    cfp.plot(x=pdf_dat['col1'], y=pdf_dat['col3'])
    cfp.plot(xlabel=r'$\mathrm{x}$', ylabel=ylabel, save=fig_path+"tevol_"+f"{variable}_"+file_name+'.pdf', legend_loc='best')


if __name__ == "__main__":

    # Argument parser setup
    parser = argparse.ArgumentParser(description="Create and plot PDFs of different variables from FLASH simulation data.")
    var_choices = ["injr", "ekdr", "emdr", "ekin", "emag"]
    parser.add_argument("-ov", "--overwrite", action='store_true', default=False, help="Overwrite files")
    args = parser.parse_args()

    # Start timing the process
    start_time = timeit.default_timer()

    # vars = ["ekin", "emag", "injr", "ekdr", "emdr", "etdr"]
    vars = ["ekdr"]

    # loop over simulations
    for path in sim_paths:
        print(f'Working on: {path}', color='green')
        # loop over simulation variables
        for var in vars:
            pdf_aver_file = "aver_"+var+".pdf_data"
            if args.overwrite:
                compute_pdf(path, var) # compute the PDF by calling C++ 'pdf'
                pdf_files = glob.glob(path+"Turb_hdf5_plt_cnt_????_"+var+".pdf_data")
                aver_dat, header_aver = aver_pdf(pdf_files) # average the PDFs
                write_pdf(pdf_aver_file, aver_dat, header_aver) # write averaged PDF

            # plot the PDF for example
            if not os.path.isdir(fig_path):
                cfp.run_shell_command('mkdir '+fig_path)

            pdf_dat, pdf_header = read_pdf(pdf_aver_file) # read the PDF data
            stop()

            # plot_pdf(pdf_dat)

    # End timing and output the total processing time
    stop_time = timeit.default_timer()
    print("Processing time: {:.2f}s".format(stop_time - start_time))
