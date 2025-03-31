#!/usr/bin/env python
# -*- coding: utf-8 -*-
# written by Edward Troccoli and Christoph Federrath, 2025

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
        vmin = -1e4
        vmax = +1e4
        bw = 20
    elif variable == "emdr":
        vmin = -1e4
        vmax = +1e4
        bw = 20
    elif variable == "emag":
        vmin = -1e4
        vmax = +1e4
        bw = 20
    elif variable == "ekin":
        vmin = -1e4
        vmax = +1e4
        bw = 20
    elif variable == "injr":
        vmin = -1e4
        vmax = +1e4
        bw = 20
    else:
        print("Variable not implemented.", error=True)
    # loop over FLASH dump files
    for d in range(20, 101):
        filename = "Turb_hdf5_plt_cnt_{:04d}".format(d)
        cfp.run_shell_command(f'mpirun -np 8 pdfs {path+filename} -dset {variable} -vmin {vmin} -vmax {vmax} -bw {bw}')


def plot_pdf(pdf_dat):
    if var == "ekdr":
        cfp.plot(x=pdf_dat['col1'], y=pdf_dat['col3'])
        cfp.plot(xlabel='Kinetic Energy Dissipation Rate', ylabel='PDF of Kinetic Energy Dissipation Rate', xlim = [-1e4,1e4], ylog = True, save=fig_path+"aver_"+var+'.pdf', legend_loc='best')
    elif var == "emdr":
        cfp.plot(x=pdf_dat['col1'], y=pdf_dat['col3'])
        cfp.plot(xlabel="Magnetic Energy Dissipation Rate", ylabel='PDF of Magnetic Energy Dissipation Rate', xlim = [-1e4,1e4], ylog = True, save=fig_path+"aver_"+var+'.pdf', legend_loc='best') 
    elif var == "injr":
        cfp.plot(x=pdf_dat['col1'], y=pdf_dat['col3'])
        cfp.plot(xlabel="Energy Injection Rate", ylabel='PDF of Energy Injection Rate', xlim = [-1e4,1e4], ylog = True, save=fig_path+"aver_"+var+'.pdf', legend_loc='best') 
    elif var == "emag":
        cfp.plot(x=pdf_dat['col1'], y=pdf_dat['col3'])
        cfp.plot(xlabel="Magnetic Energy", ylabel='PDF of Magnetic Energy', xlim = [0,230], ylog = True, save=fig_path+"aver_"+var+'.pdf', legend_loc='best') 
    elif var == "ekin":
        cfp.plot(x=pdf_dat['col1'], y=pdf_dat['col3'])
        cfp.plot(xlabel="Kinetic Energy", ylabel='PDF of Kinetic Energy', xlim = [0,1550], ylog = True, save=fig_path+"aver_"+var+'.pdf', legend_loc='best')  

if __name__ == "__main__":

    # Argument parser setup
    parser = argparse.ArgumentParser(description="Create and plot PDFs of different variables from FLASH simulation data.")
    var_choices = ["injr", "ekdr", "emdr", "ekin", "emag"]
    parser.add_argument("-ov", "--overwrite", action='store_true', default=False, help="Overwrite files")
    args = parser.parse_args()

    # Start timing the process
    start_time = timeit.default_timer()

    vars = ["ekin", "emag", "injr", "ekdr", "emdr"]

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
            plot_pdf(pdf_dat)

    # End timing and output the total processing time
    stop_time = timeit.default_timer()
    print("Processing time: {:.2f}s".format(stop_time - start_time))
