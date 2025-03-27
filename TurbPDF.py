#!/usr/bin/env python
# -*- coding: utf-8 -*-
# written by Edward Troccoli, 2025

import argparse
import numpy as np
import timeit
import os
import cfpack as cfp
from cfpack.defaults import *
from turblib import read_pdf
from Globals import *
import glob


def compute_pdf(variable,file_name):

    if not os.path.isdir(fig_path):
        cfp.run_shell_command('mkdir '+fig_path)
    
    for i, path in enumerate(sim_paths):

        print(f'Working on: {path}', color='green')

        if variable == "emdr":
            vmin = 1
            vmax = 1e4
            stop()
            cfp.run_shell_command(f'pdfs {file_name} -dset {variable} -bw 100 -vmin vmin -vmax vmax')
            
            pdf_file = glob.glob(f'{file_name}_{variable}*')
            var = read_pdf(pdf_file[0])
            data = var[0]
        else:
            print("Variable not implemented.", error=True)
        

    return data
        

def plot_pdf(pdf):
    y_label = 'PDF'

    cfp.plot(x=pdf['col1'], y=pdf['col2'])
    cfp.plot(xlabel=r'$\mathrm{x}$', ylabel=ylabel, save=fig_path+"tevol_"+f"{variable}_"+file_name+'.pdf', legend_loc='best')


if __name__ == "__main__":

    # Argument parser setup
    parser = argparse.ArgumentParser(description="Plot different variables from simulation data.")
    var_choices = ["vstd", "ekin", "emag", "injr", "ekdr", "emdr", "etdr"]
    name_choices = ["../Mach5-n128/AlfvenMach6/","../Mach5-n128/AlfvenMach10/","..."]
    parser.add_argument("-v", "--variable", nargs='*', choices=var_choices, required=True, help="Variable to plot; choice of "+str(var_choices))
    parser.add_argument("-n", "--filename", required=True, help="Give the plotfile a name "+str(name_choices))
    # Parse the command-line arguments
    args = parser.parse_args()

    # Start timing the process
    start_time = timeit.default_timer()


    for var in args.variable:
        plot_pdf(compute_pdf(var, args.filename))

    # End timing and output the total processing time
    stop_time = timeit.default_timer()
    print("Processing time: {:.2f}s".format(stop_time - start_time))
