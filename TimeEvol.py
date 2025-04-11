#!/usr/bin/env python
# -*- coding: utf-8 -*-
# written by Edward Troccoli, 2025

import argparse
import numpy as np
import timeit
import os
import cfpack as cfp
from cfpack.defaults import *
from Globals import *


def plot_var(variable):

    if not os.path.isdir(fig_path):
        cfp.run_shell_command('mkdir '+fig_path)

    for i, path in enumerate(sim_paths):

        print(f'Working on: {path}', color='green')

        dat = cfp.read_ascii(path+"Turb.dat")
        time = dat['01_time'] / t_turb[i]
        
        if variable == "vstd":
            var = dat['#14_rms_velocity']
            ylabel = r'$\sigma_v$'
            cfp.plot(x=time, y=var, label=path[3:-1], color=color[i])
        elif variable == "ekin":
            var = dat['#10_E_kinetic']
            ylabel = r'$E_\mathrm{kin}$'
            cfp.plot(x=time, y=var, label=path[3:-1], color=color[i])
        elif variable == "emag":
            var = dat['#12_E_magnetic']
            ylabel = r'$E_\mathrm{mag}$'
            cfp.plot(x=time, y=var, label=path[3:-1], color=color[i])
        elif variable == "injr":
            var = dat['#41_injection_rate']
            ylabel = r'injection rate'
            cfp.plot(x=time, y=var, label=path[3:-1], color=color[i])
        elif variable == "ekdr":
            var = dat['#42_ekin_diss_rate']
            ylabel = r'$\varepsilon_{\textrm{kin}}$'
            cfp.plot(x=time, y=var, label=path[3:-1], color=color[i])
        elif variable == "emdr":
            var = dat['#43_emag_diss_rate']
            ylabel = r'magnetic energy dissipation rate'
            cfp.plot(x=time, y=var, label=path[3:-1], color=color[i])
        elif variable == "etdr":
            var = dat['#42_ekin_diss_rate'] + dat['#43_emag_diss_rate']
            ylabel = r'total energy dissipation rate'
            cfp.plot(x=time, y=var, label=path[3:-1], color=color[i])
        elif variable == "ired":
            injr = dat['#41_injection_rate']
            ekdr = dat['#42_ekin_diss_rate']
            ylabel = r'$\varepsilon_{\textrm{kin}}$ and Injection Rate'
            cfp.plot(x=time, y=injr, label=r'$\varepsilon_{\textrm{inj}}$', color=color[i])
            cfp.plot(x=time, y=ekdr, label=r'$\varepsilon_{\textrm{kin}}$', color=color[i+1])
        
    out_path = path + "TimeEvol/"
    if not os.path.isdir(out_path):
        cfp.run_shell_command('mkdir '+out_path)
    cfp.plot(xlabel=r'$t/t_\mathrm{turb}$', ylabel=ylabel, save=out_path+"tevol_"+f"{variable}_.pdf", legend_loc='best')


if __name__ == "__main__":

    # Argument parser setup
    parser = argparse.ArgumentParser(description="Plot different variables from simulation data.")
    var_choices = ["vstd", "ekin", "emag", "injr", "ekdr", "emdr", "etdr", "ired"]
    parser.add_argument("-v", "--variable", nargs='*', choices=var_choices, required=True, help="Variable to plot; choice of "+str(var_choices))
    # Parse the command-line arguments
    args = parser.parse_args()

    # Start timing the process
    start_time = timeit.default_timer()

    for var in args.variable:
        plot_var(var)

    # End timing and output the total processing time
    stop_time = timeit.default_timer()
    print("Processing time: {:.2f}s".format(stop_time - start_time))
