#!/usr/bin/env python
# -*- coding: utf-8 -*-
# written by Edward Troccoli & Christoph Federrath, 2025

import glob
import numpy as np
import os
import timeit
import cfpack as cfp
import cfpack.hdfio as hdfio
from cfpack.defaults import *
from Globals import *
import argparse
import matplotlib.pyplot as plt


def plot_variable(files, out_path, variable, tturb):

    # loop over plot files
    for i, filen in enumerate(files):

        # defaults
        log = True
        cmap = 'afmhot'
        xlabel = None
        ylabel = None
        if '0p2' in out_path:
            M = '0.2'
            ekdr_max = 0.25
            dens_min = 0.95
            dens_max = 1.035
            ekin_min = 0
            ekin_max = 1e-1
            log = False
        elif '5' in out_path:
            M = '5'
            ekdr_max = 2e2
            dens_min = 1e-2
            dens_max = 1e2
            ekin_min = 1e-2
            ekin_max = 1e2
            log = True
        if variable == "dens":
            remove_x_ticks = False
            log = log
            vmin = dens_min
            vmax = dens_max
            var = hdfio.read(filen, "dens_slice_xy")
            if '0p2' in out_path:
                cmap_label = None
                ylabel = r'$y$'
                xlabel = r'$x$'
                remove_y_ticks = False
            elif '5' in out_path:
                cmap_label = r'Density $\rho/\langle\rho\rangle$'
                xlabel = r'$x$'
                remove_y_ticks = True
        elif variable == "ekin":
            remove_x_ticks = True
            cmap = 'RdPu'
            vmin = ekin_min
            vmax = ekin_max
            dens = hdfio.read(filen, "dens_slice_xy")
            velx = hdfio.read(filen, "velx_slice_xy")
            vely = hdfio.read(filen, "vely_slice_xy")
            velz = hdfio.read(filen, "velz_slice_xy")
            var = 0.5 * dens * (velx ** 2 + vely ** 2 + velz ** 2)
            if '0p2' in out_path:
                cmap_label = None
                ylabel = r'$y$'
                remove_y_ticks = False
                log = False
            elif '5' in out_path:
                cmap_label = r"Kinetic energy $E_{\textrm{kin}}/\langle\rho\rangle\, c_{\textrm{s}}^2$"
                remove_y_ticks = True
        elif variable == "ekdr":
            remove_x_ticks = True
            cmap = 'BuPu'
            log = False
            vmin = 0
            vmax = ekdr_max
            if '0p2' in out_path:
                cmap_label = None
                ylabel = r'$y$'
                remove_y_ticks = False
                var = hdfio.read(filen, "ekdr_slice_xy")*tturb
            elif '5' in out_path:
                cmap_label = r"Dissipation rate $\varepsilon_{\textrm{kin}}//\langle\rho\rangle\, c_{\textrm{s}}^2\,t_{\textrm{turb}}^{-1}$"
                var = hdfio.read(filen, "ekdr_slice_xy")*tturb
                remove_y_ticks = True
        elif variable == "emag":
            magx = hdfio.read(filen, "magx_slice_xy")
            magy = hdfio.read(filen, "magy_slice_xy")
            magz = hdfio.read(filen, "magz_slice_xy")
            var = (magx ** 2 + magy ** 2 + magz ** 2) / (8 * np.pi)
            cmap_label = r"Magnetic energy density"
        elif variable == "emdr":
            var = hdfio.read(filen, "emdr_slice_xy")
            cmap_label = r"Magnetic energy dissipation rate"
            log = False
            vmin = 0
            vmax = 5
        else:
            print("Variable not implemented.", error=True)

        # Define formatted filename correctly
        #out_file = out_path+f"frame_{variable}_{i:06d}.png"
        out_file = out_path+f"frame_{variable}_000250.pdf"

        # Plot and save the image
        ret = cfp.plot_map(var, log=log, cmap_label=cmap_label, cmap=cmap, xlim=[0,1], ylim=[0,1], aspect_data='equal', vmin=vmin, vmax=vmax)
        ax = ret.ax()[0]
        ax.text(0.05, 0.95, r"$\mathcal{M} = $"+M, transform=ax.transAxes,
        fontsize=14, color='white', verticalalignment='top',
        bbox=dict(boxstyle="round,pad=0.3", facecolor='gray', alpha=0.5))
        if remove_x_ticks == True:
            ax.set_xticklabels([])
        if remove_y_ticks == True:
            ax.set_yticklabels([])
        #time = hdfio.read(filen, "time")[0] / tturb
        #time_str = cfp.round(time, 3, str_ret=True)
        cfp.plot(ax=ret.ax()[0], x=0.05, y=0.925, xlabel=xlabel, ylabel=ylabel,  color='white', normalised_coords=True, save=out_file)#text=r"$t/t_\mathrm{turb}="+time_str+r"$",

if __name__ == "__main__":

    # Argument parser setup
    parser = argparse.ArgumentParser(description="Plot different variables from simulation data.")
    var_choices = ["dens", "ekin", "ekdr", "emag", "emdr"]
    parser.add_argument("-v", "--variable", choices=var_choices, required=True, nargs= '*', help="Variable to plot; choice of "+str(var_choices))
    # Parse the command-line arguments
    args = parser.parse_args()

    # Start timing the process
    start_time = timeit.default_timer()

    # loop over simulations
    for i, path in enumerate(sim_paths):

        print(f'Working on: {path}', color='green')

        for var in args.variable:
            out_path = path + "movie_files/"
            if not os.path.isdir(out_path):
                cfp.run_shell_command('mkdir '+out_path)

            # Get all files matching the pattern Turb_slice_xy_*
            #files = sorted(glob.glob(out_path+"Turb_slice_xy_*"))
            files = sorted(glob.glob(out_path+"Turb_slice_xy_000250"))
            
            # call plot function
            plot_variable(files, out_path, var, t_turb[i])

    # End timing and output the total processing time
    stop_time = timeit.default_timer()
    print("Processing time: {:.2f}s".format(stop_time - start_time))

