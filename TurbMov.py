#!/usr/bin/env python
# -*- coding: utf-8 -*-
# written by Edward Troccoli & Christoph Federrath, 2025

import glob
from tkinter import Variable
import numpy as np
import os
import timeit
import cfpack as cfp
import cfpack.hdfio as hdfio
from cfpack.defaults import *
from Globals import *
import argparse


def plot_variable(files, out_path, variable, tturb):

    # loop over plot files
    for i, filen in enumerate(files):

        # defaults
        log = True
        vmin = 1e-2
        vmax = 1e2

        if variable == "dens":
            var = hdfio.read(filen, "dens_slice_xy")
            cmap_label = r"Density ($\rho/\langle\rho\rangle$)"
        elif variable == "ekin":
            dens = hdfio.read(filen, "dens_slice_xy")
            velx = hdfio.read(filen, "velx_slice_xy")
            vely = hdfio.read(filen, "vely_slice_xy")
            velz = hdfio.read(filen, "velz_slice_xy")
            var = 0.5 * dens * (velx ** 2 + vely ** 2 + velz ** 2)
            cmap_label = r"Kinetic energy density"
        elif variable == "ekdr":
            var = hdfio.read(filen, "ekdr_slice_xy")
            cmap_label = r"Kinetic energy dissipation rate"
            log = False
            vmin = -1e4
            vmax = 1e4
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
            vmin = -1e4
            vmax = 1e4
        else:
            print("Variable not implemented.", error=True)

        # Define formatted filename correctly
        out_file = out_path+f"frame_{variable}_{i:06d}.png"

        # Plot and save the image
        ret = cfp.plot_map(var, log=log, cmap_label=cmap_label, cmap='afmhot', xlabel=r"$x$", ylabel=r"$y$", xlim=[0,1], ylim=[0,1], aspect_data='equal', vmin=vmin, vmax=vmax)

        # Roundabout way to add the time label
        time = hdfio.read(filen, "time")[0] / tturb
        time_str = cfp.round(time, 3, str_ret=True)
        cfp.plot(ax=ret.ax()[0], x=0.05, y=0.925, xlabel=r"$x$", ylabel=r"$y$", text=r"$t/t_\mathrm{turb}="+time_str+r"$", color='white', normalised_coords=True, save=out_file)


if __name__ == "__main__":

    # Argument parser setup
    parser = argparse.ArgumentParser(description="Plot different variables from simulation data.")
    var_choices = ["dens", "ekin", "ekdr", "emag", "emdr"]
    parser.add_argument("-v", "--variable", choices=var_choices, required=True, help="Variable to plot; choice of "+str(var_choices))
    # Parse the command-line arguments
    args = parser.parse_args()

    # Start timing the process
    start_time = timeit.default_timer()

    # loop over simulations
    for i, path in enumerate(sim_paths):

        out_path = path + "MovieFrames/"
        if not os.path.isdir(out_path):
            cfp.run_shell_command('mkdir '+out_path)

        # Get all files matching the pattern Turb_slice_xy_*
        files = sorted(glob.glob(path+"Turb_slice_xy_*"))
        #files = sorted(glob.glob(path+"Turb_slice_xy_00000"))

        # call plot function
        plot_variable(files, out_path, args.variable, t_turb[i])

    # End timing and output the total processing time
    stop_time = timeit.default_timer()
    print("Processing time: {:.2f}s".format(stop_time - start_time))

