#!/usr/bin/env python
# -*- coding: utf-8 -*-
# written by Christoph Federrath, 2021-2025
# Modified by Edward Troccoli

import glob
from tkinter import Variable
import numpy as np
import timeit
import cfpack as cfp
import os
import cfpack.hdfio as hdfio
import flashplotlib as fpl 
from cfpack.defaults import *
from Globals import *
import argparse


def plot_variable(files, save_dir, variable):
    for i, filen in enumerate(files):
        # Read different variables from the file
        magx = hdfio.read(filen, "magx_slice_xy")
        magy = hdfio.read(filen, "magy_slice_xy")
        magz = hdfio.read(filen, "magz_slice_xy")
        emag = (magx ** 2 + magy ** 2 + magz ** 2) / (8 * np.pi)
        dens = hdfio.read(filen, "dens_slice_xy")
        emdr = hdfio.read(filen, "emdr_slice_xy")

        # Create a dictionary to select the correct variable
        data_dict = {
            "emag": emag,
            "dens": dens,
            "emdr": emdr
        }

        if variable not in data_dict:
            raise ValueError(f"Invalid variable name '{variable}'. Choose from {list(data_dict.keys())}")
        
        if variable == "emag":
            cmap_label = r"Electromagnetic Energy Density (erg/cm$^3$)"
            log = True
            vmin=10 ** -2
            vmax=10 ** 2
        elif variable == "dens":
            cmap_label = r"Density $\rho/\langle \rho\rangle$"
            log = True
            vmin=10 ** -2
            vmax=10 ** 2
        elif variable == "emdr":
            cmap_label = r"Kinetic Energy Dissipation Rate"
            log = False
            vmin=-10 ** 4
            vmax=10 ** 4
        else:
            return
            
        # Retrieve the correct data array
        data_to_plot = data_dict[variable]

        # Define formatted filename correctly
        save_filename = f"test_{variable}_{i:06d}.png"

        # Plot and save the image
        ret = cfp.plot_map(data_to_plot, log=log, cmap_label=cmap_label, vmin=vmin, vmax=vmax)

        # Roundabout way to add the time label
        cfp.plot(ax=ret.ax()[0], x=0.05, y=0.925, text=f"Time = {i/100}"+r"$t_\mathrm{turb}$", color='white', normalised_coords=True, save=os.path.join(save_dir, save_filename))


if __name__ == "__main__":
    # Start timing the process
    start_time = timeit.default_timer()
    path = "../Mach5-n128/"
    save_dir = os.path.join(path, "MovieFrames/") 

    # Argument parser setup
    parser = argparse.ArgumentParser(description="Plot different variables from simulation data.")
    parser.add_argument("--variable", type=str, required=True, help="Variable to plot (e.g., emag, dens, emdr)")
    # Parse the command-line arguments
    args = parser.parse_args()
     # Store the variable argument
    variable = args.variable

    # Get all files matching the pattern Turb_slice_xy_*
    files = sorted(glob.glob(path+"Turb_slice_xy_*"))
    
    #call plot function
    plot_variable(files,save_dir,variable)
    
    # End timing and output the total processing time
    stop_time = timeit.default_timer()
    print("Processing time: {:.2f}s".format(stop_time - start_time))

