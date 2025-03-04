#!/usr/bin/env python
# -*- coding: utf-8 -*-
# written by Christoph Federrath, 2021-2025
# Modified by Edward Troccoli

import glob
import numpy as np
import timeit
import cfpack as cfp
import os
import cfpack.hdfio as hdfio
import flashplotlib as fpl 

if __name__ == "__main__":
    # Start timing the process
    start_time = timeit.default_timer()
    path = "../Mach5-n128/"
    save_dir = os.path.join(path, "MovieFrames/") 
    # Get all files matching the pattern Turb_slice_xy_*
    files = sorted(glob.glob(path+"Turb_slice_xy_*"))
    #dat = cfp.read_ascii("Turb.dat")
    #time = dat['01_time']
    #t_turb = 0.1
    #normalised_time = np.round(time / t_turb, 4)
    # Loop over each file
    def emag_dens(files,save_dir):
        for i, filen in enumerate(files):

            # Read the magnetic field slices from the file using hdfio
            magx = hdfio.read(filen, "magx_slice_xy")
            magy = hdfio.read(filen, "magy_slice_xy")
            magz = hdfio.read(filen, "magz_slice_xy")

            # Compute the magnetic energy density (or magnitude) emag
            emag = (magx ** 2 + magy ** 2 + magz ** 2) / (8 * np.pi)

            # Define formatted filename correctly
            save_filename = "test_emag_dens_{:06d}.png".format(i)
            # Plot and save the image
            ret = cfp.plot_map(emag, log=True, cmap_label="Electromagnetic Energy Density (erg/m^3)", vmin=10 ** -2, vmax=10 ** 2)
            #Roundabout way to add the time label
            cfp.plot(ax=ret.ax()[0], x=0, y=1.025, text=f"Time = {i/100}$t_{{\\mathrm{{turb}}}}$", normalised_coords=True, save=os.path.join(save_dir, save_filename))
    def dens(files,save_dir):
        for i, filen in enumerate(files):

            # Read the magnetic field slices from the file using hdfio
            dens = hdfio.read(filen, "dens_slice_xy")
            # Compute the magnetic energy density (or magnitude) emag
            save_filename = "test_dens_{:06d}.png".format(i)
            # Plot and save the image
            ret = cfp.plot_map(dens, log=True, cmap_label="Density (kg/m^3)", vmin=10 ** -2, vmax=10 ** 2)
            #Roundabout way to add the time label
            cfp.plot(ax=ret.ax()[0], x=0, y=1.025, text=f"Time = {i/100}$t_{{\\mathrm{{turb}}}}$", normalised_coords=True, save=os.path.join(save_dir, save_filename))
            
    def emdr(files,save_dir):
        for i, filen in enumerate(files):

            # Read the magnetic field slices from the file using hdfio
            emdr = hdfio.read(filen, "emdr_slice_xy")
            # Compute the magnetic energy density (or magnitude) emag
            save_filename = "test_emdr_{:06d}.png".format(i)
            # Plot and save the image
            ret = cfp.plot_map(emdr, log=True, cmap_label="Power (W)", vmin=10 ** -2, vmax=10 ** 2)
            #Roundabout way to add the time label
            cfp.plot(ax=ret.ax()[0], x=0, y=1.025, text=f"Time = {i/100}$t_{{\\mathrm{{turb}}}}$", normalised_coords=True, save=os.path.join(save_dir, save_filename))
            
    # End timing and output the total processing time
    stop_time = timeit.default_timer()
    print("Processing time: {:.2f}s".format(stop_time - start_time))

