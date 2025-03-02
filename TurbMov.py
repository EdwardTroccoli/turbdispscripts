#!/usr/bin/env python
# -*- coding: utf-8 -*-
# written by Christoph Federrath, 2021-2025
# Modified by Edward Troccoli

import glob
import numpy as np
import timeit
import cfpack as cfp
import cfpack.hdfio as hdfio
import flashplotlib as fpl  # flashplotlib is assumed to handle plotting and saving

if __name__ == "__main__":
    # Start timing the process
    start_time = timeit.default_timer()
    path = "../Mach5-n128/"
    # Get all files matching the pattern Turb_slice_xy_*
    files = sorted(glob.glob(path+"Turb_slice_xy_*"))

    # Loop over each file
    for i, filen in enumerate(files):
        dat = cfp.read_ascii("Turb.dat")
        time = dat['01_time']
        t_turb = 0.1
        normalised_time = np.round(time/t_turb,4)
        # Read the magnetic field slices from the file using hdfio
        magx = hdfio.read(filen, "magx_slice_xy")
        magy = hdfio.read(filen, "magy_slice_xy")
        magz = hdfio.read(filen, "magz_slice_xy")        
        # Compute the magnetic energy density (or magnitude) emag
        emag = (magx**2 + magy**2 + magz**2) / (8 * np.pi)
       
        # Plot the computed emag using flashplotlib’s plotting routine.
        # The argument log=True tells the plot_map function to use a logarithmic scale.
        ret = cfp.plot_map(emag, log=True,cmap_label="Electromagnetic Energy Density (erg/m^3)",vmin=10**-2,vmax=10**2)
        cfp.plot(ax=ret.ax()[0], x=0, y=-0.05, text=normalised_time[i], normalised_coords=True, save=path/MovieFrames+f"test_{:04d}".format(i))
    # End timing and output the total processing time
    stop_time = timeit.default_timer()
    print("Processing time: {:.2f}s".format(stop_time - start_time))

