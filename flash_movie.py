#!/usr/bin/env python
# -*- coding: utf-8 -*-
# written by Christoph Federrath, 2021-2025

import glob
from cfpack import stop, print
import timeit
import flashplotlib as fpl
import argparse


# ===== the following applies in case we are running this in script mode =====
if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Create a movie out of projections/slices.')
    #parser.add_argument("-i", "--i", dest='files', nargs='*', help="HDF5 input file(s)")
    args0 = parser.parse_args()

    # time the script
    start_time = timeit.default_timer()

    # files
    files = sorted(glob.glob("Turb_hdf5_plt_cnt_00??"))

    # animate
    values = fpl.animate(nframes=10, value_beg=0.0, value_end=360)

    # animate camera view an move
    #view_start = np.array([45, 55, 1.05]) # azimuth, elevation, distance
    #view_end = np.array([20, 65, 0.90])
    #views = fpl.animate(nframes=nframes, value_beg=view_start, value_end=view_end)

    #move_start = np.array([0.0, -0.015]) # right, up
    #move_end = np.array([-0.01, -0.025])
    #moves = fpl.animate(nframes=nframes, value_beg=move_start, value_end=move_end)

    # loop over files
    for i, filen in enumerate(files):
        # Parse arguments for flashplotlib
        args = fpl.parse_args()
        args.zoom = 1                    # Fixed zoom level
        args.direction = 'z'               # Fixed viewing direction
        args.rot_angle = 0                # Fixed rotation angle (change if needed)
        args.rot_axis = [0, 0, 1]          # Rotation around X-axis
        args.view_angle = 10               # View angle in degrees
        args.outtype = ['pdf']             # Save output as PDF
        args.outname = 'test_{:04d}'.format(i)  # Unique filename per frame
        fpl.process_file(filen, args)

    # time the script
    stop_time = timeit.default_timer()
    total_time = stop_time - start_time
    print("***************** time to finish = "+str(total_time)+"s *****************")
