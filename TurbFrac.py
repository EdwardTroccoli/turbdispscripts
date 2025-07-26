#!/usr/bin/env python
# -*- coding: utf-8 -*-
# written by Edward Troccoli and Christoph Federrath, 2025

from cfpack.defaults import *
from Globals import *
import timeit
import argparse
import cfpack as cfp
import numpy as np
import matplotlib.pyplot as plt
import flashlib as fl
from scipy.stats import linregress

def compute_ekdr_size_fractal_dim(filename, lmax=1.0):
    gg = fl.FlashGG(filename)
    centre = gg.GetMaxLoc("ekdr")
    N = gg.NMax[0]
    D = gg.D[0][0]
    sum = []
    size = []
    for l in range(1, N, 2):
        print("Working on size l = "+str(l*D), color="cyan")
        box_bounds = np.array([centre-l*D/2, centre+l*D/2]).T
        co = gg.GetCells("ekdr", box_bounds=box_bounds)
        sum.append(co.cell_dat.sum())
        size.append(l*D)
    sum = np.log(np.array(sum))
    size = np.log(np.array(size))
    slope, intercept, r_value, p_value, std_err = linregress(size, sum)
    stop()
    return slope

if __name__ == "__main__":

    # Argument parser setup
    parser = argparse.ArgumentParser(description="Compute fractal dimension from FLASH simulation data.")
    parser.add_argument("-ov", "--overwrite", action='store_true', default=False, help="Overwrite files")

    args = parser.parse_args()

    # Start timing the process
    start_time = timeit.default_timer()

    # loop over simulations
    for i, path in enumerate(sim_paths):

        print(f'\nWorking on: {path}', color='cyan')

        # creates the file output dir
        out_path = path + "FracDim/"
        if not os.path.isdir(out_path):
            cfp.run_shell_command('mkdir '+out_path)

        frac_dims = []
        for d in range(20, 101):
            filename = "Turb_hdf5_plt_cnt_{:04d}".format(d)
            d = compute_ekdr_size_fractal_dim(path+filename)
            frac_dims.append(d)

        print(f"Estimated fractal dimension: {np.mean(frac_dims):.3f}")
