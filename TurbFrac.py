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


def extract_ekdr(filename, variable='ekdr', overwrite=False, norm=None):
    fname_pkl = f"{out_path+os.path.basename(filename)}_{variable}_frac_dim_M{MachStr}.pkl"
    if not os.path.isfile(fname_pkl) or overwrite:
        # read data
        gg = fl.FlashGG(filename)
        print(f"Reading {variable} data from plot files...", color="red")
        field = gg.GetUniformGrid(dset=variable)
        max_index_flat = np.argmax(field)
        max_loc = np.unravel_index(max_index_flat, field.shape)
        ret = {
            'field': field,
            'max_loc': max_loc,
        }
        with open(fname_pkl, "wb") as fobj:
            print(f"Writing '{fname_pkl}'", color="magenta")
            dill.dump(ret, fobj)
    else:
        print("Read '"+fname_pkl+"'", color="green")
        ret = dill.load(open(fname_pkl, "rb"))
    return ret


def compute_mass_radius_fractal_dim(field, center, r_max=10):
    shape = field.shape
    radii = np.arange(1, r_max + 1)
    masses = []

    for r in radii:
        xmin = max(0, center[0] - r)
        xmax = min(shape[0], center[0] + r + 1)
        ymin = max(0, center[1] - r)
        ymax = min(shape[1], center[1] + r + 1)
        zmin = max(0, center[2] - r)
        zmax = min(shape[2], center[2] + r + 1)

        xg = np.arange(xmin, xmax)
        yg = np.arange(ymin, ymax)
        zg = np.arange(zmin, zmax)
        xx, yy, zz = np.meshgrid(xg, yg, zg, indexing='ij')

        dist = np.sqrt((xx - center[0])**2 + (yy - center[1])**2 + (zz - center[2])**2)
        filter = dist <= r

        subfield = field[xmin:xmax, ymin:ymax, zmin:zmax]
        mass = subfield[filter].sum()
        masses.append(mass)

    log_r = np.log(radii)
    log_m = np.log(masses)
    slope, intercept = np.polyfit(log_r, log_m, 1)
    D = slope

    return D, radii, masses



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
    sum = np.array(sum)
    size = np.array(size)
    stop()


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
            compute_ekdr_size_fractal_dim(path+filename)
            #D, radii, masses = compute_mass_radius_fractal_dim(field, center, r_max=10)
            #frac_dims.append(D)

        print(f"Estimated fractal dimension: {np.mean(frac_dims):.3f}")
