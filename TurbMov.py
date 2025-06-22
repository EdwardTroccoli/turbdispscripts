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
cfp.import_matplotlibrc(fontscale=0.8)

# function that is called to create the plot_maps of the various variables, saves having to call the same script over and over
def plot_frame(var, log, cmap_label, cmap, vmin, vmax, remove_x_ticks, remove_y_ticks, xlabel, ylabel, out_file):
    ret = cfp.plot_map(var, log=log, cmap_label=cmap_label, cmap=cmap, xlim=[0,1], ylim=[0,1], aspect_data='equal', vmin=vmin, vmax=vmax)
    ax = ret.ax()[0]
    # create box for Mach number label
    ax.text(0.05, 0.95, r"$\mathcal{M}=$ "+str(Mach), transform=ax.transAxes, color='white', verticalalignment='top',
            bbox=dict(boxstyle="round,pad=0.3", facecolor='gray', alpha=0.5))
    if remove_x_ticks == True:
        ax.set_xticklabels([])
    if remove_y_ticks == True:
        ax.set_yticklabels([])
    #time = hdfio.read(filen, "time")[0] / tturb
    #time_str = cfp.round(time, 3, str_ret=True)
    cfp.plot(ax=ret.ax()[0], xlabel=xlabel, ylabel=ylabel, color='white', normalised_coords=True, save=out_file)#text=r"$t/t_\mathrm{turb}="+time_str+r"$",

def plot_variable(slices, plot_files, out_path, variable, t_turb, Mach, N):

    #set defaults for quick visualisation plots
    log = True
    cmap = 'afmhot'
    cmap_label = None
    xlabel, ylabel = r'$x$', r'$y$'
    vmin, vmax = None, None
    remove_x_ticks, remove_y_ticks = None, None

    if '0p2' in out_path:
        MachNum = '0p2'
    elif '5' in out_path:
        MachNum = '5'
    #need to handle vorticity and other variables separately as they come from plot files and slice files respectively.

    #loop over plot files for vorticity
    if variable == "vort":
        cmap_label = r"Vorticity $|\nabla\times\mathbf{v}|/(\mathcal{M}c_{\textrm{s}}\Delta x^{-1})$"
        for i, filen in enumerate(plot_files):
                #call vorticity computer from Globals to make the slice files
                compute_vort(filen, overwrite=False)
                #need to construct the filenames appropiately as they are stored in the movie_files directory
                vortx = os.path.dirname(filen)+'/movie_files/'+os.path.basename(filen+'_vorticity_x_slice_z.h5')
                vorty = os.path.dirname(filen)+'/movie_files/'+os.path.basename(filen+'_vorticity_y_slice_z.h5')
                vortz = os.path.dirname(filen)+'/movie_files/'+os.path.basename(filen+'_vorticity_z_slice_z.h5')
                vortx,vorty,vortz = hdfio.read(vortx, "vorticity_x_slice"), hdfio.read(vorty, "vorticity_y_slice"), hdfio.read(vortz, "vorticity_z_slice")

                var = np.sqrt(vortx**2+vorty**2+vortz**2)*((1/N)/ Mach)
                out_file = out_path+f"frame_{variable}_{i:06d}.png"

                plot_frame(var, log, cmap_label, cmap, vmin, vmax, remove_x_ticks, remove_y_ticks, xlabel, ylabel, out_file)
        return

    #loop over movie files for other variables
    for i, filen in enumerate(slices):

        if variable == "dens":
            cmap_label = r'Density $\rho/\langle\rho\rangle$'
            var = hdfio.read(filen, "dens_slice_xy")
        elif variable == "ekin":
            cmap = 'RdPu'
            dens = hdfio.read(filen, "dens_slice_xy")
            velx,vely,velz = hdfio.read(filen, "velx_slice_xy"), hdfio.read(filen, "vely_slice_xy"), hdfio.read(filen, "velz_slice_xy")
            var = 0.5 * dens * (velx ** 2 + vely ** 2 + velz ** 2)*(1/Mach**2)
            cmap_label = r"Kinetic energy $E_{\textrm{kin}}/(\langle\rho\rangle\,\mathcal{M}^2\, c_{\textrm{s}}^2)$"
        elif variable == "ekdr":
            cmap = 'BuPu'
            var = hdfio.read(filen, "ekdr_slice_xy")*(t_turb/Mach**2)
            cmap_label = r"Dissipation rate $\varepsilon_{\textrm{kin}}/(\langle\rho\rangle\,\mathcal{M}^2\, c_{\textrm{s}}^2\,t_{\textrm{turb}}^{-1}$)"

        # Define formatted filename correctly
        #out_file = out_path+f"frame_{variable}_{i:06d}.png"
        out_file = out_path+f"frame_{variable}_000250_M{MachNum}.pdf"
        plot_frame(var, log, cmap_label, cmap, vmin, vmax, remove_x_ticks, remove_y_ticks, xlabel, ylabel, out_file)

def make_paper_plots(slices, plot_files, out_path, t_turb, Mach, N):

    #variables that are relevant for the paper to plot
    variables = ['dens', 'ekdr', 'vort']

    #set default parameters
    log = True
    cmap = 'afmhot'
    ekdr_min, ekdr_max = 1e-4, 1e3

    #same logic for handling the paper plots
    for variable in variables:

        vmin,vmax = None,None
        remove_x_ticks,remove_y_ticks = None,None
        cmap_label = None
        xlabel, ylabel = None, None

        if 'M0p2' in out_path: MachNum = '0p2'
        if 'M5' in out_path: MachNum = '5'

        if variable == "vort":
            cmap_label = r"Vorticity $|\nabla\times\mathbf{v}|/(\mathcal{M}c_{\textrm{s}}\Delta x^{-1})$"
            for i, filen in enumerate(plot_files):
                    compute_vort(filen,overwrite=False)
                    vortx = os.path.dirname(filen)+'/movie_files/'+os.path.basename(filen+'_vorticity_x_slice_z.h5')
                    vorty = os.path.dirname(filen)+'/movie_files/'+os.path.basename(filen+'_vorticity_y_slice_z.h5')
                    vortz = os.path.dirname(filen)+'/movie_files/'+os.path.basename(filen+'_vorticity_z_slice_z.h5')
                    vortx,vorty,vortz = hdfio.read(vortx, "vorticity_x_slice"), hdfio.read(vorty, "vorticity_y_slice"), hdfio.read(vortz, "vorticity_z_slice")

                    var = np.sqrt(vortx**2+vorty**2+vortz**2)*((1/N) / Mach)
                    out_file = out_path+f"frame_{variable}_{i:06d}.png"

                    plot_frame(var, log, cmap_label, cmap, vmin, vmax, remove_x_ticks, remove_y_ticks, xlabel, ylabel, out_file)
            return
        for i, filen in enumerate(slices):
            if variable == "dens":
                if '0p2' in out_path:
                    vmin = 0.95
                    vmax = 1.05
                    ylabel = r'$y$'
                elif '5' in out_path:
                    vmin = 1e-2
                    vmax = 1e2
                    cmap_label = r'Density $\rho/\langle\rho\rangle$'
                    remove_y_ticks = True
                var = hdfio.read(filen, "dens_slice_xy")
                xlabel = r'$x$'
            elif variable == "ekdr":
                remove_x_ticks = True
                cmap = 'BuPu'
                vmin,vmax = ekdr_min,ekdr_max
                var = hdfio.read(filen, "ekdr_slice_xy")*(t_turb/Mach**2)
                if '0p2' in out_path:
                    ylabel = r'$y$'
                elif '5' in out_path:
                    cmap_label = r"Dissipation rate $\varepsilon_{\textrm{kin}}/(\langle\rho\rangle\,\mathcal{M}^2\, c_{\textrm{s}}^2\,t_{\textrm{turb}}^{-1}$)"
                    remove_y_ticks = True

        out_file = out_path+f"frame_{variable}_000250_M{MachNum}.pdf"
        plot_frame(var, log, cmap_label, cmap, vmin, vmax, remove_x_ticks, remove_y_ticks, xlabel, ylabel, out_file)


if __name__ == "__main__":

    # Argument parser setup
    parser = argparse.ArgumentParser(description="Plot different variables from simulation data.")
    var_choices = ["dens", "ekin", "ekdr", "vort"]
    parser.add_argument("-v", "--variable", choices=var_choices, nargs= '*', help="Variable to plot; choice of "+str(var_choices))
    parser.add_argument("-paper", "--paper_plots", action='store_true', default=False, help="Runs all movie frame plots at paper level quality")
    # Parse the command-line arguments
    args = parser.parse_args()

    # Start timing the process
    start_time = timeit.default_timer()

    # loop over simulations
    for i, path in enumerate(sim_paths):

        print(f'Working on: {path}', color='green')
        
        t_turb = params(path).t_turb
        Mach = params(path).Mach
        N = params(path).N

        out_path = path + "movie_files/"
        if not os.path.isdir(out_path):
            cfp.run_shell_command('mkdir '+out_path)

        # Get all files matching the pattern Turb_slice_xy_*
        #files = sorted(glob.glob(out_path+"Turb_slice_xy_*"))
        slices = sorted(glob.glob(out_path+"Turb_slice_xy_000250"))
        plot_files = sorted(glob.glob(path+"Turb_hdf5_plt_cnt_0050"))

        if args.paper_plots == False:
            for var in args.variable:
                # call plot function
                plot_variable(slices, plot_files, out_path, var, t_turb, Mach, N)
        else:
            make_paper_plots(slices, plot_files, out_path, t_turb, Mach, N)

    # End timing and output the total processing time
    stop_time = timeit.default_timer()
    print("Processing time: {:.2f}s".format(stop_time - start_time))