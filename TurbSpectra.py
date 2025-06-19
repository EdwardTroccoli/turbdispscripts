#!/usr/bin/env python
# -*- coding: utf-8 -*-
# written by Edward Troccoli and Christoph Federrath, 2025

import argparse
import timeit
import os, time
import cfpack as cfp
from cfpack.defaults import *
from turblib import aver_spect, write_spect, read_spect
from Globals import *
import glob
cfp.import_matplotlibrc(fontscale=0.8)

# plotting function for spectras
def plot_spectra(dat, var, Mach):
    ylabel = None
    if var == "vels": 
        xlabel = None
        if '0p2' in out_path: ylabel=r'Power spectrum of $E_\mathrm{kin}$'
    if var == "ekdr": 
        xlabel = r'Wavenumber $k$'
        if '0p2' in out_path: ylabel=r'Power spectrum of $\varepsilon_\mathrm{kin}$'
    y = 10**dat['col6']
    sigylo = y - 10**(dat['col6']-dat['col7'])
    sigyup = 10**(dat['col6']+dat['col7']) - y
    ret = cfp.plot(
    x=dat['col1'],
    y=y,
    yerr=[sigylo, sigyup],
    shaded_err=True)
    ax = ret.ax()
    ax.text(0.95, 0.95, rf"$\mathcal{{M}} = {Mach}$", transform=ax.transAxes,
        fontsize=14, color='black', verticalalignment='top',
        bbox=dict(boxstyle="round,pad=0.3", facecolor='gray', alpha=0))
    cfp.plot(ax=ret.ax(), xlabel=xlabel, ylabel=ylabel, xlog=True, ylog=True,
            save=out_path+'aver_spectra'+ "_" + var + "_" + "M" +MachNum +'.pdf')

if __name__ == "__main__":

    # Argument parser setup
    parser = argparse.ArgumentParser(description="Create and plot PDFs of different variables from FLASH simulation data.")
    parser.add_argument("-ov", "--overwrite", action='store_true', default=False, help="Overwrite files")

    args = parser.parse_args()

    # Start timing the process
    start_time = timeit.default_timer()

    # loop over simulations
    for i, path in enumerate(sim_paths):

        print(f'Working on: {path}', color='green')

        # creates the figure output dir
        if not os.path.isdir(fig_path):
            cfp.run_shell_command('mkdir '+fig_path)

        # create output dir for spectra
        out_path = path + "spectra/"
        if not os.path.isdir(out_path):
            cfp.run_shell_command('mkdir '+out_path)

        Mach = params(path).Mach
        if 'M0p2' in out_path: MachNum = '0p2'
        if 'M5' in out_path: MachNum = '5'

        # loop over simulation variables
        vars = ['ekdr', 'vels']
        for var in vars:
            spectra_aver_file = out_path+"aver_spectra_"+var+"_M"+MachNum+".dat"
            if not os.path.isfile(spectra_aver_file) or args.overwrite:
                for d in range(20, 101, 1):
                    filename = "Turb_hdf5_plt_cnt_{:04d}".format(d)
                    compute_spectra(args.overwrite) # compute the spectra by calling C++ 'spectra'
                if var == 'vels':
                    spectra_files = sorted(glob.glob(out_path+"*_spect_vels.dat"))
                elif var == 'ekdr':
                    spectra_files = sorted(glob.glob(out_path+"*_spect_dset_ekdr.dat"))
                aver_dat, header_aver = aver_spect(spectra_files) # average the spectras
                write_spect(spectra_aver_file, aver_dat, header_aver) # write averaged spectra

            # plot the spectra
            spectra_dat, spectra_header = read_spect(spectra_aver_file) # read the PDF data
            plot_spectra(spectra_dat, var, Mach)

    # End timing and output the total processing time
    stop_time = timeit.default_timer()
    print("Processing time: {:.2f}s".format(stop_time - start_time))
