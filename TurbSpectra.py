#!/usr/bin/env python
# -*- coding: utf-8 -*-
# written by Edward Troccoli and Christoph Federrath, 2025

import argparse
import numpy as np
import timeit
import os, time
import cfpack as cfp
from cfpack.defaults import *
from turblib import aver_spect, write_spect, read_spect
from Globals import *
import glob
cfp.import_matplotlibrc(fontscale=0.8)

# plotting function for spectras
def plot_spectra(dat, var, Mach,save=False):
    xlabel, ylabel = None, None
    if 0.2 == Mach: MachNum = '0p2'
    if 5 == Mach: MachNum = '5'
    if var == "vels": 
        if 0.2 == Mach: ylabel=r'Power spectrum of $\mathcal{M}$'
    if var == "ekdr": 
        xlabel = r'Wavenumber $k$'
        if  0.2 == Mach: ylabel=r'Power spectrum of $\varepsilon_\mathrm{kin}$'
    if var == "sqrtrho": 
        if  0.2 == Mach: ylabel=r'Power spectrum of $E_\mathrm{kin}$'
    y = 10**dat['col6']
    sigylo = y - 10**(dat['col6']-dat['col7'])
    sigyup = 10**(dat['col6']+dat['col7']) - y
    ret = cfp.plot(
    x=dat['col1'],
    y=y,
    yerr=[sigylo, sigyup],
    shaded_err=True)
    ax = ret.ax()
    ax.text(0.75, 0.95, rf"$\mathcal{{M}} = {Mach}$", transform=ax.transAxes,
        fontsize=14, color='black', verticalalignment='top',
        bbox=dict(boxstyle="round,pad=0.3", facecolor='gray', alpha=0))
    cfp.plot(ax=ret.ax(), xlabel=xlabel, ylabel=ylabel, xlog=True, ylog=True,
            save=False)

def make_paper_plots():
    # loop over figures
    vars = ["sqrtrho", "ekdr", "vels"]
    for ivar, var in enumerate(vars):
        # loop over Mach numbers
        machs = [0.2, 5]
        for imach, mach in enumerate(machs):
            if mach == 0.2:
                sims = ["N2048M0p2HDRe2500HP", "N1024M0p2HDRe2500", "N512M0p2HDRe2500", "N256M0p2HDRe2500"]
                MachNum = '0p2'
                MachSim = 'Sub'
                compensation = '5/3'
            if mach == 5:
                sims = ["N2048M5HDRe2500HP", "N1024M5HDRe2500", "N512M5HDRe2500", "N256M5HDRe2500"]
                MachNum = '5'
                MachSim = 'Sup'
                compensation = '2'   
            color = ['black', 'magenta', 'green', 'grey']
            linestyle = ['solid', 'dashed', 'dashdot', 'dotted']
            dx = [0, 0.24, 0.24, 0.235]
            # loop over simulations
            for isim, sim in enumerate(sims):
                # get sim parameters
                N = params(sim).N
                Mach = params(sim).Mach
                t_turb = params(sim).t_turb
                ylabel = ''
                if var == 'vels': 
                    spectra_files = sorted(glob.glob("../"+sim+"/spectra/"+"*_spect_vels.dat"))
                    norm = 1
                    ylabel = r'$P_{\mathcal{M}}(k)\,k^{'+compensation+'}$'
                if var == 'ekdr': 
                    spectra_files = sorted(glob.glob("../"+sim+"/spectra/"+"*_spect_dset_ekdr.dat"))
                    norm = (t_turb/Mach**2)**2
                    if Mach == 0.2:
                        ylabel = r'$P_{\varepsilon_\mathrm{kin}}(k)\,/\,(\mathcal{M}^2t_\mathrm{turb}^{-1})^2$'
                if var == 'sqrtrho': 
                    spectra_files = sorted(glob.glob("../"+sim+"/spectra/"+"*_spect_sqrtrho.dat"))
                    norm = 1.0/(Mach*N)
                    ylabel = r'$P_{\rho^{1/2}v}(k)\,k^{'+compensation+'}$'
                # read data
                dat = "../"+sim+"/spectra/aver_spectra_"+var+"_M"+MachNum+".dat"
                if not os.path.isfile(dat):
                    aver_dat, header_aver = aver_spect(spectra_files) # average the spectras
                    write_spect(dat, aver_dat, header_aver) # write averaged spectra
                spectra_dat, spectra_header = read_spect(dat) 
                # plot
                xpos, ypos, length = 0.012, 0.085, 1.4
                lf = cfp.legend_formatter(pos=(xpos+isim*dx[isim], ypos), length=length)
                ylim = None
                if var in ['vels','sqrtrho']: 
                    compensation_factor = spectra_dat['col1']**eval(compensation)
                    xlabel = None
                    axes_format=['', None]
                    legend_formatter = None
                    label = None
                    if Mach == 0.2: 
                        ylim = [1e-10,1e-1]
                if var == 'ekdr': 
                    axes_format=[None, None]
                    compensation_factor = 1
                    xlabel = r'Wavenumber $k$'
                    legend_formatter = lf
                    label=MachSim+str(N)
                    ylim = [5e-5,5e0]
                y = 10**spectra_dat['col6']*norm*compensation_factor
                sigylo = y - 10**(spectra_dat['col6']-spectra_dat['col7'])*norm*compensation_factor
                sigyup = 10**(spectra_dat['col6']+spectra_dat['col7'])*norm*compensation_factor - y
                ret = cfp.plot(x=spectra_dat['col1'], y=y, label=label,
                                color=color[isim], linestyle=linestyle[isim], legend_formatter=legend_formatter,
                                yerr=[sigylo, sigyup], shaded_err=[color[isim],0.2])                
            # add Mach label
            #cfp.plot(x=0.75, y=0.95, text=rf"$\mathcal{{M}} = {mach}$", normalised_coords=True)
            #create final plot
            cfp.plot(axes_format=axes_format, xlabel=xlabel, ylabel=ylabel, xlog=True, ylim=ylim, ylog=True, save=fig_path+"spectra_"+var+"_M"+MachNum+".pdf")


if __name__ == "__main__":

    # Argument parser setup
    parser = argparse.ArgumentParser(description="Create and plot PDFs of different variables from FLASH simulation data.")
    var_choices = ["sqrtrho", "vels", "ekdr"]
    parser.add_argument("-v", "--variable", nargs='*', choices=var_choices, help="Variable to plot; choice of "+str(var_choices))
    parser.add_argument("-paper", "--paper_plots", action='store_true', default=False, help="Runs all movie frame plots at paper level quality")
    args = parser.parse_args()

    # Start timing the process
    start_time = timeit.default_timer()

    # create paper plots
    if args.paper_plots:
        make_paper_plots()

    if args.variable is not None:
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
            for var in args.variable:
                spectra_aver_file = out_path+"aver_spectra_"+var+"_M"+MachNum+".dat"
                if not os.path.isfile(spectra_aver_file):
                    if var == 'vels': spectra_files = sorted(glob.glob(out_path+"*_spect_vels.dat"))
                    if var == 'ekdr': spectra_files = sorted(glob.glob(out_path+"*_spect_dset_ekdr.dat"))
                    if var == 'sqrtrho': spectra_files = sorted(glob.glob(out_path+"*_spect_sqrtrho.dat"))
                    aver_dat, header_aver = aver_spect(spectra_files) # average the spectras
                    write_spect(spectra_aver_file, aver_dat, header_aver) # write averaged spectra

                # plot the spectra
                spectra_dat, spectra_header = read_spect(spectra_aver_file) # read the PDF data
                plot_spectra(spectra_dat, var, Mach,save=fig_path+'aver_spectra'+ "_" + var + "_" + "M" +MachNum +'.pdf')

    # End timing and output the total processing time
    stop_time = timeit.default_timer()
    print("Processing time: {:.2f}s".format(stop_time - start_time))
