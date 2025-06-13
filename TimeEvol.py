#!/usr/bin/env python
# -*- coding: utf-8 -*-
# written by Edward Troccoli, 2025

import argparse
import numpy as np
import timeit
import os
import cfpack as cfp
from cfpack.defaults import *
from Globals import *


def plot_var(path, dat, variable):
    t_turb = params(path).t_turb
    Mach = params(path).Mach
    color = params(path).color
    ylim = None
    xlabel = None
    remove_x_ticks = False
    time = dat['01_time'] / t_turb
    if variable == "mach":
        ylabel = r'$\mathcal{M}$'
        remove_x_ticks = True
        if Mach == 0.2:
            ylim = [0, 0.25]
        if Mach == 5:
            ylabel = None
            ylim = [0, 6.5]
        var = dat['#14_rms_velocity']
    elif variable == "ekin":
        var = dat['#10_E_kinetic'] / Mach**2
        ylim = [0, 0.65]
        remove_x_ticks = True
        if Mach == 0.2:
            ylabel = r"Kinetic energy $E_{\textrm{kin}}/(\langle\rho\rangle\,\mathcal{M}^2\, c_{\textrm{s}}^2)$"
        if Mach == 5:
            ylabel = None
    elif variable == "injr":
        var = dat['#41_injection_rate']
        ylabel = r'injection rate'
    elif variable == "ekdr":
        var = dat['#42_ekin_diss_rate']
        ylabel = r'$\varepsilon_{\textrm{kin}}$'
    # create plot
    ret = cfp.plot(x=time, y=var, color=color)
    # add Mach label
    ax = ret.ax()
    ax.text(0.05, 0.95, rf"$\mathcal{{M}} = {Mach}$", transform=ax.transAxes, fontsize=14, color='black',
            verticalalignment='top', bbox=dict(boxstyle="round,pad=0.3", facecolor='gray', alpha=0))
    if remove_x_ticks:
        ax.set_xticklabels([])
    # create final plot
    cfp.plot(xlabel=xlabel, ylabel=ylabel, ylim=ylim, save=out_path+"tevol_"+f"{variable}.pdf", legend_loc='upper left')


def make_paper_plots(): # please implement properly
    variable = 'test'
    t_turb = params(path).t_turb
    Mach = params(path).Mach
    color = params(path).color
    time = dat['01_time'] / t_turb
    injr = dat['#41_injection_rate'] * t_turb / Mach**2
    ekdr = dat['#42_ekin_diss_rate'] * t_turb / Mach**2
    xlabel = r'$t/t_\mathrm{turb}$'
    ylim = [0, 1.5]
    if Mach == 0.2:
        ylabel = r'$\varepsilon_{\textrm{kin}}$ and $\varepsilon_{\textrm{inj}} / (\langle\rho\rangle\,\mathcal{M}^2\, c_{\textrm{s}}^2\, t_{\textrm{turb}}^{-1})$'
        remove_x_ticks = False
        dat1 = cfp.read_ascii("../N512M0p2HDRe2500/Turb.dat_cleaned")
        dat2 = cfp.read_ascii("../N256M0p2HDRe2500/Turb.dat_cleaned")
    if Mach == 5:
        dat1 = cfp.read_ascii("../N512M5HDRe2500/Turb.dat_cleaned")
        dat2 = cfp.read_ascii("../N256M5HDRe2500/Turb.dat_cleaned")
        ylabel = None
        remove_x_ticks = False
    time1 = dat1['01_time'] / t_turb
    time2 = dat2['01_time'] / t_turb
    ekdr1 = dat1['#42_ekin_diss_rate'] * t_turb / Mach**2
    ekdr2 = dat2['#42_ekin_diss_rate'] * t_turb / Mach**2
    if args.interpolate:
        npts = 1001
        npts_per_tturb = int(npts / 10)
        time_int = np.linspace(0, 10, npts) # interpolated time axis
        ekdr_int = np.interp(time_int, time, ekdr)
        injr_int = np.interp(time_int, time, injr)
        tshifts = []
        L2s = []
        for ishift in range(1, 3*npts_per_tturb):
            tshifts.append(ishift/npts_per_tturb)
            L2s.append((np.std(injr_int[:-ishift]-ekdr_int[ishift:])/np.mean(injr_int[:-ishift]))**2)
        tshifts = np.array(tshifts); L2s=np.array(L2s)
        ret = cfp.plot(x=tshifts, y=L2s)
        ax = ret.ax()
        ax.text(0.05, 0.95, rf"$\mathcal{{M}} = {Mach}$", transform=ax.transAxes,
            fontsize=14, color='black', verticalalignment='top',
            bbox=dict(boxstyle="round,pad=0.3", facecolor='gray', alpha=0))
        cfp.plot(ylim=[0,0.35], xlabel=r'$t / t_{\textrm{turb}}$',
                    ylabel=r'$\ell^2$ norm of $\varepsilon_{\textrm{kin}}-\varepsilon_{\textrm{inj}}$',
                    save=out_path+"tevol_"+f"{variable}_time_correlation.pdf")
        # get optimal time shift for max correlation
        tshift_max_correlation = tshifts[L2s==L2s.min()]
        print('time shift for maximum eps_kin to eps_inj correlation (in t_turb) = ', tshift_max_correlation,color = 'green')
    cfp.plot(x=time2, y=ekdr2, label=r'$256^3$', color='pink')
    cfp.plot(x=time, y=ekdr, label=r'$1024^3$', color='black')
    cfp.plot(x=time1, y=ekdr1, label=r'$512^3$', color='green')
    # reorder by desired index
    order = [0, 2, 1, 3] 
    ret = cfp.plot(x=time, y=injr, label=r'$\varepsilon_{\textrm{inj}}$', color='b')
    time1 = dat1['01_time'] / t_turb[i]
    ax = ret.ax()
    handles, labels = ax.get_legend_handles_labels()
    ax.legend([handles[i] for i in order], [labels[i] for i in order], bbox_to_anchor=(0.005, 0.93))
    # add Mach label
    ax.text(0.05, 0.95, rf"$\mathcal{{M}} = {Mach}$", transform=ax.transAxes, fontsize=14, color='black',
            verticalalignment='top', bbox=dict(boxstyle="round,pad=0.3", facecolor='gray', alpha=0))
    if remove_x_ticks:
        ax.set_xticklabels([])
    # create final plot
    cfp.plot(xlabel=xlabel, ylabel=ylabel, ylim=ylim, save=fig_path+"tevol_"+f"{variable}.pdf", legend_loc='upper left')



if __name__ == "__main__":

    # Argument parser setup
    parser = argparse.ArgumentParser(description="Plot different time evolution variables from simulation data.")
    var_choices = ["mach", "ekin", "injr", "ekdr"]
    parser.add_argument("-v", "--variable", nargs='*', choices=var_choices, required=True, help="Variable to plot; choice of "+str(var_choices))
    parser.add_argument("-ov", "--overwrite", action='store_true', default=False, help="Overwrite files")
    args = parser.parse_args()

    # loop over simulation paths (in Global.py)
    for i, path in enumerate(sim_paths):
        print(f'Working on: {path}', color='green')
        # clean the .dat file
        datfile = path+'Turb.dat'
        if not os.path.isfile(datfile+'_cleaned') or args.overwrite:
            cfp.run_shell_command('flashlib.py datfile -clean '+datfile)
        out_path = path + "TimeEvol/"
        if not os.path.isdir(out_path):
            cfp.run_shell_command('mkdir '+out_path)
        # read cleaned time evolution file
        dat = cfp.read_ascii(datfile+'_cleaned')
        # loop over variables to plot
        for var in args.variable:
            plot_var(path, dat, var)
