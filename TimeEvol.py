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
cfp.import_matplotlibrc(fontscale=0.8)
from matplotlib import rcParams


def plot_var(path, dat, variable):
    t_turb = params(path).t_turb
    Mach = params(path).Mach
    time = dat['01_time'] / t_turb
    xlabel=r'$t/t_{\textrm{turb}}$'
    if variable == "mach":
        var = dat['#14_rms_velocity']
        ylabel = r'$\mathcal{M}$'
    elif variable == "ekin":
        var = dat['#10_E_kinetic']
        ylabel = r"Kinetic energy $E_{\textrm{kin}}/(\langle\rho\rangle\,\mathcal{M}^2\, c_{\textrm{s}}^2)$"
    elif variable == "injr":
        var = dat['#41_injection_rate']
        ylabel = r'injection rate'
    elif variable == "ekdr":
        var = dat['#42_ekin_diss_rate']
        ylabel = r'$\varepsilon_{\textrm{kin}}$'
    # create plot
    ret = cfp.plot(x=time, y=var)
    # add Mach label
    ax = ret.ax()
    ax.text(0.05, 0.95, rf"$\mathcal{{M}} = {Mach}$", transform=ax.transAxes, fontsize=14, color='black',
            verticalalignment='top', bbox=dict(boxstyle="round,pad=0.3", facecolor='gray', alpha=0))
    # create final plot
    cfp.plot(xlabel=xlabel, ylabel=ylabel, save=out_path+"tevol_"+variable+".pdf", legend_loc='upper left')


def make_paper_plots(): # please implement properly
    # loop over figures
    figs = ["Time_evol_Mach", "Time_evol_injr_ekdr", "Time_correlation"]
    for fig in figs:
        # loop over Mach numbers
        machs = [0.2, 5]
        for mach in machs:
            if mach == 0.2:
                sims = ["N256M0p2HDRe2500", "N512M0p2HDRe2500", "N1024M0p2HDRe2500", "N2048M0p2HDRe2500HP"]
                MachStr = '0p2'
            if mach == 5:
                sims = ["N256M5HDRe2500", "N512M5HDRe2500", "N1024M5HDRe2500", "N2048M5HDRe2500HP"]
                MachStr = '5'
            color = ['grey', 'green', 'magenta', 'black']
            linestyle = ['dotted', 'dashdot', 'dashed', 'solid']
            # loop over simulations
            for isim, sim in enumerate(sims):
                # get sim parameters
                N = params(sim).N
                Mach = params(sim).Mach
                t_turb = params(sim).t_turb
                # read data
                clean_datfile("../"+sim+'/')
                dat = cfp.read_ascii("../"+sim+"/Turb.dat_cleaned")
                time = dat['01_time'] / t_turb
                # plot
                remove_x_ticks = False
                xlabel, ylabel = None, None
                if fig == 'Time_evol_Mach':
                    ylabel = r'$\mathcal{M}$'
                    remove_x_ticks = True
                    if Mach == 0.2: ylim = [0, 0.25]
                    if Mach == 5: ylim = [0, 6.5]
                    ret = cfp.plot(x=time, y=dat['#14_rms_velocity'], color=color[isim], linestyle=linestyle[isim])
                if fig == 'Time_evol_injr_ekdr':
                    injr = dat['#41_injection_rate'] * t_turb / Mach**2
                    ekdr = dat['#42_ekin_diss_rate'] * t_turb / Mach**2
                    xlabel = r'$t/t_\mathrm{turb}$'
                    ylim = [0, 1.5]
                    if mach == 0.2:
                        ylabel = r'$\varepsilon_{\textrm{kin}}$ and $\varepsilon_{\textrm{inj}} / (\langle\rho\rangle\,\mathcal{M}^2\, c_{\textrm{s}}^2\, t_{\textrm{turb}}^{-1})$'
                    xpos, ypos, dx, length = 0.012, 0.84, 0.19, 1.4
                    lf = cfp.legend_formatter(pos=(xpos+isim*dx, ypos), length=length)
                    ret = cfp.plot(x=time, y=ekdr, label=r'$'+str(N)+'^3$', color=color[isim], linestyle=linestyle[isim], legend_formatter=lf)
                    if N==2048:
                        lf = cfp.legend_formatter(pos=(xpos+4*dx, ypos), length=length)
                        cfp.plot(x=time, y=injr, label=r'$\varepsilon_{\textrm{inj}}$', color='blue', legend_formatter=lf)
                if fig == 'Time_correlation':
                    if N==2048:
                        injr = dat['#41_injection_rate'] * t_turb / Mach**2
                        ekdr = dat['#42_ekin_diss_rate'] * t_turb / Mach**2
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
                        ret = cfp.plot(x=tshifts, y=L2s, color='black')
                        ylim = [0, 0.35]
                        xlabel = r'$\Delta t/t_{\textrm{turb}}$'
                        ylabel = r'$\ell^2$ norm of $\varepsilon_{\textrm{kin}}-\varepsilon_{\textrm{inj}}$'
                        # print optimal time shift for max correlation
                        tshift_max_correlation = tshifts[L2s==L2s.min()]
                        print('time shift for maximum eps_kin to eps_inj correlation (in t_turb) = ', tshift_max_correlation, color='yellow')
            # add Mach label
            cfp.plot(x=0.04, y=0.91, text=rf"$\mathcal{{M}} = {mach}$", normalised_coords=True)
            if remove_x_ticks:
                ax = ret.ax()
                ax.set_xticklabels([])
            # create final plot
            cfp.plot(xlabel=xlabel, ylabel=ylabel, ylim=ylim, save=fig_path+fig+"_M"+MachStr+".pdf")


if __name__ == "__main__":

    # Argument parser setup
    parser = argparse.ArgumentParser(description="Plot different time evolution variables from simulation data.")
    var_choices = ["mach", "ekin", "injr", "ekdr"]
    parser.add_argument("-v", "--variable", nargs='*', choices=var_choices, help="Variable to plot; choice of "+str(var_choices))
    parser.add_argument("-ov", "--overwrite", action='store_true', default=False, help="Overwrite files")
    parser.add_argument("-paper", "--paper_plots", action='store_true', default=False, help="Runs all movie frame plots at paper level quality")
    args = parser.parse_args()

    # create paper plots
    if args.paper_plots:
        make_paper_plots()

    # flexible plotting of variable in selected sims based on Globals
    if args.variable is not None:
        # loop over simulation paths (in Global.py)
        for path in sim_paths:
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
