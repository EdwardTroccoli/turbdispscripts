#!/usr/bin/env python
# -*- coding: utf-8 -*-
# written by Edward Troccoli & Christoph Federrath, 2025

import argparse
import numpy as np
import timeit
import os
import cfpack as cfp
from Globals import *
from matplotlib import rcParams
cfp.load_plot_style() # to load cfpack plot style

def plot_var(path, dat, variable):
    t_turb = params(path).t_turb
    Mach = params(path).Mach
    time = dat['01_time'] / t_turb
    xlabel=r'$t/t_{\mathrm{turb}}$'
    if variable == "mach":
        var = dat['#14_rms_velocity']
        ylabel = r'$\mathcal{M}$'
    elif variable == "ekin":
        var = dat['#10_E_kinetic']
        ylabel = r"Kinetic energy $E_{\mathrm{kin}}/(\langle\rho\rangle\,\mathcal{M}^2\, c_{\mathrm{s}}^2)$"
    elif variable == "injr":
        var = dat['#41_injection_rate']
        ylabel = r'injection rate'
    elif variable == "ekdr":
        var = dat['#42_ekin_diss_rate']
        ylabel = r'$\varepsilon_{\mathrm{kin}}$'
    # create plot
    ret = cfp.plot(x=time, y=var)
    # add Mach label
    ax = ret.ax()
    ax.text(0.05, 0.95, rf"$\mathcal{{M}} = {Mach}$", transform=ax.transAxes, fontsize=14, color='black',
            verticalalignment='top', bbox=dict(boxstyle="round,pad=0.3", facecolor='gray', alpha=0))
    # create final plot
    cfp.plot(xlabel=xlabel, ylabel=ylabel, save=out_path+"tevol_"+variable+".pdf", legend_loc='upper left')


def make_paper_plots():
    # loop over figures
    figs = ["Time_evol_correlation", "Time_evol_injr_ekdr", "Time_evol_Mach"]
    for fig in figs:
        # loop over Mach numbers
        machs = [0.2, 5]
        for mach in machs:
            if mach == 0.2:
                sims = ["N2048M0p2HDRe2500HP", "N1024M0p2HDRe2500", "N512M0p2HDRe2500", "N256M0p2HDRe2500"]
                MachStr = '0p2'
                MachSim = 'Sub'
            if mach == 5:
                sims = ["N2048M5HDRe2500HP", "N1024M5HDRe2500", "N512M5HDRe2500", "N256M5HDRe2500"]
                MachStr = '5'
                MachSim = 'Sup'
            color = ['black', 'magenta', 'green', 'grey']
            linestyle = ['solid', 'dashed', 'dashdot', 'dotted']
            dx = [0,0.225,0.225,0.22]
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
                    if Mach == 5: ylim = [0, 6]
                    ret = cfp.plot(x=time, y=dat['#14_rms_velocity'], color=color[isim], linestyle=linestyle[isim])
                if fig == 'Time_evol_injr_ekdr':
                    injr = dat['#41_injection_rate'] * t_turb / Mach**2
                    ekdr = dat['#42_ekin_diss_rate'] * t_turb / Mach**2
                    xlabel = r'$t/t_\mathrm{turb}$'
                    ylim = [0, 1.4]
                    if mach == 0.2:
                        ylabel = r'$\varepsilon_\mathrm{kin}$ and $\varepsilon_\mathrm{inj}$\quad$[\langle\rho\rangle\,\mathcal{M}^2\,c_{\mathrm{s}}^2\,t_{\mathrm{turb}}^{-1}]$'
                    xpos, ypos, length = 0.072, 0.91, 1.4
                    lf = cfp.legend_formatter(pos=(xpos+isim*dx[isim], ypos), length=length)
                    ret = cfp.plot(x=time, y=ekdr, label=MachSim+str(N), color=color[isim], linestyle=linestyle[isim], legend_formatter=lf)
                    if N==2048:
                        lf = cfp.legend_formatter(pos=(xpos, ypos-0.1), length=length)
                        cfp.plot(x=time, y=injr, label = MachSim+str(N), color='blue', legend_formatter=lf)#label=r'$\varepsilon_{\mathrm{inj}}$'
                        cfp.plot(x=0.022, y=ypos-0.1, text=r"$\varepsilon_{\mathrm{inj}}$:", normalised_coords=True)
                        cfp.plot(x=0.022, y=ypos, text=r"$\varepsilon_{\mathrm{kin}}$:", normalised_coords=True)
                if fig == 'Time_evol_correlation':
                    if N==2048:
                        injr = dat['#41_injection_rate'] * t_turb / Mach**2
                        ekdr = dat['#42_ekin_diss_rate'] * t_turb / Mach**2
                        npts = 1001
                        npts_per_tturb = int(npts / 10)
                        time_int = np.linspace(0, 10, npts) # interpolated time axis
                        ekdr_int = np.interp(time_int, time, ekdr)
                        injr_int = np.interp(time_int, time, injr)

                        time_lags = [] # Overall list of the minima of L2 for each interation
                        for i in range(1,6,4): # Loop such that we look at the region (1,3), (3,5) and so on
                            tshifts = []
                            L2s = []
                            time_int_scan = time_int[i*npts_per_tturb : (i+4)*npts_per_tturb] # Restrict the domains (remove data that was previously used)
                            ekdr_int_scan = ekdr_int[i*npts_per_tturb : (i+4)*npts_per_tturb]
                            injr_int_scan = injr_int[i*npts_per_tturb : (i+4)*npts_per_tturb]

                            # Compute L2 norm
                            for ishift in range(1, 2*npts_per_tturb):
                                tshifts.append(ishift/npts_per_tturb)
                                L2s.append((np.std(injr_int_scan[:-ishift]-ekdr_int_scan[ishift:])/np.mean(injr_int_scan[:-ishift]))**2)
                        
                            # Append results to time_lags and find the minimal L2 in the tshift range.
                            tshifts = np.array(tshifts); L2s=np.array(L2s)
                            time_lags.append(tshifts[L2s==L2s.min()][0])

                        ret = cfp.plot(x=tshifts, y=L2s, color='black')
                        ax = ret.ax()
                        ax.text(0.05, 0.95, rf"$\mathcal{{M}} = {Mach}$", transform=ax.transAxes, fontsize=14, color='black',
                                verticalalignment='top', bbox=dict(boxstyle="round,pad=0.3", facecolor='gray', alpha=0))
                        ylim = [0, 0.35]
                        xlabel = r'$\Delta t/t_{\mathrm{turb}}$'
                        if Mach == 0.2: ylabel = r'$\vert\vert(\varepsilon_{\mathrm{kin}},\varepsilon_{\mathrm{inj}})\vert\vert_{\ell^2}$'
                        # print optimal time shift for max correlation
                        time_lags = np.array(time_lags)
                        print(
                            f"time shift for maximum eps_kin → eps_inj correlation (in t_turb) = {time_lags}\n"
                            f"Result: {time_lags.mean():.3g} ± {time_lags.std():.3g}"
                        )
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
