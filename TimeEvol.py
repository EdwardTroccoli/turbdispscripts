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


def plot_var(variable,path,Mach):

    if not os.path.isdir(fig_path):
        cfp.run_shell_command('mkdir '+fig_path)
    out_path = path + "TimeEvol/"
    if not os.path.isdir(out_path):
        cfp.run_shell_command('mkdir '+out_path)
    dat = cfp.read_ascii(path+"Turb.dat")
    time = dat['01_time'] / t_turb[i]
    xlabel = None
    if variable == "vstd":
        if '0p2' in out_path:
            ylabel = r'$\mathcal{M}$'
            ylim = [0,0.25]
            remove_x_ticks = True
        elif '5' in out_path:
            ylabel = None
            ylim = [0,6.5]
            remove_x_ticks = True
        var = dat['#14_rms_velocity']
        ret = cfp.plot(x=time, y=var, color=color[i])
    elif variable == "ekin":
        var = dat['#10_E_kinetic']*(1/Mach**2)
        ylim = [0,0.65]
        if '0p2' in out_path:
            ylabel = r"Kinetic energy $E_{\textrm{kin}}/(\langle\rho\rangle\,\mathcal{M}^2\, c_{\textrm{s}}^2)$" 
            remove_x_ticks = True
        elif '5' in out_path:
            ylabel = None
            remove_x_ticks = True
        ret = cfp.plot(x=time, y=var, color=color[i])
    elif variable == "emag":
        var = dat['#12_E_magnetic']
        ylabel = r'$E_\mathrm{mag}$'
        cfp.plot(x=time, y=var, label=path[3:-1], color=color[i])
    elif variable == "injr":
        var = dat['#41_injection_rate']
        ylabel = r'injection rate'
        cfp.plot(x=time, y=var, label=path[3:-1], color=color[i])
    elif variable == "ekdr":
        var = dat['#42_ekin_diss_rate']
        ylabel = r'$\varepsilon_{\textrm{kin}}$'
        ret = cfp.plot(x=time, y=var, color=color[i])
    elif variable == "emdr":
        var = dat['#43_emag_diss_rate']
        ylabel = r'magnetic energy dissipation rate'
        cfp.plot(x=time, y=var, label=path[3:-1], color=color[i])
    elif variable == "etdr":
        var = dat['#42_ekin_diss_rate'] + dat['#43_emag_diss_rate']
        ylabel = r'total energy dissipation rate'
        ret = cfp.plot(x=time, y=var, label=Mach, color=color[i])
    elif variable == "ired":
        injr = dat['#41_injection_rate']*(t_turb[i]/Mach**2)
        ekdr = dat['#42_ekin_diss_rate']*(t_turb[i]/Mach**2)
        xlabel = r'$t/t_\mathrm{turb}$'
        ylim=[0,1.5]
        if '0p2' in out_path:
            ylabel = r'$\varepsilon_{\textrm{kin}}$ and $\varepsilon_{\textrm{inj}} / (\langle\rho\rangle\,\mathcal{M}^2\, c_{\textrm{s}}^2\, t_{\textrm{turb}}^{-1})$'
            remove_x_ticks = False
            dat1 = cfp.read_ascii("../N512M0p2HDRe2500/Turb.dat")
            dat2 = cfp.read_ascii("../N256M0p2HDRe2500/Turb.dat")
        elif '5' in out_path:
            dat1 = cfp.read_ascii("../N512M5HDRe2500/Turb.dat")
            dat2 = cfp.read_ascii("../N256M5HDRe2500/Turb.dat")
            ylabel = None
            remove_x_ticks = False
        time1   = dat1['01_time']/t_turb[i]
        time2   = dat2['01_time']/t_turb[i]
        ekdr1   = dat1['#42_ekin_diss_rate'] * (t_turb[i]/Mach**2)
        ekdr2   = dat2['#42_ekin_diss_rate'] * (t_turb[i]/Mach**2)
        if args.interpolate == True:
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
            cfp.plot(x=tshifts, y=L2s, ylim=[0,0.3], xlabel=r'$\Delta t / t_{\textrm{turb}}$',
                     ylabel=r'L2 norm of $\varepsilon_{\textrm{kin}}-\varepsilon_{\textrm{inj}}$',
                     save=out_path+"tevol_"+f"{variable}_M"+MachNumber[i]+"_time_correlation.pdf")
            # get optimal time shift for max correlation
            tshift_max_correlation = tshifts[L2s==L2s.min()]
            print('time shift for maximum eps_kin to eps_inj correlation (in t_turb) = ', tshift_max_correlation)
        cfp.plot(x=time2, y=ekdr2, label=r'$256^3$', color='pink')
        cfp.plot(x=time, y=ekdr, label=r'$1024^3$', color='black')
        cfp.plot(x=time1, y=ekdr1, label=r'$512^3$', color='green')
        #Reorder by desired index
        order = [1, 3, 2, 0]  # example: ε_inj, 1024³, 512³, 256³
        ret = cfp.plot(x=time, y=injr, label=r'$\varepsilon_{\textrm{inj}}$', color='b')
        time1   = dat1['01_time']/t_turb[i]
    ax = ret.ax()
    handles, labels = ax.get_legend_handles_labels()
    ax.legend([handles[i] for i in order], [labels[i] for i in order],bbox_to_anchor=(0.005, 0.93))
    ax.text(0.05, 0.95, rf"$\mathcal{{M}} = {Mach}$", transform=ax.transAxes,
        fontsize=14, color='black', verticalalignment='top',
        bbox=dict(boxstyle="round,pad=0.3", facecolor='gray', alpha=0))
    if remove_x_ticks == True:
        ax.set_xticklabels([])

    cfp.plot(xlabel=xlabel, ylabel=ylabel, ylim=ylim, save=out_path+"tevol_"+f"{variable}_M"+MachNumber[i]+".pdf", legend_loc='upper left')


if __name__ == "__main__":

    # Argument parser setup
    parser = argparse.ArgumentParser(description="Plot different variables from simulation data.")
    var_choices = ["vstd", "ekin", "emag", "injr", "ekdr", "emdr", "etdr", "ired"]
    parser.add_argument("-v", "--variable", nargs='*', choices=var_choices, required=True, help="Variable to plot; choice of "+str(var_choices))
    parser.add_argument("-int", "--interpolate", action='store_true', default=False, help="Allows interpolation of dissipation and injection rates")
    # Parse the command-line arguments
    args = parser.parse_args()

    # Start timing the process
    start_time = timeit.default_timer()

    for i, path in enumerate(sim_paths):

        print(f'Working on: {path}', color='green')
        if '0p2' in MachNumber[i]:
            Mach = 0.2
        else:
            Mach = 5
        for var in args.variable:
            plot_var(var,path,Mach)

    # End timing and output the total processing time
    stop_time = timeit.default_timer()
    print("Processing time: {:.2f}s".format(stop_time - start_time))
