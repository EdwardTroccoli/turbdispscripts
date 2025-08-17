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


def make_paper_plots():
        machs = [0.2, 5]
        for imach, mach in enumerate(machs):
            if mach == 0.2:
                sims = ["N1024M0p2HDRe2500", "N512M0p2HDRe2500", "N256M0p2HDRe2500"]
                MachNum = '0p2'
                MachSim = 'Sub'
            if mach == 5:
                sims = ["N1024M5HDRe2500", "N512M5HDRe2500", "N256M5HDRe2500"]
                MachNum = '5'
                MachSim = 'Sup' 
            color = ['black', 'magenta', 'green', 'grey']
            linestyle = ['solid', 'dashed', 'dashdot', 'dotted']
            dy = [0.1, 0.1, 0.1]
            # loop over simulations
            for isim, sim in enumerate(sims):
                # get sim parameters
                N = params(sim).N
                Mach = params(sim).Mach
                t_turb = params(sim).t_turb
                if Mach == 0.2: MachStr = '0p2'
                if Mach == 5:   MachStr = '5'
                # read data
                ylabel = r"Dissipation rate $\varepsilon_{\textrm{kin}}/(\langle\rho\rangle\,\mathcal{M}^2\, c_{\textrm{s}}^2\,t_{\textrm{turb}}^{-1}$)"
                bsdat = get_ekdr_size_fractal_dim("../"+sim+"/")

                # plot
                xpos, ypos, length = 0.7, 0.1, 1.4
                lf = cfp.legend_formatter(pos=(xpos, ypos+isim*dy[isim]), length=length)
                ret = cfp.plot(x=bsdat.x, y=bsdat.y, label=MachSim+str(N),
                                color=color[isim], linestyle=linestyle[isim], legend_formatter=lf
                                )                
            # add Mach label
            #cfp.plot(x=0.75, y=0.95, text=rf"$\mathcal{{M}} = {mach}$", normalised_coords=True)
            #create final plot
            for i in range(1,4):
                xpos, ypos, length, dy = 0.55, 0.1, 1.4, 0.1
                lf = cfp.legend_formatter(pos=(xpos, ypos+i*dy), length=length)
                cfp.plot(x=bsdat.x, y=bsdat.x**i, label=rf'$\propto r^{i}$', legend_formatter=lf)
            cfp.plot(xlabel="Radius", ylabel=ylabel, xlog=True, ylog=True, save=fig_path+"ekdr_vs_size_frac_dim_M"+MachStr+".pdf")

if __name__ == "__main__":

    # Argument parser setup
    parser = argparse.ArgumentParser(description="Plot fractal dimension.")
    parser.add_argument("-paper", "--paper_plots", action='store_true', default=False, help="Runs all movie frame plots at paper level quality")
    parser.add_argument("-ov", "--overwrite", action='store_true', default=False, help="Overwrite files")
    args = parser.parse_args()


    # Start timing the process
    start_time = timeit.default_timer()

    # loop over simulations
    if args.paper_plots:
        make_paper_plots()

    if args.overwrite:
        out_path = sim_paths[0] + "FracDim/"
        dump_range = [20,100]
        ydat = []
        for d in range(dump_range[0], dump_range[1]+1, 1):
            filename = out_path+"Turb_hdf5_plt_cnt_{:04d}_ekdr_vs_size.pkl".format(d)
            bsdat = dill.load(open(filename, "rb"))
            print(filename,bsdat.y[0])
    else:
        for i, path in enumerate(sim_paths):

            print(f'\nWorking on: {path}', color='cyan')

            N = params(path).N
            t_turb = params(path).t_turb
            Mach = params(path).Mach
            if Mach == 0.2: MachStr = '0p2'
            if Mach == 5:   MachStr = '5'

            # Read data and plot
            ylabel = r"Dissipation rate $\varepsilon_{\textrm{kin}}/(\langle\rho\rangle\,\mathcal{M}^2\, c_{\textrm{s}}^2\,t_{\textrm{turb}}^{-1}$)"
            bsdat = get_ekdr_size_fractal_dim(path)
            cfp.plot(x=bsdat.x, y=bsdat.y, type='scatter', label='sim')
            cfp.plot(x=0.85, y=0.05, text=rf"$\mathcal{{M}} = {Mach}$", normalised_coords=True)
            for i in range(1,4):
                cfp.plot(x=bsdat.x, y=bsdat.x**i, label=rf'$\propto r^{i}$')
            cfp.plot(xlog=True, ylog=True, xlabel=r"Distance $r$", ylabel=ylabel,
                    save=fig_path+"ekdr_vs_size_frac_dim_M"+MachStr+"_N"+str(N)+".pdf")

    # End timing and output the total processing time
    stop_time = timeit.default_timer()
    print("Processing time: {:.2f}s".format(stop_time - start_time))
