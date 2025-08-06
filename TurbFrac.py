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


if __name__ == "__main__":

    # Argument parser setup
    parser = argparse.ArgumentParser(description="Plot fractal dimension.")
    args = parser.parse_args()

    # loop over simulations
    for i, path in enumerate(sim_paths):

        print(f'\nWorking on: {path}', color='cyan')

        N = params(path).N
        t_turb = params(path).t_turb
        Mach = params(path).Mach
        if Mach == 0.2: MachStr = '0p2'
        if Mach == 5:   MachStr = '5'


        # Read data and plot
        outfile = path+'FracDim/aver_ekdr_vs_size'+MachStr+'.pdf'
        ylabel = r"Dissipation rate $\varepsilon_{\textrm{kin}}/(\langle\rho\rangle\,\mathcal{M}^2\, c_{\textrm{s}}^2\,t_{\textrm{turb}}^{-1}$)"
        bsdat = get_ekdr_size_fractal_dim(path)
        scaling = [1e2, 1e-3]
        cfp.plot(x=bsdat.x, y=bsdat.y, type='scatter', label='sim')
        cfp.plot(x=bsdat.x, y=scaling[i]*bsdat.x**1, label=r'$\propto r^1$')
        cfp.plot(x=bsdat.x, y=scaling[i]*bsdat.x**2, label=r'$\propto r^2$')
        cfp.plot(x=bsdat.x, y=scaling[i]*bsdat.x**3, label=r'$\propto r^3$')
        cfp.plot(xlog=True, ylog=True, xlabel=r"Distance $r$", ylabel=ylabel)
        cfp.plot(x=0.75, y=0.95, text=rf"$\mathcal{{M}} = {Mach}$", save=fig_path+"ekdr_vs_size_frac_dim_M"+MachStr+".pdf")

        # fit
        slope, intercept, r_value, p_value, std_err=linregress(bsdat.x, bsdat.y)

