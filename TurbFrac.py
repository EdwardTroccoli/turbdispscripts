#!/usr/bin/env python
# -*- coding: utf-8 -*-
# written by Edward Troccoli and Christoph Federrath, 2025

from Globals import *
import timeit
import argparse
import cfpack as cfp
import numpy as np
import matplotlib.pyplot as plt
import flashlib as fl
cfp.load_plot_style() # to load cfpack plot style

def make_paper_plots():
        machs = [0.2, 5]
        for imach, mach in enumerate(machs):
            Re = 2500
            k_turb = 2
            Re_coe_sub = np.array([0.10, 0.01])
            Re_coe_sup = np.array([0.33, 0.01])
            if mach == 0.2:
                sims = ["N2048M0p2HDRe2500HP","N1024M0p2HDRe2500", "N512M0p2HDRe2500", "N256M0p2HDRe2500"]
                MachNum = '0p2'
                MachSim = 'Sub'
                l_nu = 1/(np.array([Re_coe_sub[0], Re_coe_sub[0]+Re_coe_sub[1], Re_coe_sub[0]-Re_coe_sub[1]])*Re**(3/4)*k_turb)
                l_nu_pos_x = 0.57
                l_nu_pos_y = 0.12
                fit_pos_x = np.array([0.60, 0.8])
                fit_pos_y = np.array([0.49, 0.67])
                ylabel = r"Dissipation rate $\varepsilon_{\textrm{kin}}/(\langle\rho\rangle\,\mathcal{M}^2\, c_{\textrm{s}}^2\,t_{\textrm{turb}}^{-1}$)"
            if mach == 5:
                sims = ["N2048M5HDRe2500HP", "N1024M5HDRe2500", "N512M5HDRe2500", "N256M5HDRe2500"]
                MachNum = '5'
                MachSim = 'Sup'
                phi = np.array([0.42, 0.09, 0.12])
                psup = np.array([0.49, 0.01, 0.01])
                l_nu = 1/(np.array([Re_coe_sup[0], Re_coe_sup[0]+Re_coe_sup[1], Re_coe_sup[0]-Re_coe_sup[1]])*Re**(2/3)*k_turb)
                l_s = [(1/k_turb)*phi[0]*mach**(-1/psup[0]), (1/k_turb)*(phi[0]-phi[1])*mach**(-1/(psup[0]-psup[1])), (1/k_turb)*(phi[0]+phi[2])*mach**(-1/(psup[0]+psup[2]))]
                l_nu_pos_x = 0.53
                l_nu_pos_y = 0.12
                fit_pos_x = np.array([0.55, 0.78])
                fit_pos_y = np.array([0.61, 0.74])
                ylabel = None
            color = ['black', 'magenta', 'green', 'grey']
            linestyle = ['solid', 'dashed', 'dashdot', 'dotted']
            dy = [0.1, 0.1, 0.1, 0.1]
            # loop over simulations
            for isim, sim in enumerate(sims):
                # get sim parameters
                N = params(sim).N
                Mach = params(sim).Mach
                t_turb = params(sim).t_turb
                if Mach == 0.2: MachStr = '0p2'
                if Mach == 5:   MachStr = '5'

                # Read pkl file from relevant sim directory
                bsdat = get_ekdr_size_fractal_dim("../"+sim+"/")

                # Plot the curve of dissipation against \ell/L
                xpos, ypos, length = 0.7, 0.1, 1.4
                lf = cfp.legend_formatter(pos=(xpos, ypos+isim*dy[isim]), length=length)
                ret = cfp.plot(x=bsdat.x, y=bsdat.y, label=MachSim+str(N),
                                color=color[isim], linestyle=linestyle[isim], legend_formatter=lf
                                )
                if N == 2048: # Only run for this sim res

                    # Line fitting around the dissipation scale (linear relation - fitting is done in log-log space)
                    def model(x, a,b):
                        return a*x+b

                    # Sample relevant data
                    factor = np.sqrt(3) # Goes sqrt(3) below and above of what l_nu is
                    centre = [l_nu[0], 0.3/np.sqrt(3)] # Define the centres for the two regions to find fractal dimension
                    for i in range(0,2): # Loop over the two fitting regions
                        good_ind = ((bsdat.x >= centre[i]/factor) & (bsdat.x <= centre[i]*factor)) # Sample only relevant region around a factor of sqrt(3) around the centre
                        fitting_range_x = bsdat.x[good_ind] # Extract relevant x data
                        fitting_range_y = bsdat.y[good_ind] # Extract relevant y data
                        guesses = {'a': [1, 1.0, 3],'b': [-10.0, 1.0, 10.0]} # Define some bounds on guesses (needed for cfp.fit)

                        # Call fitting function
                        res = cfp.fit(model, np.log(fitting_range_x), np.log(fitting_range_y), params=guesses)

                        # Plot the fitted region, with a scaling label
                        scale_factor = [[1.5e-1, 8e-1], [6e-1, 2.5e-1]] # Define some scaling parameters to move the fits just beneath simulation curve
                        ret = cfp.plot(x = fitting_range_x, y=scale_factor[i][imach]*fitting_range_x**res.popt[0], color = color[isim])
                        cfp.plot(
                            x = fit_pos_x[i], y = fit_pos_y[i], ax=ret.ax(),
                            text=fr'$\propto \ell^{{{res.popt[0]:.2f}\,\pm\,{res.pstd[0]:.2f}}}$',
                            color=color[isim], normalised_coords=True
                        )

            # Plot a vertical line at l_nu with 16th and 84th percentile shaded (both sub and sup)
            ax = ret.ax()
            ax.axvline(x=l_nu[0], color='#d62728', linewidth = 0.9, linestyle = "dotted")
            ax.axvspan(l_nu[1], l_nu[2], color='#d62728', alpha=0.3, linewidth=0)
            cfp.plot(x=l_nu_pos_x, y=l_nu_pos_y, ax=ret.ax(),
                    text=r'$\ell_{\nu}$', color='#d62728', normalised_coords=True)

            # Plot a vertical line at l_s with 16th and 84th percentile shaded (only sup)
            if mach > 1:
                ax.axvline(x=l_s[0], color='darkgoldenrod', linewidth = 0.9, linestyle = "dotted")
                ax.axvspan(l_s[1], l_s[2], color='darkgoldenrod', alpha=0.1, linewidth=0)
                cfp.plot(x=0.39, y=0.12, ax=ret.ax(),
                    text=r'$\ell_{s}$', color = "darkgoldenrod", normalised_coords=True)

            # Define parameters for plotting ell^1, ell^2, ell^3 in the corner
            r = np.logspace(-3.5, -3)
            r0, y0 = 1e-3, 10**(-0.5)
            dy = 0.075

            # Plot the sample curves
            for i in range(1,4):
                cfp.plot(x=0.04, y=0.93 - i*dy, ax=ret.ax(),
                        text=rf'$\propto \ell^{i}$',normalised_coords=True)
                C = y0 / (r0**i)
                y = C * r**i
                cfp.plot(x=r, y=y, color = "#17becf")

            # Final plot in log-log space
            ret = cfp.plot(xlabel=r"$\ell/L$", ylabel=ylabel, ylim=[1e-8,1e0], xlog=True, ylog=True, save=fig_path+"ekdr_vs_size_frac_dim_M"+MachStr+".pdf")

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
