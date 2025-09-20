#!/usr/bin/env python
# -*- coding: utf-8 -*-
# written by Edward Troccoli and Christoph Federrath, 2025

import argparse
import numpy as np
import timeit
import os
import cfpack as cfp
from turblib import aver_pdf, write_pdf, read_pdf
from Globals import *
import glob
cfp.load_plot_style() # to load cfpack plot style

# computes 1d_pdfs using C++ pdfs function
def compute_1d_pdf(filename, variable):
    if variable == "ekdr":
        vmin = -1e4
        vmax = +1e4
        bw = 20
    elif variable == "dens":
        vmin = 1e-1
        vmax = 1e1
        bw = 0.002
    elif variable == "ekin":
        vmin = 0.00001
        vmax = +1e4
        bw = 0.009
    elif variable == "injr":
        vmin = -1e4
        vmax = +1e4
        bw = 20
    else:
        print("Variable not implemented.", error=True)
    if variable in ['ekin']:
        cfp.run_shell_command(f'mpirun -np 8 pdfs {filename} -dset {variable} -vmin {vmin} -vmax {vmax} -bw {bw} -log')
    elif variable in ["injr", "ekdr", "dens"]:
        cfp.run_shell_command(f'mpirun -np 8 pdfs {filename} -dset {variable} -vmin {vmin} -vmax {vmax} -bw {bw}')

# plotting function for 1d pdfs
def plot_1d_pdf(pdf_dat):
    if var == "ekdr":
        cfp.plot(x=pdf_dat['col1'], y=pdf_dat['col2'])
        cfp.plot(xlabel='Kinetic Energy Dissipation Rate', ylabel='PDF of Kinetic Energy Dissipation Rate', xlim=[0,1e2],
                 yerr=np.array([pdf_dat['col3'], pdf_dat['col3']]), shaded_err=True, ylog=True, save=out_path+'aver_1DPDF_'+var+'.pdf')
    elif var == "dens":
        cfp.plot(x=pdf_dat['col1'], y=pdf_dat['col2'])
        cfp.plot(xlabel=r"Density ($\rho/\langle\rho\rangle$)", ylabel=r"PDF of Density ($\rho/\langle\rho\rangle$)",
                 yerr=np.array([pdf_dat['col3'], pdf_dat['col3']]), shaded_err=True, xlim=[1e-1,1e1], xlog=True, ylog=True, save=out_path+'aver_1DPDF_'+var+'.pdf')
    elif var == "injr":
        cfp.plot(x=pdf_dat['col1'], y=pdf_dat['col2'])
        cfp.plot(xlabel="Energy Injection Rate", ylabel='PDF of Energy Injection Rate', xlim=[-1e2,1e2],
                 yerr=np.array([pdf_dat['col3'], pdf_dat['col3']]), shaded_err=True, ylog=True, save=out_path+'aver_1DPDF_'+var+'.pdf')
    elif var == "ekin":
        cfp.plot(x=pdf_dat['col1'], y=pdf_dat['col2'], xlabel="Kinetic Energy", ylabel='PDF of Kinetic Energy', xlim=[1e-3,1e2],
                 yerr=np.array([pdf_dat['col3'], pdf_dat['col3']]), shaded_err=True, xlog=True, ylog=True, save=out_path+'aver_1DPDF_'+var+'.pdf')

def plot_2Dpdf(outfile, pdat, do_fit=False, by_hand_fit=None, fit_xlim=None, fit_ylim=None):
    remove_y_ticks = False
    if pdat.variables[1] == "ekdr":
        ylabel = r"Dissipation rate $\varepsilon_{\textrm{kin}}/(\langle\rho\rangle\,\mathcal{M}^2\, c_{\textrm{s}}^2\,t_{\textrm{turb}}^{-1}$)"
    if pdat.variables[0] == "dens":
        xlabel = r"Density $\rho/(\langle\rho\rangle)$"
    if pdat.variables[0] == "vort":
        xlabel = r"Vorticity $|\nabla\times\mathbf{v}|/(\mathcal{M}c_{\textrm{s}}\Delta x^{-1})$"
    if 'M0p2' in outfile:
        MachStr = '0.2'
        cmap_label = None
    if 'M5' in outfile:
        MachStr = '5'
        remove_y_ticks = True
        ylabel = None
        cmap_label = 'PDF'
    ret = cfp.plot_map(pdat.pdf, xedges=pdat.x_edges, yedges=pdat.y_edges, xlabel=xlabel, ylabel=ylabel,
                       log=True, xlog=True, ylog=True, cmap_label=cmap_label)
    ax = ret.ax()[0]
    cfp.plot(x=0.04, y=0.91, ax=ax, text=r"$\mathcal{M}="+MachStr+"$", normalised_coords=True) # Mach number label
    if remove_y_ticks: ax.set_yticklabels([])
    if do_fit or by_hand_fit is not None:
        line_fitting(pdat, xlabel, ylabel, outfile, xlim=fit_xlim, ylim=fit_ylim, by_hand_fit=by_hand_fit)
    else:
        cfp.plot(ax=ax, xlabel=xlabel, ylabel=ylabel, normalised_coords=True, save=outfile)

def line_fitting(po, xlabel, ylabel, save_output, xlim=None, ylim=None, by_hand_fit=None):
    # compute bin centres (take care of symlog as well)
    # --- x ---
    ineg = po.x_edges <= 0
    ipos = po.x_edges >= 0
    x_mid = np.concatenate([-np.sqrt(po.x_edges[ineg][:-1] * po.x_edges[ineg][1:]), np.sqrt(po.x_edges[ipos][:-1] * po.x_edges[ipos][1:])])
    indz = np.where(po.x_edges == 0)[0]
    if len(indz) > 0:
        indz = indz[0]
        x_mid[indz-1] = 0.5*(po.x_edges[indz-1]+po.x_edges[indz+0])
        x_mid[indz+0] = 0.5*(po.x_edges[indz+0]+po.x_edges[indz+1])
    # --- y ---
    ineg = po.y_edges <= 0
    ipos = po.y_edges >= 0
    y_mid = np.concatenate([-np.sqrt(po.y_edges[ineg][:-1] * po.y_edges[ineg][1:]), np.sqrt(po.y_edges[ipos][:-1] * po.y_edges[ipos][1:])])
    indz = np.where(po.y_edges == 0)[0]
    if len(indz) > 0:
        indz = indz[0]
        y_mid[indz-1] = 0.5*(po.y_edges[indz-1]+po.y_edges[indz+0])
        y_mid[indz+0] = 0.5*(po.y_edges[indz+0]+po.y_edges[indz+1])
    # extend to 2D arrays
    X, Y = np.meshgrid(x_mid, y_mid, indexing='ij')
    # define weights based on PDF
    pdf = np.copy(po.pdf)
    weights = pdf.ravel()
    ind = weights > 0
    weights /= weights[ind].min()
    weights[ind] = np.log10(weights[ind])
    ind = weights > 0
    weights = 10**weights[ind]
    # select only relevant data based on the PDF
    xdat = X.ravel()[ind]
    ydat = Y.ravel()[ind]
    # select fit range
    if xlim is None: xlim = [xdat.min(), xdat.max()]
    if ylim is None: ylim = [ydat.min(), ydat.max()]
    indlim = ((xlim[0] <= xdat) & (xdat <= xlim[1])) & ((ylim[0] <= ydat) & (ydat <= ylim[1]))
    # linear function for power-law fit in log-log space
    def linear_func(x, m, t):
        return m * x + t
    # perform fit
    if by_hand_fit is None:
        xfit = xdat[indlim]
        ineg = xfit <= 0
        ipos = xfit  > 0
        xfit[ipos] = +np.log10(+xfit[ipos])
        xfit[ineg] = -np.log10(-xfit[ineg])
        yfit = ydat[indlim]
        ineg = yfit <= 0
        ipos = yfit  > 0
        yfit[ipos] = +np.log10(+yfit[ipos])
        yfit[ineg] = -np.log10(-yfit[ineg])
        fit_result = cfp.fit(linear_func, xfit, yfit, weights=weights[indlim], params={'m': [0, 1, 5], 't': [-10, 1, 10]})
        yfit = 10**linear_func(np.log10(xlim), *fit_result.popt)
    else:
        print("Drawing by-hand/eye fit with m, t = ", by_hand_fit, color='yellow')
        yfit = 10**linear_func(np.log10(xlim), by_hand_fit[0], by_hand_fit[1]) # by-hand fit
    # plot fit
    cfp.plot(x=xlim, y=yfit, xlabel=xlabel, ylabel=ylabel, color='black', linewidth=1, linestyle=(0, (1, 5)), save=save_output)


if __name__ == "__main__":

    # Argument parser setup
    parser = argparse.ArgumentParser(description="Create and plot PDFs of different variables from FLASH simulation data.")
    parser.add_argument("-ov", "--overwrite", action='store_true', default=False, help="Overwrite files")
    parser.add_argument("-p1", "--pdf1d", action='store_true', help="Compute and plot 1D PDFs")
    parser.add_argument("-p2", "--pdf2d", action='store_true', help="Compute and plot 2D PDFs")

    args = parser.parse_args()

    # Start timing the process
    start_time = timeit.default_timer()

    # loop over simulations
    for i, path in enumerate(sim_paths):

        N = params(path).N
        t_turb = params(path).t_turb
        Mach = params(path).Mach
        if Mach == 0.2: MachStr = '0p2'
        if Mach == 5:   MachStr = '5'

        print(f'\nWorking on: {path}', color='cyan')

        # creates the file output dir
        out_path = path + "PDFs/"
        if not os.path.isdir(out_path):
            if myPE == 0: cfp.run_shell_command('mkdir '+out_path)

        # 1D PDFs
        if args.pdf1d:
            # loop over simulation variables
            vars_1Dpdf = ["dens", "ekin", "injr", "ekdr"]
            for var in vars_1Dpdf:
                pdf_aver_file = out_path+"aver_1DPDF_"+var+".pdf_data"
                if args.overwrite:
                    for d in range(20, 101):
                        filename = "Turb_hdf5_plt_cnt_{:04d}".format(d)
                        compute_1d_pdf(path+filename, var) # compute the PDF by calling C++ 'pdf'
                    if vars_1Dpdf in ['emag', 'ekin']:
                        pdf_files = glob.glob(path+"Turb_hdf5_plt_cnt_????_"+var+".pdf_data_log")
                        aver_dat, header_aver = aver_pdf(pdf_files) # average the PDFs
                        if myPE == 0: write_pdf(pdf_aver_file, aver_dat, header_aver) # write averaged PDF
                    elif vars_1Dpdf in ["injr", "ekdr", "emdr", "dens"]: 
                        pdf_files = glob.glob(path+"Turb_hdf5_plt_cnt_????_"+var+".pdf_data")
                        aver_dat, header_aver = aver_pdf(pdf_files) # average the PDFs
                        if myPE == 0: write_pdf(pdf_aver_file, aver_dat, header_aver) # write averaged PDF
                pdf_dat, pdf_header = read_pdf(pdf_aver_file) # read the PDF data
                plot_1d_pdf(pdf_dat)

        # 2D PDFs
        if args.pdf2d:
            # loop over 2D PDF variables
            for vars in vars_2Dpdf:
                # set defaults
                do_fit = False
                by_hand_fit = None
                fit_xlim = None
                fit_ylim = None
                if vars[0] == "dens" and Mach == 5:
                    by_hand_fit = [1.5, -1.0] # exponent and normalisation of power-law line to draw
                    fit_xlim = [1e-2, 1e1]
                if vars[0] == "vort":
                    by_hand_fit = [2.0, 2.75] # exponent and normalisation of power-law line to draw
                    fit_xlim = [1e-2, 3e-1]
                pdat = get_2d_pdf(path, vars)
                # Plot 2D-PDF
                outfile = out_path+"aver_2Dpdf_"+vars[0]+"_"+vars[1]+'_M'+MachStr+'.pdf'
                plot_2Dpdf(outfile, pdat, do_fit=do_fit, by_hand_fit=by_hand_fit, fit_xlim=fit_xlim, fit_ylim=fit_ylim)

    # End timing and output the total processing time
    stop_time = timeit.default_timer()
    print("Processing time: {:.2f}s".format(stop_time - start_time))
