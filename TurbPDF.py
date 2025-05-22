#!/usr/bin/env python
# -*- coding: utf-8 -*-
# written by Edward Troccoli and Christoph Federrath, 2025

import argparse
import numpy as np
import timeit
import os
import cfpack as cfp
from cfpack.defaults import *
from turblib import aver_pdf, write_pdf, read_pdf
from Globals import *
import glob
import flashlib as fl
from scipy.stats import binned_statistic_2d
import dill
import h5py
from pathlib import Path
import gc
import matplotlib.pyplot as plt


# function to compute divv and vort, also checks if they already exist and runs only if they dont
def compute_divv_vort(filename, overwrite=False):
    # first check if we need to generate anything
    run_derivative_var = False
    with h5py.File(filename, "r") as f:
        if not all(x in f for x in ['vorticity_x', 'vorticity_y', 'vorticity_z', 'divv']):
            run_derivative_var = True
    if run_derivative_var or overwrite:
        # if plot file doesn't already contain -vort and -divv, this will generate them.
        cfp.run_shell_command(f'mpirun -np 8 derivative_var {filename} -vort -divv')

# computes 1d_pdfs using C++ pdfs function
def compute_1d_pdf(filename, variable):
    if variable == "ekdr":
        vmin = -1e4
        vmax = +1e4
        bw = 20
    elif variable == "emdr":
        vmin = -1e4
        vmax = +1e4
        bw = 20
    elif variable == "dens":
        vmin = 1e-1
        vmax = 1e1
        bw = 0.002
    elif variable == "emag":
        vmin = 0.00001
        vmax = +1e4
        bw = 0.009
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
    if variable in ['emag', 'ekin']:
        cfp.run_shell_command(f'mpirun -np 8 pdfs {filename} -dset {variable} -vmin {vmin} -vmax {vmax} -bw {bw} -log')
    elif variable in ["injr", "ekdr", "emdr", "dens"]:
        cfp.run_shell_command(f'mpirun -np 8 pdfs {filename} -dset {variable} -vmin {vmin} -vmax {vmax} -bw {bw}')


# plotting function for 1d pdfs
def plot_1d_pdf(pdf_dat):
    if var == "ekdr":
        cfp.plot(x=pdf_dat['col1'], y=pdf_dat['col2'])
        cfp.plot(xlabel='Kinetic Energy Dissipation Rate', ylabel='PDF of Kinetic Energy Dissipation Rate', xlim=[0,1e2],
                 yerr=np.array([pdf_dat['col3'], pdf_dat['col3']]), shaded_err=True, ylog=True, save=out_path+'aver_1DPDF_'+var+'.pdf')
    elif var == "emdr":
        cfp.plot(x=pdf_dat['col1'], y=pdf_dat['col2'])
        cfp.plot(xlabel="Magnetic Energy Dissipation Rate", ylabel='PDF of Magnetic Energy Dissipation Rate',
                 yerr=np.array([pdf_dat['col3'], pdf_dat['col3']]), shaded_err=True, xlim=[-1e4,1e4], ylog=True, save=out_path+'aver_1DPDF_'+var+'.pdf')
    elif var == "dens":
        cfp.plot(x=pdf_dat['col1'], y=pdf_dat['col2'])
        cfp.plot(xlabel=r"Density ($\rho/\langle\rho\rangle$)", ylabel=r"PDF of Density ($\rho/\langle\rho\rangle$)",
                 yerr=np.array([pdf_dat['col3'], pdf_dat['col3']]), shaded_err=True, xlim=[1e-1,1e1], xlog=True, ylog=True, save=out_path+'aver_1DPDF_'+var+'.pdf')
    elif var == "injr":
        cfp.plot(x=pdf_dat['col1'], y=pdf_dat['col2'])
        cfp.plot(xlabel="Energy Injection Rate", ylabel='PDF of Energy Injection Rate', xlim=[-1e2,1e2],
                 yerr=np.array([pdf_dat['col3'], pdf_dat['col3']]), shaded_err=True, ylog=True, save=out_path+'aver_1DPDF_'+var+'.pdf')
    elif var == "emag":
        cfp.plot(x=pdf_dat['col1'], y=pdf_dat['col2'])
        cfp.plot(xlabel="Magnetic Energy", ylabel='PDF of Magnetic Energy', xlim=[0.000001,250],
                 yerr=np.array([pdf_dat['col3'], pdf_dat['col3']]), shaded_err=True, ylog=True, save=out_path+'aver_1DPDF_'+var+'.pdf')
    elif var == "ekin":
        cfp.plot(x=pdf_dat['col1'], y=pdf_dat['col2'], xlabel="Kinetic Energy", ylabel='PDF of Kinetic Energy', xlim=[1e-3,1e2],
                 yerr=np.array([pdf_dat['col3'], pdf_dat['col3']]), shaded_err=True, xlog=True, ylog=True, save=out_path+'aver_1DPDF_'+var+'.pdf')


def compute_2d_pdf(filename, variables, bins, overwrite=False, norm=None):
    #fname_pkl = filename + "_2Dpdf_" + variables[0] + "_" + variables[1] + "_" + "M" +MachNumber[i] + ".pkl"
    fname_pkl = out_path + f"{Path(filename).stem}_2Dpdf_{variables[1]}_{variables[0]}_M{MachNumber[i]}.pkl"
    if not os.path.isfile(fname_pkl) or overwrite:
        # read data
        gg = fl.FlashGG(filename)
        print("reading x data...", color="red")
        x = gg.ReadVar(dsets=variables)[1].flatten()*norm
        print("reading y data...", color="red")
        y = gg.ReadVar(dsets=variables)[0].flatten()*(t_turb[i]/Mach[i]**2)
        print("computing binned_statistic_2d...", color="blue")
        counts_, x_edges_, y_edges_, binnum_ = binned_statistic_2d(x, y, np.ones_like(x, dtype=np.float32), statistic='count', bins=bins)
        del x; del y; del binnum_
        gc.collect()
        class ret:
            # compute 2D counts
            counts, x_edges, y_edges, = counts_, x_edges_, y_edges_
            # compute bin areas for normalization
            dx = np.diff(x_edges)[0]  # bin width in x
            dy = np.diff(y_edges)[0]  # bin width in y
            bin_area = dx * dy
            # normalize to get PDF (probability density)
            pdf = counts / (np.sum(counts) * bin_area)  # ensures sum(pdf * bin_area) = 1
        # save the data to file
        with open(fname_pkl, "wb") as fobj:
            print("Writing '"+fname_pkl+"'", color="magenta")
            dill.dump(ret, fobj, protocol = 4)
    else:
        print("Read '"+fname_pkl+"'", color="green")
        ret = dill.load(open(fname_pkl, "rb"))
    return ret

def plot_2Dpdf(po, MachNumber, do_fit=False, by_hand_fit=False, fit_xlim=None, fit_ylim=None):
    out_path = path + "PDFs/"
    save_output = out_path+'averaged_2Dpdf_' + var[1] + "_" + var[0] + "_" + "M" +MachNumber[i] +'.pdf'
    if not os.path.isdir(out_path):
        cfp.run_shell_command('mkdir '+out_path)
    remove_x_ticks,remove_y_ticks=None,None
    if po.variables[0] == "ekdr":
        ylabel = r"Dissipation rate $\varepsilon_{\textrm{kin}}/(\langle\rho\rangle\,\mathcal{M}^2\, c_{\textrm{s}}^2\,t_{\textrm{turb}}^{-1}$)"
    if po.variables[1] == "dens":
        xlabel = r"$\rho/(\langle\rho\rangle)$"
        remove_x_ticks = True
    if po.variables[1] == "divv":
        xlabel = r"$\nabla\cdot\mathbf{v}/(\mathcal{M}c_{\textrm{s}}\Delta x^{-1})$"
        remove_x_ticks = True
    if po.variables[1] == "vorticity":
        xlabel = r"$|\nabla\times\mathbf{v}|/(\mathcal{M}c_{\textrm{s}}\Delta x^{-1})$"
        remove_x_ticks = False
    if '0p2' == MachNumber[i]:
        Mach = '0.2'
        remove_y_ticks = False
        remove_x_ticks = False
        cmap_label = None
    elif '5' ==  MachNumber[i]:
        Mach = '5'
        remove_y_ticks = True
        remove_x_ticks = False
        ylabel = None
        cmap_label = 'PDF'
    ret = cfp.plot_map(po.pdf, xedges=po.x_edges, yedges=po.y_edges, xlabel=xlabel, ylabel=ylabel,
                 log=True, xlog=True, ylog=True, cmap_label=cmap_label)
    ax = ret.ax()[0]
    ax.text(0.05, 0.95, rf"$\mathcal{{M}} = {Mach}$", transform=ax.transAxes,
        fontsize=14, color='black', verticalalignment='top',
        bbox=dict(boxstyle="round,pad=0.3", facecolor='gray', alpha=0.5))
    if po.variables[1] == "divv":
        ax.set_xticks([-1e0,-1e-1,-1e-2,-1e-3, 0, 1e-3,1e-2,1e-1,1e0])
        ax.set_xticklabels(['$-10^0$','$-10^{-1}$','$-10^{-2}$','$-10^{-3}$','$0$', '$10^{-3}$','$10^{-2}$','$10^{-1}$','$10^0$'])
    cfp.plot(ax=ret.ax()[0], normalised_coords=True)
    if remove_x_ticks == True:
        ax.set_xticklabels([])
    if remove_y_ticks == True:
        ax.set_yticklabels([])
    if do_fit:
        line_fitting(po, xlabel, ylabel, save_output, xlim=fit_xlim, ylim=fit_ylim, by_hand_fit=by_hand_fit)
    else:
        cfp.plot(ax=ret.ax()[0], xlabel=xlabel, ylabel=ylabel, normalised_coords=True, save=save_output)
    cfp.run_shell_command(f'shellutils.py pdf_compress -i {save_output} -overwrite')


def line_fitting(po, xlabel, ylabel, save_output, xlim=None, ylim=None, by_hand_fit=False):
    geo_x = np.sqrt(po.x_edges[:-1] * po.x_edges[1:])
    geo_y = np.sqrt(po.y_edges[:-1] * po.y_edges[1:])
    pdf = np.copy(po.pdf)

    if po.variables[1]=='divv':
        x_edges = x_edges[x_edges < 0]
        geo_x = -np.sqrt(x_edges[:-1] * x_edges[1:])
        pdf = pdf[:len(x_edges), :]

    X, Y = np.meshgrid(geo_x, geo_y, indexing='ij')


    def linear_func(x, m, t):
        return m * x + t

    # define weights based on PDF
    weights = np.copy(pdf).ravel()
    ind = weights > 0
    weights /= weights[ind].min()
    weights[ind] = np.log10(weights[ind])
    ind = weights > 0
    weights = weights[ind]

    # selecting only relevant data based on the PDF
    xdat = X.ravel()[ind]
    ydat = Y.ravel()[ind]

    # selecting fit range
    if xlim is None: xlim = [xdat.min(), xdat.max()]
    if ylim is None: ylim = [ydat.min(), ydat.max()]
    indlim = ((xlim[0] <= xdat) & (xdat <= xlim[1])) & ((ylim[0] <= ydat) & (ydat <= ylim[1]))

    if not by_hand_fit:
        fit_result = cfp.fit(linear_func, np.log10(xdat[indlim]), np.log10(ydat[indlim]),
                            weights=weights[indlim], params={'m': [0, 1, 5], 't': [-10, 1, 10]})
        y_fit = 10**linear_func(np.log10(xlim), *fit_result.popt)

    # plot the mode
    # indmax = np.array(np.where(pdf==pdf.max())).flatten()
    # cfp.plot(x=X[indmax[0],indmax[1]], y=Y[indmax[0],indmax[1]], type='scatter', color='magenta')

    if by_hand_fit:
        if po.variables[1]=='vorticity':
            cfp.plot(x=xlim, y=10**linear_func(np.log10(xlim), 2.0, 2.5), xlabel=xlabel, ylabel=ylabel,
                     color='black', linewidth=0.5, linestyle=(0, (1, 5)), save=save_output)
    else:
        cfp.plot(x=xlim, y=y_fit, xlabel=xlabel, ylabel=ylabel,
                 color='black', linewidth=0.5, linestyle='dashed', save=save_output)

    if po.variables[1]=='dens':
        cfp.plot(x=geo_x, y=y_fit, xlabel=xlabel, ylabel=ylabel, color='black',
            save=save_output)
    elif po.variables[1]=='divv':
        cfp.plot(x=geo_x[100:200], y=y_fit[100:200], xlabel=xlabel, ylabel=ylabel, color='black',
            save=save_output)

    stop()

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
        if '1024' in path:
            N=1024
        elif '512' in path:
            N=512
        elif '256' in path:
            N=256
        print(f'Working on: {path}', color='green')

        # creates the file output dir
        out_path = path + "PDFs/"
        if not os.path.isdir(out_path):
            cfp.run_shell_command('mkdir '+out_path)

        # 1D PDFs
        if args.pdf1d:
            # loop over simulation variables
            vars_1Dpdf = ["dens", "ekin", "injr", "ekdr"]
            for var in vars_1Dpdf:
                pdf_aver_file = "aver_1DPDF_"+var+ "_" + "M" +MachNumber[i] + ".pdf_data"
                if args.overwrite:
                    for d in range(20, 101):
                        filename = "Turb_hdf5_plt_cnt_{:04d}".format(d)
                        compute_1d_pdf(path+filename, var) # compute the PDF by calling C++ 'pdf'
                    if vars_1Dpdf in ['emag', 'ekin']:
                        pdf_files = glob.glob(path+"Turb_hdf5_plt_cnt_????_"+var+".pdf_data_log")
                        aver_dat, header_aver = aver_pdf(pdf_files) # average the PDFs
                        write_pdf(pdf_aver_file, aver_dat, header_aver) # write averaged PDF
                    elif vars_1Dpdf in ["injr", "ekdr", "emdr", "dens"]: 
                        pdf_files = glob.glob(path+"Turb_hdf5_plt_cnt_????_"+var+".pdf_data")
                        aver_dat, header_aver = aver_pdf(pdf_files) # average the PDFs
                        write_pdf(pdf_aver_file, aver_dat, header_aver) # write averaged PDF

                pdf_dat, pdf_header = read_pdf(pdf_aver_file) # read the PDF data
                plot_1d_pdf(pdf_dat)

        # 2D PDFs
        if args.pdf2d:
            # variables for the 2d pdf plots, can add more.
            vars_2Dpdf = [["ekdr", "dens"], ["ekdr", "vorticity"], ["ekdr", "divv"]]# ["ekdr", "dens"], ["ekdr", "vorticity"],
            # loop over simulation variables
            for var in vars_2Dpdf:
                # set defaults
                norm = 1
                do_fit = False
                by_hand_fit = False
                fit_xlim = None
                fit_ylim = None
                # set normalisation
                if var[1] == "divv" or var[1] == "vorticity":
                    norm = (1/N) / Mach[i]
                # set binning
                if var[0] == "ekdr":
                    bins_y = np.logspace(-8, 6, 250)
                if var[1] == "dens":
                    bins_x = np.logspace(-4, 3, 250)
                if var[1] == "vorticity":
                    bins_x = np.logspace(-6, 1, 250)
                    do_fit = True
                    by_hand_fit = True
                    fit_xlim = [1e-2, 1]
                if var[1] == "divv":
                    bins_x = cfp.symlogspace(-4, 0, 250)
                fname_pkl = out_path+"averaged_2Dpdf_" + var[1] + "_" + var[0] + "_" + "M" +MachNumber[i] + ".pkl"
                if not os.path.isfile(fname_pkl) or args.overwrite:
                    pdf_data = []
                    for d in range(20, 101, 1):
                        filename = "Turb_hdf5_plt_cnt_{:04d}".format(d)
                        if "vorticity" in var or "divv" in var:
                            compute_divv_vort(path+filename, overwrite=False)
                        po = compute_2d_pdf(path+filename, var, bins=[bins_x,bins_y],overwrite = args.overwrite, norm=norm)
                        pdf_data.append(po.pdf)
                    # setup a class to store edges and the averaged pdf data.
                    class ret:
                        pdf = np.mean(np.stack(pdf_data, axis=0), axis=0)
                        x_edges = po.x_edges
                        y_edges = po.y_edges
                        variables = var
                    with open(fname_pkl, "wb") as fobj:
                        print("Writing '"+fname_pkl+"'", color="magenta")
                        dill.dump(ret, fobj)
                else:
                    print("Read '"+fname_pkl+"'", color="green")
                    ret = dill.load(open(fname_pkl, "rb"))
                # Plot
                plot_2Dpdf(ret, MachNumber, do_fit=do_fit, by_hand_fit=by_hand_fit, fit_xlim=fit_xlim, fit_ylim=fit_ylim)

    # End timing and output the total processing time
    stop_time = timeit.default_timer()
    print("Processing time: {:.2f}s".format(stop_time - start_time))
