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
import gc
import copy
import matplotlib.pyplot as plt
cfp.import_matplotlibrc(fontscale=0.8)

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


@cfp.timer_decorator
def compute_2d_pdf(out_path, filename, variables, bins, norm=[1.0,1.0], overwrite=False):
    from tqdm import tqdm
    from cfpack.mpi import MPI, comm, nPE, myPE
    print("Total number of MPI ranks = "+str(nPE))
    if MPI: comm.Barrier()
    fname_pkl = out_path+os.path.basename(filename)+"_2Dpdf_"+variables[0]+"_"+variables[1]+".pkl"
    if not os.path.isfile(fname_pkl) or overwrite:
        # read data
        gg = fl.FlashGG(filename, verbose=1)
        print("reading x and y data and computing binned_statistic_2d...", color="red")
        dsets = copy.deepcopy(variables)
        if dsets[0] == 'vort': dsets[0] = 'vorticity'
        MyBlocks = gg.GetMyBlocks(myPE, nPE) # domain decomposition
        counts_loc = []
        for b in tqdm(MyBlocks, disable=(myPE!=0), desc=f"[{myPE}]"): # loop over local list of block for the MPI rank
            x, y = gg.ReadBlockVar(b, dsets=dsets)
            x = x.flatten()*norm[0]
            y = y.flatten()*norm[1]
            counts_, x_edges_, y_edges_, binnum_ = binned_statistic_2d(x, y, np.ones_like(x, dtype=np.float32), statistic='count', bins=bins)
            counts_loc.append(counts_)
            del x; del y; del binnum_; del counts_
            gc.collect()
        # MPI reduction operation
        counts_loc = np.sum(np.array(counts_loc), axis=0)
        if MPI:
            counts_all = np.zeros(np.shape(counts_loc))
            comm.Reduce(counts_loc, counts_all, op=MPI.SUM) # master gets total sum
        else:
            counts_all = counts_loc
        if myPE == 0: # only the master PE writes
            class ret:
                # compute 2D counts
                counts, x_edges, y_edges, = counts_all, x_edges_, y_edges_
                # compute bin areas for normalization
                dx = np.diff(x_edges)[0]  # bin width in x
                dy = np.diff(y_edges)[0]  # bin width in y
                bin_area = dx * dy
                # normalize to get PDF (probability density)
                pdf = counts / (np.sum(counts) * bin_area)  # ensures sum(pdf * bin_area) = 1
            # save the data to file
            with open(fname_pkl, "wb") as fobj:
                print("Writing '"+fname_pkl+"'", color="magenta")
                dill.dump(ret, fobj)
    # read from file
    if MPI: comm.Barrier()
    print("Read '"+fname_pkl+"'", color="green")
    ret = dill.load(open(fname_pkl, "rb"))
    return ret


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
    cfp.run_shell_command('shellutils.py pdf_compress -i '+outfile+' -ov')


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
    cfp.plot(x=xlim, y=yfit, xlabel=xlabel, ylabel=ylabel, color='black', linewidth=0.5, linestyle=(0, (1, 5)), save=save_output)


if __name__ == "__main__":

    compute_2d_pdf("../N256M0p2HDRe2500/PDFs/", "../N256M0p2HDRe2500/Turb_hdf5_plt_cnt_0060", variables=["dens","ekdr"], bins=[np.logspace(-4,3,250),np.logspace(-8,6,250)], norm=[1.0,1.0], overwrite=True)
    exit()

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
            cfp.run_shell_command('mkdir '+out_path)

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
            vars_2Dpdf = [["dens", "ekdr"], ["vort", "ekdr"]]
            norms = [[1.0, t_turb/Mach**2], [1.0/(Mach*N), t_turb/Mach**2]]
            # loop over simulation variables
            for ivar, var in enumerate(vars_2Dpdf):
                # set defaults
                do_fit = False
                by_hand_fit = None
                fit_xlim = None
                fit_ylim = None
                # set binning
                if var[1] == "ekdr":
                    bins_y = np.logspace(-8, 6, 250)
                if var[0] == "dens":
                    bins_x = np.logspace(-4, 3, 250)
                    if Mach == 5:
                        by_hand_fit = [1.5, -1.0] # exponent and normalisation of power-law line to draw
                        fit_xlim = [1e-2, 1e1]
                if var[0] == "vort":
                    bins_x = np.logspace(-6, 1, 250)
                    by_hand_fit = [2.0, 2.5] # exponent and normalisation of power-law line to draw
                    fit_xlim = [1e-2, 3e-1]
                fname_pkl = out_path+"aver_2Dpdf_"+var[0]+"_"+var[1]+".pkl"
                if not os.path.isfile(fname_pkl) or args.overwrite:
                    pdf_data = []
                    for d in range(60, 61, 1):
                        filename = "Turb_hdf5_plt_cnt_{:04d}".format(d)
                        if "vort" in var:
                            compute_vort_file(path+filename, overwrite=False)
                        po = compute_2d_pdf(out_path, path+filename, var, bins=[bins_x,bins_y], norm=norms[ivar], overwrite=args.overwrite)
                        pdf_data.append(po.pdf)
                    # setup a class to store edges and the averaged pdf data.
                    class pdat:
                        pdf = np.mean(np.stack(pdf_data, axis=0), axis=0)
                        x_edges = po.x_edges
                        y_edges = po.y_edges
                        variables = var
                    with open(fname_pkl, "wb") as fobj:
                        print("Writing '"+fname_pkl+"'", color="magenta")
                        dill.dump(pdat, fobj)
                else:
                    print("Read '"+fname_pkl+"'", color="green")
                    pdat = dill.load(open(fname_pkl, "rb"))
                # Plot 2D-PDF
                outfile = fname_pkl[:-4]+'_M'+MachStr+'.pdf'
                plot_2Dpdf(outfile, pdat, do_fit=do_fit, by_hand_fit=by_hand_fit, fit_xlim=fit_xlim, fit_ylim=fit_ylim)

    # End timing and output the total processing time
    stop_time = timeit.default_timer()
    print("Processing time: {:.2f}s".format(stop_time - start_time))
