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
from matplotlib.ticker import SymmetricalLogLocator
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


def compute_2d_pdf(filename, variables, bins, overwrite=False):
    #fname_pkl = filename + "_2Dpdf_" + variables[0] + "_" + variables[1] + "_" + "M" +MachNumber[i] + ".pkl"
    fname_pkl = out_path + f"{Path(filename).stem}_2Dpdf_{variables[1]}_{variables[0]}_M{MachNumber[i]}.pkl"
    if not os.path.isfile(fname_pkl) or overwrite:
        # read data
        gg = fl.FlashGG(filename)
        print("reading x data...", color="red")
        x = gg.ReadVar(dsets=variables)[0].flatten()
        print("reading y data...", color="red")
        y = gg.ReadVar(dsets=variables)[1].flatten()
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

def plot_2Dpdf(po, MachNumber):
    out_path = path + "PDFs/"
    save_output = out_path+'averaged_2Dpdf_' + var[1] + "_" + var[0] + "_" + "M" +MachNumber[i] +'.pdf'
    if not os.path.isdir(out_path):
        cfp.run_shell_command('mkdir '+out_path)
    if po.variables[0] == "ekdr":
        xlabel=r'$\varepsilon_{\textrm{kin}}$'
    if po.variables[1] == "dens":
        ylabel = r"$\rho/\langle\rho\rangle$"
        remove_x_ticks = True
        xlabel = None
    if po.variables[1] == "divv":
        ylabel = r"$\nabla\cdot\mathbf{v}$"
        remove_x_ticks = True
        xlabel = None
    if po.variables[1] == "vorticity":
        ylabel = r"$|\nabla\times\mathbf{v}|$"
        remove_x_ticks = False
    if '0p2' == MachNumber[i]:
        Mach = '0.2'
        remove_y_ticks = False
    elif '5' ==  MachNumber[i]:
        Mach = '5'
        remove_y_ticks = True
        ylabel = None
        #if po.variables[1] == "divv":
        #    xlim = [-1e3,1e3]
    ret = cfp.plot_map(po.pdf, xedges=po.x_edges, yedges=po.y_edges, xlabel=xlabel, ylabel=ylabel,
                 log=True, xlog=True, ylog=True)
    if po.variables[0] == "divv":
        ax.set_xticks([-1e4,-1e3,-1e2,-10, 0, 10,1e2,1e3,1e4])
        ax.set_xticklabels([r,'$-10^{-4}$','$-10^{-3}$','$-10^{-2}$','$-10^{-1}$', '$0$', '$10^1$','$10^2$','$10^3$','$10^4$'])
    ax = ret.ax()[0]
    ax.text(0.05, 0.95, rf"$\mathcal{{M}} = {Mach}$", transform=ax.transAxes,
        fontsize=14, color='black', verticalalignment='top',
        bbox=dict(boxstyle="round,pad=0.3", facecolor='gray', alpha=0.5))
    cfp.plot(ax=ret.ax()[0], normalised_coords=True)
    if remove_x_ticks == True:
        ax.set_xticklabels([])
    if remove_y_ticks == True:
        ax.set_yticklabels([])
    if ('0p2' in path) and ( (po.variables[1] == 'dens') or (po.variables[1] == 'divv')):
        cfp.plot(ax=ret.ax()[0], xlabel=xlabel, ylabel=ylabel, normalised_coords=True, 
             save=save_output)
    else:
        line_fitting(po,xlabel,ylabel,save_output)
    cfp.run_shell_command(f'shellutils.py pdf_compress -i {save_output} -overwrite')
    
def line_fitting(po,xlabel,ylabel,save_output):
    x_edges=po.x_edges
    y_edges=po.y_edges
    geometric_mean_x=np.sqrt(x_edges[:-1] * x_edges[1:])
    geometric_mean_y=np.sqrt(y_edges[:-1] * y_edges[1:])

    ridge_x = []
    ridge_y = []
    ridge_weights =[]

    if po.variables[1]=='divv':
        x_edges=po.x_edges[-250:]
        geometric_mean_x=np.sqrt(x_edges[:-1] * x_edges[1:])
        y_edges=po.y_edges[-250:]
        geometric_mean_y=-np.sqrt(y_edges[:-1] * y_edges[1:])
        po.pdf = po.pdf[:, :249]  
        cfp.get_2d_coords() 
    elif po.variables[1]=='vorticity':
        geometric_mean_x=geometric_mean_x[150:225]
        geometric_mean_y=geometric_mean_y[150:225]
        po.pdf = po.pdf[:, 100:225] 

    for k in range(len(geometric_mean_x)):
        col = po.pdf[k, :]
        j = np.argmax(col)
        ridge_x.append(geometric_mean_x[k])
        ridge_y.append(geometric_mean_y[j])
        ridge_weights.append(j)

    def power_law(x, A, B):
        return A * x**B

    fit_result = cfp.fit(power_law, xdat=geometric_mean_x, ydat=geometric_mean_y,
    params={'A': [-500, 1e-1, 500], 'B': [-2, 0, 2]},
    fit_method='ls',
    plot_fit=False, weights=ridge_weights)

    y_fit = power_law(geometric_mean_x, *fit_result.popt)

    if po.variables[1]=='dens':
        cfp.plot(x=geometric_mean_x[50:200], y=y_fit[50:200], xlabel=xlabel, ylabel=ylabel, color='black',
            save=save_output)
    elif po.variables[1]=='divv':
        cfp.plot(x=geometric_mean_x[100:200], y=y_fit[100:200], xlabel=xlabel, ylabel=ylabel, color='black',
            save=save_output)    
    elif po.variables[1]=='vorticity':
        cfp.plot(x=geometric_mean_x, y=y_fit, xlabel=xlabel, ylabel=ylabel, color='black',
            save=save_output)   
        

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
            vars_2Dpdf = [["ekdr", "vorticity"]]# ,["ekdr", "dens"],["ekdr", "vorticity"]
            # loop over simulation variables
            for var in vars_2Dpdf:
                if '0p2' in out_path:
                    if var[1] == "dens":
                        bins_y = np.logspace(-4, 3, 250)
                    if var[1] == "vorticity":
                        bins_y = np.logspace(-4, 3, 250)
                    if var[1] == "divv":
                        bins_y = cfp.symlogspace(-2, 4, 250)
                    if var[0] == "ekdr":
                        bins_x = np.logspace(-7, 2, 250)
                elif '5' in out_path:
                    if var[1] == "dens":
                        bins_y = np.logspace(-4, 3, 250)
                    if var[1] == "vorticity":
                        bins_y = np.logspace(-2, 4, 250)
                    if var[1] == "divv":
                        bins_y = cfp.symlogspace(-2, 4, 250)
                    if var[0] == "ekdr":
                        bins_x = np.logspace(-6, 7, 250)
                fname_pkl = out_path+"averaged_2Dpdf_" + var[1] + "_" + var[0] + "_" + "M" +MachNumber[i] + ".pkl"
                if not os.path.isfile(fname_pkl) or args.overwrite:
                    pdf_data = []
                    for d in range(20, 101, 1):
                        filename = "Turb_hdf5_plt_cnt_{:04d}".format(d)
                        if "vorticity" in var or "divv" in var:
                            compute_divv_vort(path+filename, overwrite=False)
                        po = compute_2d_pdf(path+filename, var, bins=[bins_x,bins_y],overwrite = args.overwrite)
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
                plot_2Dpdf(ret, MachNumber)

    # End timing and output the total processing time
    stop_time = timeit.default_timer()
    print("Processing time: {:.2f}s".format(stop_time - start_time))
