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

# function to compute divv and vort, also checks if they already exist and runs only if they dont
def derivative(filename, var):
    full_path = path + filename

    # first check if we need to generate anything
    with h5py.File(full_path, "r") as f:
        if not all(x in f for x in ['vort', 'divv']):
            # if plot file doesn't already contain -vort and -divv, this will generate them.
            cfp.run_shell_command(f'derivative_var {full_path} -vort -divv')

    # now reopen the file to see new content
    with h5py.File(full_path, "r+") as f:
        # runs only if vort is missing from the plot file
        if 'vort' in var and 'vort' not in f:
            try:
                vx = f["vorticity_x"][:]
                vy = f["vorticity_y"][:]
                vz = f["vorticity_z"][:]
            except KeyError as e:
                raise KeyError(f"Expected vorticity component missing after derivative_var: {e}")
            vort = np.sqrt(vx**2 + vy**2 + vz**2) # compute the magnitude of vorticity
            f.create_dataset("vort", data=vort)

# computes 1d_pdfs using C++ pdfs function
def compute_1d_pdf(path, variable):
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
    # loop over FLASH dump files
    for d in range(20, 101):
        filename = "Turb_hdf5_plt_cnt_{:04d}".format(d)
        if variable in ['emag', 'ekin']:
            cfp.run_shell_command(f'mpirun -np 8 pdfs {path+filename} -dset {variable} -vmin {vmin} -vmax {vmax} -bw {bw} -log')
        elif variable in ["injr", "ekdr", "emdr", "dens"]:
            cfp.run_shell_command(f'mpirun -np 8 pdfs {path+filename} -dset {variable} -vmin {vmin} -vmax {vmax} -bw {bw}')

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
    fname_pkl = filename + "_2Dpdf_" + variables[0] + "_" + variables[1] + "_" + "M" +MachNumber[i] + ".pkl"
    if not os.path.isfile(fname_pkl) or overwrite:
        # read data
        gg = fl.FlashGG(filename)
        x = gg.ReadVar(dsets=variables)[0].flatten()
        y = gg.ReadVar(dsets=variables)[1].flatten()
        class ret:
            # compute 2D counts
            counts, x_edges, y_edges, binnumber = binned_statistic_2d(x, y, np.ones_like(x, dtype=np.float32), statistic='count', bins=bins)
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
    else:
        print("Read '"+fname_pkl+"'", color="green")
        ret = dill.load(open(fname_pkl, "rb"))
    return ret


def plot_2Dpdf(po):
    out_path = path + "PDFs/"
    if not os.path.isdir(out_path):
        cfp.run_shell_command('mkdir '+out_path)
    if po.variables[0] == "dens": 
        xlabel = r"$\rho/\langle\rho\rangle$"
    if po.variables[0] == "divv": 
        xlabel = r"$\nabla\cdot\mathbf{v}$"
    if po.variables[0] == "vort": 
        xlabel = r"$|\nabla\times\mathbf{v}|$"
    if po.variables[1] == "ekdr": 
        ylabel=r'$\varepsilon_{\textrm{kin}}$'
    cfp.plot_map(po.pdf, xedges=po.x_edges, yedges=po.y_edges, xlabel=xlabel, ylabel=ylabel, cmap_label="PDF",
                 log=True, xlog=True, ylog=True, save=out_path+'averaged_2Dpdf_' + var[0] + "_" + var[1] + "_" + "M" +MachNumber[i] +'.pdf')


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
    for i,path in enumerate(sim_paths):

        print(f'Working on: {path}', color='green')

        # creates the figure output dir
        if not os.path.isdir(fig_path):
            cfp.run_shell_command('mkdir '+fig_path)

        # 1D PDFs
        if args.pdf1d:
            # loop over simulation variables
            vars_1Dpdf = ["dens", "ekin", "injr", "ekdr"]
            for var in vars_1Dpdf:
                pdf_aver_file = "aver_1DPDF_"+var+".pdf_data"
                if args.overwrite:
                    compute_1d_pdf(path, var) # compute the PDF by calling C++ 'pdf'
                    if vars_1Dpdf in ['emag', 'ekin']:
                        pdf_files = glob.glob(path+"Turb_hdf5_plt_cnt_????_"+var+".pdf_data_log")
                        aver_dat, header_aver = aver_pdf(pdf_files) # average the PDFs
                        write_pdf(pdf_aver_file, aver_dat, header_aver) # write averaged PDF
                    elif vars_1Dpdf in ["injr", "ekdr", "emdr", "dens"]: 
                        pdf_files = glob.glob(path+"Turb_hdf5_plt_cnt_????_"+var+".pdf_data")
                        aver_dat, header_aver = aver_pdf(pdf_files) # average the PDFs
                        write_pdf(pdf_aver_file, aver_dat, header_aver) # write averaged PDF

                # plot the PDF
                out_path = path + "PDFs/"
                if not os.path.isdir(out_path):
                    cfp.run_shell_command('mkdir '+out_path)

                pdf_dat, pdf_header = read_pdf(pdf_aver_file) # read the PDF data
                plot_1d_pdf(pdf_dat)

        # 2D PDFs
        if args.pdf2d:
            # variables for the 2d pdf plots, can add more.
            vars_2Dpdf = [["vort", "ekdr"],["divv", "ekdr"],["dens", "ekdr"]]
            if "M5" in path: # supersonic
                bins = np.array([np.logspace(-4, 3, 500), np.logspace(-6, 6, 500)])
            elif "M0p5" in path: # subsonic
                bins = np.array([np.logspace(-4, 3, 500), np.logspace(-6, 6, 500)])
            # loop over simulation variables
            for var in vars_2Dpdf:
                fname_pkl = "averaged_2Dpdf_" + var[0] + "_" + var[1] + "_" + "M" +MachNumber[i] + ".pkl"
                if not os.path.isfile(fname_pkl) or args.overwrite:
                    pdf_data = []
                    for d in range(51, 53):
                        filename = "Turb_hdf5_plt_cnt_{:04d}".format(d)
                        if "vort" in var or "divv" in var:
                            derivative(filename,var)
                        po = compute_2d_pdf(path+filename, var, bins=bins, overwrite=True)
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
                plot_2Dpdf(ret)

    # End timing and output the total processing time
    stop_time = timeit.default_timer()
    print("Processing time: {:.2f}s".format(stop_time - start_time))
