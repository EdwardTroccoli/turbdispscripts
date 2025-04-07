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


def compute_pdf(path, variable):
    if variable == "ekdr":
        vmin = -1e4
        vmax = +1e4
        bw = 20
    elif variable == "emdr":
        vmin = -1e4
        vmax = +1e4
        bw = 20
    elif variable == "dens":
        vmin = 10**-0.5
        vmax = +10**0.5
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
        elif variable in ["injr", "ekdr", "emdr","dens"]: 
            cfp.run_shell_command(f'mpirun -np 8 pdfs {path+filename} -dset {variable} -vmin {vmin} -vmax {vmax} -bw {bw}')

def compute_2d_pdf(filename, variables, bins, overwrite=False):
    out_filename = filename + "_2Dpdf_" + variables[0][0] + "_" + variables[0][1] + ".pkl"
    if not os.path.isfile(out_filename) or overwrite:
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
        with open(out_filename, "wb") as fobj:
            print("Writing '"+out_filename+"'", color="magenta")
            dill.dump(ret, fobj)
    else:
        print("Read '"+out_filename+"'", color="green")
        ret = dill.load(open(out_filename, "rb"))
    return ret


def plot_pdf(pdf_dat):
    if var == "ekdr":
        cfp.plot(x=pdf_dat['col1'], y=pdf_dat['col2'])
        cfp.plot(xlabel='Kinetic Energy Dissipation Rate', ylabel='PDF of Kinetic Energy Dissipation Rate', xlim=[0,1e2], yerr=np.array([pdf_dat['col3'], pdf_dat['col3']]), shaded_err=True, ylog=True,  save=out_path+'aver_1DPDF_'+var+'.pdf')
    elif var == "emdr":
        cfp.plot(x=pdf_dat['col1'], y=pdf_dat['col2'])
        cfp.plot(xlabel="Magnetic Energy Dissipation Rate", ylabel='PDF of Magnetic Energy Dissipation Rate', yerr=np.array([pdf_dat['col3'], pdf_dat['col3']]), shaded_err=True, xlim=[-1e4,1e4], ylog=True,  save=out_path+'aver_1DPDF_'+var+'.pdf') 
    elif var == "dens":
        cfp.plot(x=pdf_dat['col1'], y=pdf_dat['col2'])
        cfp.plot(xlabel=r"Density ($\rho/\langle\rho\rangle$)", ylabel=r"PDF of Density ($\rho/\langle\rho\rangle$)", yerr=np.array([pdf_dat['col3'], pdf_dat['col3']]), shaded_err=True, xlim=[1e-1,1e1], xlog=True, ylog=True,  save=out_path+'aver_1DPDF_'+var+'.pdf') 
    elif var == "injr":
        cfp.plot(x=pdf_dat['col1'], y=pdf_dat['col2'])
        cfp.plot(xlabel="Energy Injection Rate", ylabel='PDF of Energy Injection Rate', xlim=[-1e2,1e2], yerr=np.array([pdf_dat['col3'], pdf_dat['col3']]), shaded_err=True, ylog=True,  save=out_path+'aver_1DPDF_'+var+'.pdf') 
    elif var == "emag":
        cfp.plot(x=pdf_dat['col1'], y=pdf_dat['col2'])
        cfp.plot(xlabel="Magnetic Energy", ylabel='PDF of Magnetic Energy', xlim=[0.000001,250], yerr=np.array([pdf_dat['col3'], pdf_dat['col3']]), shaded_err=True, ylog=True,  save=out_path+'aver_1DPDF_'+var+'.pdf') 
    elif var == "ekin":
        cfp.plot(x=pdf_dat['col1'], y=pdf_dat['col2'], xlabel="Kinetic Energy", ylabel='PDF of Kinetic Energy', xlim=[1e-3,1e2], yerr=np.array([pdf_dat['col3'], pdf_dat['col3']]), shaded_err=True, xlog=True, ylog=True, save=out_path+'aver_1DPDF_'+var+'.pdf')  

def plot_2Dpdf(po):
    out_path = path + "PDFs/"
    if not os.path.isdir(out_path):
        cfp.run_shell_command('mkdir '+out_path)
    cfp.plot_map(po.pdf, xedges=po.x_edges, yedges=po.y_edges, xlabel=r"Density ($\rho/\langle\rho\rangle$)", ylabel=r'$\varepsilon_{\textrm{kin}}$', log=True, xlog=True, ylog=True, save=out_path+'aver_2DPDF_'+variables[0][0]+"_"+variables[0][1]+'.pdf')
     

if __name__ == "__main__":

    # Argument parser setup
    parser = argparse.ArgumentParser(description="Create and plot PDFs of different variables from FLASH simulation data.")
    #var_choices = ["injr", "ekdr", "emdr", "ekin", "emag"]
    parser.add_argument("-ov", "--overwrite", action='store_true', default=False, help="Overwrite files")
    parser.add_argument("-p1", "--pdf1d", action='store_true', help="Compute and plot 1D PDFs")
    parser.add_argument("-p2", "--pdf2d", action='store_true', help="Compute and plot 2D PDFs")

    args = parser.parse_args()

    # Start timing the process
    start_time = timeit.default_timer()
    
    #variables for the 1d pdf plots. yet to implement dens.
    vars = ["dens","ekin", "injr", "ekdr"]
    #variables for the 2d pdf plots, can add more.
    variables = [["dens", "ekdr"]]

    # loop over simulations
    if args.pdf1d:
        for path in sim_paths:
            print(f'Working on: {path}', color='green')
            # loop over simulation variables
            for var in vars:
                pdf_aver_file = "aver_1DPDF_"+var+".pdf_data"
                if args.overwrite:
                    compute_pdf(path, var) # compute the PDF by calling C++ 'pdf'
                    if var in ['emag', 'ekin']:
                        pdf_files = glob.glob(path+"Turb_hdf5_plt_cnt_????_"+var+".pdf_data_log")
                        aver_dat, header_aver = aver_pdf(pdf_files) # average the PDFs
                        write_pdf(pdf_aver_file, aver_dat, header_aver) # write averaged PDF
                    elif var in ["injr", "ekdr", "emdr", "dens"]: 
                        pdf_files = glob.glob(path+"Turb_hdf5_plt_cnt_????_"+var+".pdf_data")
                        aver_dat, header_aver = aver_pdf(pdf_files) # average the PDFs
                        write_pdf(pdf_aver_file, aver_dat, header_aver) # write averaged PDF

                # plot the PDF
                out_path = path + "PDFs/"
                if not os.path.isdir(out_path):
                    cfp.run_shell_command('mkdir '+out_path)

                pdf_dat, pdf_header = read_pdf(pdf_aver_file) # read the PDF data
                plot_pdf(pdf_dat)

    if args.pdf2d:
        for path in sim_paths:
            print(f'Working on: {path}', color='green')
            if "M5" in path:
                bins=np.array([np.logspace(-4, 5, 500), np.logspace(-4, 4, 500)])
            elif "M0p5" in path:
                bins=np.array([np.logspace(-4, 5, 500), np.logspace(-4, 4, 500)])
            # loop over simulation variables
            for var in variables:
                if args.overwrite:
                    pdf_data = []
                    for d in range(50, 51):
                        filename = "Turb_hdf5_plt_cnt_{:04d}".format(d)
                        po = compute_2d_pdf(path+filename, variables, bins=bins, overwrite=True)
                        pdf_data.append(po.pdf)   
                    #setup a class to store edges and the averaged pdf data.    
                    class PO:
                        pass
                    po_avg = PO()
                    po_avg.pdf = np.mean(np.stack(pdf_data, axis=0), axis=0)
                    po_avg.x_edges = po.x_edges
                    po_avg.y_edges = po.y_edges
                    # Save to file
                    out_filename = "averaged_2Dpdf_" + variables[0][0] + "_" + variables[0][1] + ".pkl"
                    with open(out_filename, "wb") as f:
                        dill.dump(po_avg, f)
                # plot the PDF
                if not os.path.isdir(fig_path):
                    cfp.run_shell_command('mkdir '+fig_path)
                with open("averaged_2Dpdf_" + variables[0][0] + "_" + variables[0][1] + ".pkl", "rb") as f:
                    po_loaded = dill.load(f)
                # Plot
                plot_2Dpdf(po_loaded)   

    # End timing and output the total processing time
    stop_time = timeit.default_timer()
    print("Processing time: {:.2f}s".format(stop_time - start_time))
