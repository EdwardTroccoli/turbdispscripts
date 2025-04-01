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
    elif variable == "emag":
        vmin = -1e4
        vmax = +1e4
        bw = 20
    elif variable == "ekin":
        vmin = -1e4
        vmax = +1e4
        bw = 20
    elif variable == "injr":
        vmin = -1e4
        vmax = +1e4
        bw = 20
    else:
        print("Variable not implemented.", error=True)
    # loop over FLASH dump files
    for d in range(20, 101):
        filename = "Turb_hdf5_plt_cnt_{:04d}".format(d)
        cfp.run_shell_command(f'mpirun -np 8 pdfs {path+filename} -dset {variable} -vmin {vmin} -vmax {vmax} -bw {bw}')


def compute_2d_pdf(filename, variables, bins=np.array([np.linspace(0,10,20),np.linspace(0,10,20)]), overwrite=False):
    out_filename = filename+"_2Dpdf_"+variables[0]+"_"+variables[1]+".pkl"
    if not os.path.isfile(out_filename) or overwrite:
        # read data
        gg = fl.FlashGG(filename)
        x = gg.ReadVar(dsets=variables)[0].flatten()
        y = gg.ReadVar(dsets=variables)[1].flatten()
        class ret:
            # compute 2D counts
            counts, x_edges, y_edges, binnumber = binned_statistic_2d(x, y, None, statistic='count', bins=bins)
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
        cfp.plot(x=pdf_dat['col1'], y=pdf_dat['col3'])
        cfp.plot(xlabel='Kinetic Energy Dissipation Rate', ylabel='PDF of Kinetic Energy Dissipation Rate', xlim = [-1e4,1e4], ylog = True, save=fig_path+"aver_"+var+'.pdf', legend_loc='best')
    elif var == "emdr":
        cfp.plot(x=pdf_dat['col1'], y=pdf_dat['col3'])
        cfp.plot(xlabel="Magnetic Energy Dissipation Rate", ylabel='PDF of Magnetic Energy Dissipation Rate', xlim = [-1e4,1e4], ylog = True, save=fig_path+"aver_"+var+'.pdf', legend_loc='best') 
    elif var == "injr":
        cfp.plot(x=pdf_dat['col1'], y=pdf_dat['col3'])
        cfp.plot(xlabel="Energy Injection Rate", ylabel='PDF of Energy Injection Rate', xlim = [-1e4,1e4], ylog = True, save=fig_path+"aver_"+var+'.pdf', legend_loc='best') 
    elif var == "emag":
        cfp.plot(x=pdf_dat['col1'], y=pdf_dat['col3'])
        cfp.plot(xlabel="Magnetic Energy", ylabel='PDF of Magnetic Energy', xlim = [0,230], ylog = True, save=fig_path+"aver_"+var+'.pdf', legend_loc='best') 
    elif var == "ekin":
        cfp.plot(x=pdf_dat['col1'], y=pdf_dat['col3'])
        cfp.plot(xlabel="Kinetic Energy", ylabel='PDF of Kinetic Energy', xlim = [0,1550], ylog = True, save=fig_path+"aver_"+var+'.pdf', legend_loc='best')  

if __name__ == "__main__":

    # Argument parser setup
    parser = argparse.ArgumentParser(description="Create and plot PDFs of different variables from FLASH simulation data.")
    var_choices = ["injr", "ekdr", "emdr", "ekin", "emag"]
    parser.add_argument("-ov", "--overwrite", action='store_true', default=False, help="Overwrite files")
    args = parser.parse_args()

    # Start timing the process
    start_time = timeit.default_timer()

    filename = sim_paths[0]+"Turb_hdf5_plt_cnt_{:04d}".format(50)
    variables = ["dens", "ekdr"]
    xbins = np.logspace(-4, 3, 500)
    ybins = np.linspace(-0.5e4, 0.5e4, 500)
    ybins = cfp.symlogspace(-2, 4, num=250)
    bins = [xbins, ybins]
    po = compute_2d_pdf(filename, variables, bins, overwrite=True)

    cfp.plot_map(po.pdf, xedges=po.x_edges, yedges=po.y_edges, xlim=[po.x_edges.min(),po.x_edges.max()], ylim=[po.y_edges.min(),po.y_edges.max()], log=True, xlog=True, show=True)

    stop()

    exit()

    vars = ["ekin", "emag", "injr", "ekdr", "emdr"]

    # loop over simulations
    for path in sim_paths:
        print(f'Working on: {path}', color='green')
        # loop over simulation variables
        for var in vars:
            pdf_aver_file = "aver_"+var+".pdf_data"
            if args.overwrite:
                compute_pdf(path, var) # compute the PDF by calling C++ 'pdf'
                pdf_files = glob.glob(path+"Turb_hdf5_plt_cnt_????_"+var+".pdf_data")
                aver_dat, header_aver = aver_pdf(pdf_files) # average the PDFs
                write_pdf(pdf_aver_file, aver_dat, header_aver) # write averaged PDF

            # plot the PDF
            if not os.path.isdir(fig_path):
                cfp.run_shell_command('mkdir '+fig_path)

            pdf_dat, pdf_header = read_pdf(pdf_aver_file) # read the PDF data
            plot_pdf(pdf_dat)

    # End timing and output the total processing time
    stop_time = timeit.default_timer()
    print("Processing time: {:.2f}s".format(stop_time - start_time))
