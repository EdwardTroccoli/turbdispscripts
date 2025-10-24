#!/usr/bin/env python
# -*- coding: utf-8 -*-
# written by Edward Troccoli, 2025

import cfpack as cfp
import argparse

# Argument parser setup

parser = argparse.ArgumentParser(description="Run analysis scripts for TurbSims")
parser.add_argument("-ov", "--overwrite", action='store_true', default=False, help="Overwrite/Run all analysis on a given TurbSim")
# Parse the command-line arguments
args = parser.parse_args()

def all_figures(overwrite):
    if overwrite: overwrite = "-ov"
    else: overwrite = ""

    cfp.print('Running Globals.py', color='yellow')
    cfp.run_shell_command(f'Globals.py -clean_datfile -compute_vort -compute_2dpdfs -compute_spectra -compute_ekdr_size {overwrite}')
    cfp.print('Running TurbPDF', color='yellow')
    cfp.run_shell_command(f'TurbPDF.py -p2 {overwrite}')
    cfp.print('Running TimeEvol', color='yellow')
    cfp.run_shell_command('TimeEvol.py -paper')
    cfp.print('Running TurbSpectra', color='yellow')
    cfp.run_shell_command('TurbSpectra.py -paper')
    cfp.print('Running TurbFrac', color='yellow')
    cfp.run_shell_command('TurbFrac.py -paper')

all_figures(args.overwrite)
