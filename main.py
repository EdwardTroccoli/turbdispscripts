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

if args.overwrite==True:
    cfp.print('Working on TimeEvol', color='yellow')
    cfp.run_shell_command('TimeEvol.py -v vstd ekin ekdr ired -int')
    cfp.print('Working on TurbMov', color='yellow')
    cfp.run_shell_command('TurbMov.py -v dens ekin ekdr')
    cfp.print('Working on TurbPDF', color='yellow')
    cfp.run_shell_command('TurbPDF.py -p2 -ov')
    cfp.print('Working on TurbSpectra', color='yellow')
    cfp.run_shell_command('TurbSpectra.py -ov')

