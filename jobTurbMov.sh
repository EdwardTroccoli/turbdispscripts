#!/bin/bash
#PBS -P ek9
#PBS -q normal
#PBS -l walltime=2:00:00
#PBS -l ncpus=8
#PBS -l mem=32GB
#PBS -l storage=scratch/ek9+gdata/ek9
#PBS -l wd
#PBS -N TurbMovies
#PBS -j oe
#PBS -m bea
#PBS -M edward.troccoli@anu.edu.au

TurbMov.py -v dens &>TurbMov.out

# In case we need to run undersubscribed nodes to get more mem per core.
# The following for instance would give us 8GB/core instead of the normal 4GB/core
# for the current job script setting with 48 cores (1 node) requested:
# mpirun -np 24 -map-by numa:SPAN -rank-by slot ./flash4 1>shell.out00 2>&1

