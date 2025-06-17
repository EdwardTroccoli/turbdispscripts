#!/usr/bin/env python
# -*- coding: utf-8 -*-
# written by Edward Troccoli & Christoph Federrath, 2025

import cfpack as cfp
import os
import h5py

# create figure output path
fig_path = "../Figures/"
if not os.path.isdir(fig_path):
    cfp.run_shell_command('mkdir '+fig_path)

# sim_paths = ["../N1024M0p2HDRe2500HP/", "../N1024M0p2HDRe2500/", "../N1024M5HDRe2500HP/", "../N1024M5HDRe2500/"]
sim_paths = ["../N1024M0p2HDRe2500/", "../N1024M5HDRe2500/"]

def params(model_name):
    class ret:
        if 'M0p2' in model_name: Mach, MachNum = 0.2, '0p2'
        if 'M5'   in model_name: Mach,MachNum = 5, '5'
        if '1024' in model_name: N = 1024
        if '512'  in model_name: N =  512
        if '256'  in model_name: N =  256
        t_turb = 0.5 / Mach
        # do the color selection here as well, based on the model name / path; needs implementation...
        color = 'black'
    return ret

# function to compute vort, also checks if they already exist and runs only if they dont
def compute_vort(filename, overwrite=False):

    # first check if we need to generate anything
    run_derivative_var = False
    with h5py.File(filename, "r") as f:
        if not all(x in f for x in ['vorticity_x', 'vorticity_y', 'vorticity_z']):
            run_derivative_var = True
    if run_derivative_var or overwrite:
        # if plot file doesn't already contain -vort and this will generate them.
        cfp.run_shell_command(f'mpirun -np 8 derivative_var {filename} -vort')
    
    slice_x_filename = f"{filename}_vorticity_x_slice_z.h5"
    slice_y_filename = f"{filename}_vorticity_y_slice_z.h5"
    slice_z_filename = f"{filename}_vorticity_z_slice_z.h5"
    dump_location = os.path.dirname(slice_x_filename)+'/movie_files/'+os.path.basename(slice_x_filename)
    if not os.path.exists(dump_location) or overwrite:
        cfp.print(f"Generating slice file...",color='magenta')
        cfp.run_shell_command(f'mpirun -np 8 projection {filename} -dset vorticity_x -slice -pixel 1024 1024')
        cfp.run_shell_command(f'mpirun -np 8 projection {filename} -dset vorticity_y -slice -pixel 1024 1024')
        cfp.run_shell_command(f'mpirun -np 8 projection {filename} -dset vorticity_z -slice -pixel 1024 1024')
        cfp.run_shell_command(f'mv {slice_x_filename} {slice_y_filename} {slice_z_filename} '+ os.path.dirname(slice_x_filename)+'/movie_files/')
    else:
        cfp.print(f"Slice file already exists.",color='red')