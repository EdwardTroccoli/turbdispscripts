#!/usr/bin/env python
# -*- coding: utf-8 -*-
# written by Edward Troccoli, 2025

import numpy as np
import timeit
import cfpack as cfp
from cfpack.defaults import *
from Globals import *

if __name__ == "__main__":
    # Start timing the process
    start_time = timeit.default_timer()

    for i, path in enumerate(sim_paths):
        print(f'Working on: {path}')
        dat = cfp.read_ascii(path+"Turb.dat")
        time = dat['01_time']/t_turb[i]
        v_disp = dat['#14_rms_velocity']
        ekin = dat['#10_E_kinetic']
        emag = dat['#12_E_magnetic']
        injr = dat['#41_injection_rate']
        ekdr = dat['#42_ekin_diss_rate']
        emdr = dat['#43_emag_diss_rate']
        cfp.plot(x=time, y=v_disp, label=path[3:-1])
    ylabel=r'$\sigma_v$'
    cfp.plot(xlabel=r'$t/t_\mathrm{turb}$', ylabel=ylabel, save=fig_path+'v_disp.pdf')

    # End timing and output the total processing time
    stop_time = timeit.default_timer()
    print("Processing time: {:.2f}s".format(stop_time - start_time))

