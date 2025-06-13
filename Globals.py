#!/usr/bin/env python
# -*- coding: utf-8 -*-
# written by Edward Troccoli & Christoph Federrath, 2025

import cfpack as cfp
import os

# create figure output path
fig_path = "../Figures/"
if not os.path.isdir(fig_path):
    cfp.run_shell_command('mkdir '+fig_path)

# sim_paths = ["../N1024M0p2HDRe2500HP/", "../N1024M0p2HDRe2500/", "../N1024M5HDRe2500HP/", "../N1024M5HDRe2500/"]
sim_paths = ["../N1024M0p2HDRe2500/", "../N1024M5HDRe2500/"]

def params(model_name):
    class ret:
        if 'M0p2' in model_name: Mach = 0.2
        if 'M5' in model_name: Mach = 5
        t_turb = 0.5 / Mach
        # do the color selection here as well, based on the model name / path; needs implementation...
        color = 'black'
    return ret
