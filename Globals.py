#!/usr/bin/env python
# -*- coding: utf-8 -*-
# written by Edward Troccoli & Christoph Federrath, 2025

import numpy as np

sim_paths = ["../N1024M5HDRe2500/","../N1024M0p2HDRe2500/","../N512M5HDRe2500/","../N512M0p2HDRe2500/","../N256M5HDRe2500/","../N256M0p2HDRe2500/"]
MachNumber = ["5","0.2", "5","0.2", "5","0.2"]
t_turb = np.array([0.1, 2.5, 0.1, 2.5, 0.1, 2.5])
color = ['r','b','g','pink','magenta','yellow']
fig_path = "../Figures/"