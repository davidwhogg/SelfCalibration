#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Rory Holmes (MPIA) / David Hogg (NYU)
# 2012
# This script reproduces the plots in the publication on these self-calibration
# simulations: Designing Large-Scale Imaging Surveys for a Retrospective
# Relative Photometric Calibration, Rory Holmes, David W. Hogg, Hans-Walter Rix
# The plotting script itself is contained in the paper_plot.py script

# Make Python 3 compatible
from __future__ import division, print_function

# Standard Modules
import os
import sys
import numpy as np

sys.path.append('./../..')  # add simulator modules to Python path

# Multiprocessing Flag:
# False - then does not use multiprocessing
# int - uses that many separate processes
multi_proc = 4

# Custom Modules
import simulation
import analysis

# Load the default parameters
dic = eval(open('parameters.py').read())

survey_files = ['A', 'B', 'C', 'D']

# Fit true flat-field only once
dic['best_fit_params'] = analysis.best_fit_ff(
                                    dic['FoV'],
                                    dic['ff_samples'],
                                    dic['flat_field_order'],
                                    dic['stop_condition'],
                                    dic['max_iterations'],
                                    verbose=dic['verbose'])
    
parameter_dictionaries = []
for srvy in survey_files:
    dic['survey_file'] = srvy + '.txt'
    dic['data_dir'] = srvy
    parameter_dictionaries.append(dic.copy())

if multi_proc:
    os.nice(19)
    from multiprocessing import Pool
    p = Pool(processes=multi_proc)
    results = p.map(simulation.run_sim, parameter_dictionaries)
else:
    results = []
    for params in parameter_dictionaries:
        results.append(simulation.run_sim(params))



