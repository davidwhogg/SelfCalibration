#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Rory Holmes (MPIA) / David Hogg (NYU)
# 2012
# This script produces the data required to make the plots in the publication:
# Designing Large-Scale Imaging Surveys for a Retrospective
# Relative Photometric Calibration, Rory Holmes, David W. Hogg, Hans-Walter Rix

# The paper_plot.py script utilizes this data to produce the plots

# Make Python 3 compatible
from __future__ import division, print_function

# Standard Modules
import os
import sys
import numpy as np

sys.path.append('./../..')  # add simulator modules to Python path

# Custom Modules
import simulation
import analysis

# Multiprocessing Flag:
# False - then does not use multiprocessing
# int - uses that many separate processes
multi_proc = False

# The four survey directories
survey_files = ['D'] #['A', 'B', 'C', 'D']

# Load the default parameters
dic = eval(open('parameters.py').read())

# Fit true flat-field only once
dic['best_fit_params'] = analysis.best_fit_ff(
                                    dic['FoV'],
                                    dic['ff_samples'],
                                    dic['flat_field_order'],
                                    dic['stop_condition'],
                                    dic['max_iterations'],
                                    verbose=dic['verbose'])

# Create parameter list
parameter_dictionaries = []
for srvy in survey_files:
    dic['survey_file'] = srvy + '.txt'
    dic['data_dir'] = srvy
    parameter_dictionaries.append(dic.copy())

# Perform simulations
if multi_proc:
    os.nice(19)
    from multiprocessing import Pool
    p = Pool(processes=multi_proc)
    results = p.map(simulation.run_sim, parameter_dictionaries)
else:
    results = []
    for params in parameter_dictionaries:
        results.append(simulation.run_sim(params))
