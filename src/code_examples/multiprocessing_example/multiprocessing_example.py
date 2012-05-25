#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Rory Holmes (MPIA) / David Hogg (NYU)
# 2012
# This script utilizes the multiprocessing aspects of the self-calibration
# functions. It calculates the performance of self-calibration with different
# source densities and plots the output. A separate self-calibration simulation
# is run for each source density value. To increase the speed of the simulation
# the best-in-basis fit to the *true* instrument response is calculated only 
# once at the beginning and is then provided as an input to each 
# self-calibration simulation.

# Make Python 3 compatible
from __future__ import division, print_function

# Standard Modules
import os
import sys
import numpy as np
import matplotlib.pyplot as plt

# Add simulator modules to Python path
sys.path.append('./../..')  

# Multiprocessing Flag:
# False - then does not use multiprocessing
# int - uses that many separate processes
multi_proc = 4

# Custom Modules
import simulation
import analysis

# Range of sources to run self-calibration simulations with
# Sort high to low, to increase simulation speed
density_of_sources = np.logspace(0, 2, 3)[::-1]

# Load the default parameters
dic = eval(open('parameters.py').read())

# Fit the best-in-basis flat-field only once
dic['best_fit_params'] = analysis.best_fit_ff(
                                    dic['FoV'],
                                    dic['ff_samples'],
                                    dic['flat_field_order'],
                                    dic['stop_condition'],
                                    dic['max_iterations'],
                                    verbose=dic['verbose'])


parameter_dictionaries = []
for density in density_of_sources:
    dic['density_of_stars'] = density
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

iterations = np.array([])
rms = np.array([])
badness = np.array([])
best_in_basis_badness = np.array([])
chi2 = np.array([])
for run in results:
    iterations = np.append(iterations, run[0])
    rms = np.append(rms, run[1])
    badness = np.append(badness, run[2])
    best_in_basis_badness = np.append(best_in_basis_badness, run[3])
    chi2 = np.append(chi2, run[4])

fig = plt.figure()
ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
ax.loglog(density_of_sources, badness, '0.5', label='Badness')
ax.loglog(density_of_sources, best_in_basis_badness, 'ko', markersize=1, label='Best-in-Basis Badness')
ax.set_xlabel('Density of Sources (deg$^{-2}$)')
ax.set_ylabel('Badness (%)')
ax.legend(loc='lower left')
fig.savefig('badnesses.pdf')

