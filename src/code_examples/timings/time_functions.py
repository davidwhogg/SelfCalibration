#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Rory Holmes (MPIA) / David Hogg (NYU)
# 2012
# This script performs the self-calibration simulation with different source
# densities, so that we can see how the time required by the self-calibration
# procedure scales.

# Make Python 3 compatible
from __future__ import division, print_function

# Standard Modules
import os
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')  # Allow plots to be made without X11
import matplotlib.pyplot as plt

# Add simulator modules to Python path
sys.path.append('./../..')
import analysis
import simulation

# Range of sources to run self-calibration simulations with
# Sort high to low, to increase simulation speed
density_of_sources = np.logspace(0, 3, 30)[::-1]

# Multiprocessing flat: int = number of threads, False = no multiprocessing
multi_proc = 2

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
dirs = []
for density in density_of_sources:
    dic['density_of_stars'] = density
    dic['data_dir'] = '{0:.3f}'.format(density)
    dirs.append('{0:.3f}'.format(density))
    parameter_dictionaries.append(dic.copy())

if multi_proc:
    os.nice(19)
    from multiprocessing import Pool
    p = Pool(processes=multi_proc)
    p.map(simulation.run_sim, parameter_dictionaries)
else:
    for params in parameter_dictionaries:
        simulation.run_sim(params)

# Read timings from individual directories
density = np.zeros(len(dirs))
time_measurement = np.zeros(len(dirs))
time_self_cal = np.zeros(len(dirs))

for indx in range(len(dirs)):
    f = open(dirs[indx] + '/timings.txt', 'r')
    density[indx] = float(dirs[indx])
    time_measurement[indx] = float(f.readline()
                                        .strip('Measurement Catalog (s): '))
    time_self_cal[indx] = float(f.readline()
                                        .strip('Self-Calibration (s): '))

fig = plt.figure(figsize=(6, 6))
ax = fig.add_axes([0.15, 0.15, 0.8, 0.8])
ax.loglog(density, time_self_cal, 'kx')
ax.set_xlabel('Density of Sources (deg-2)')
ax.set_ylabel('Time for Self-Calibration to Find a Solution (s)')
fig.savefig('self_cal_time.png')
fig.clf()

fig = plt.figure(figsize=(6, 6))
ax = fig.add_axes([0.15, 0.15, 0.8, 0.8])
ax.loglog(density, time_measurement, 'kx')
ax.set_xlabel('Density of Sources (deg-2)')
ax.set_ylabel('Time to Perform Survey (s)')
fig.savefig('survey_time.png')
fig.clf()
