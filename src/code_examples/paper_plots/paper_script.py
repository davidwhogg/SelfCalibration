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
import pickle

sys.path.append('./../..')  # add simulator modules to Python path

# Custom Modules
import simulation
import analysis

# Multiprocessing Flag:
# False - then does not use multiprocessing
# int - uses that many separate processes
multi_proc = 4

# The four survey directories
survey_files = ['A', 'B', 'C', 'D']

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

#  Since there are lots of souces at the edge of Survey D that are only
#  measured a few times, the RMS 
f = open('rms_results.txt', 'w')
survey_area = [-4., 4., -4., 4.]
header = 'Source RMS Error with {0} survey area:'.format(survey_area)
f.write(header + '\n')
print(header)

for path in survey_files:
	
	# Load the catalogs from the simulation run
	fitted_cat = pickle.load(open((path + '/fitted_catalog.p'), mode='rb'))
	true_cat = pickle.load(open((path + '/source_catalog.p'), mode='rb'))
	
	if true_cat.k.max() > fitted_cat.k.max():
		k_max = true_cat.k.max()
	else:
		k_max = fitted_cat.k.max()
	
	error = np.array([])
	norm_error = np.array([])
	
	true_flux = np.array([])
	fitted_flux = np.array([])
	
	inside = ((true_cat.alpha > survey_area[0])
	                * (true_cat.alpha < survey_area[1])
	                * (true_cat.beta > survey_area[2])
	                * (true_cat.beta < survey_area[3]))

	for indx in range(k_max):
		indx_t = np.where(true_cat.k[inside] == indx)[0]
		indx_f = np.where(fitted_cat.k == indx)[0]
	
		if (len(indx_t == 1) and len(indx_f == 1)):
			true_flux = np.append(true_flux, true_cat.flux[inside][indx_t])
			fitted_flux = np.append(fitted_flux, fitted_cat.flux[indx_f])
	
	rms = analysis.rms_error(fitted_flux, true_flux)
	result = 'Survey {0}:\t\tRMS = {1:0.7f} %'.format(path, rms)
	f.write(result + '\n')
	print(result)

print('Simulation Complete!')
