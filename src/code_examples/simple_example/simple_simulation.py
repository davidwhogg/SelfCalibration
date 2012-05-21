#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Rory Holmes (MPIA) / David Hogg (NYU)
# 2012
# This script performs a simple self-calibration simulation with the parameters
# defined in 'parameters.py'. The output is saved to the default_output/
# directory and then used to generate plots.

# Make Python 3 compatible
from __future__ import division, print_function

# Standard Modules
import os
import sys

# Add simulator modules to Python path
sys.path.append('./../..')  

# Custom Modules
import simulation


dic = eval(open('parameters.py').read())

results = simulation.run_sim(dic)

if dic['data_dir']:
    os.system('./../../plot.py {0}'.format(dic['data_dir']))

print('Final Fitted Solution:')
print('======================')
print('Self-Calibration Iterations to Converge: {0:.5}'.format(results[0]))
print('RMS Source Error: {0:.5} %'.format(results[1]))
print('Instrument Response Badness: {0:.5} %'.format(results[2]))
print('Instrument Response Best-in-Basis Badness: {0:.5} %'.format(results[3]))
print('Chi2 of Fit: {0:.5}'.format(results[4]))
