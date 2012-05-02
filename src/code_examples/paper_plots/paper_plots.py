#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Rory Holmes (MPIA) / David Hogg (NYU)
# 2012

# This code reproduces the plots from the paper

# Make Python 3 compatible
from __future__ import division, print_function

# Standard Modules
import os
import sys
sys.path.append('./../..')  # add simulator modules to Python path

# Custom Modules
import simulation

# Load the default parameters
dic = eval(open('parameters.py').read())

# Perform the four surveys
survey_files = ['A', 'B', 'C', 'D']
for srvy in survey_files:
    dic['survey_file'] = srvy + '.txt'
    dic['data_dir'] = srvy
    simulation.run_sim(dic)

if dic['data_dir']:
    os.system('./../../plot.py {0}'.format(dic['data_dir']))


