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
from multiprocessing import Pool


# Custom Modules
import simulation

# Load the default parameters
dic = eval(open('parameters.py').read())

# Perform the four surveys (each a separate process)
survey_files = ['A', 'B', 'C', 'D']
mult_proc = []
for srvy in survey_files:
    dic['survey_file'] = srvy + '.txt'
    dic['data_dir'] = srvy
    mult_proc.append(dic)
p = Pool(4)
p.map(simulation.run_sim, mult_proc)

for srvy in survey_files:
    os.system('./../../plot.py {0}'.format(srvy))


