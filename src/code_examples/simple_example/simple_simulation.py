#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Rory Holmes (MPIA) / David Hogg (NYU)
# 2012

# Make Python 3 compatible
from __future__ import division, print_function

# Standard Modules
import os
import sys
sys.path.append('./../..')  # add simulator modules to Python path

# Custom Modules
import simulation

dic = eval(open('parameters.py').read())

print(simulation.run_sim(dic))

if dic['data_dir']:
    os.system('./../../plot.py {0}'.format(dic['data_dir']))
