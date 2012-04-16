#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Make Python 3 compatible
from __future__ import division, print_function

import simulation
import os

dic = eval(open("default_parameters.py").read())

dic['survey_file'] = "survey"
dic['verbose'] = True   
result = simulation.run_sim(dic)
