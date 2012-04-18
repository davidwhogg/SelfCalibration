#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Make Python 3 compatible
from __future__ import division, print_function

import simulation
import os

dic = eval(open("default_parameters.py").read())

dic['survey_file'] = "survey"
dic['verbose'] = True   

dic['best_fit_params'] = False
result = simulation.run_sim(dic)

os.system('./plot.py {0}'.format(dic['data_dir']))
