#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Make Python 3 compatible
from __future__ import division, print_function

import simulation
import os
import numpy as np

dic = eval(open("default_parameters.py").read())

result = simulation.run_sim(dic)

print(result)

if dic['data_dir']:
    os.system('./plot.py {0}'.format(dic['data_dir']))
