# Rory Holmes
# Feb 2012

# This file contains the master function that controls a single end-to-end
# self_calibration simulation.

# Make Python 3 compatible
from __future__ import division, print_function

# Standard Modules
import numpy as np
import os

# Custom self-cal modules
import analysis


def run_sim(dic):
    if os.path.isdir(dic['data_dir']):  # Create output directories
        os.system('rm -r {0}'.format(dic['data_dir']))
    os.mkdir(dic['data_dir'])
    os.mkdir((dic['data_dir'] + '/FF'))

    # Only find best fit to the true flat-field if not already done
    if 'best_fit_params' in dic:
        print("Using pre-existing best fit to the true flat-field.")
    else:
        dic['best_fit_params'] = analysis.best_fit_ff(dic['FoV'],
            dic['ff_samples'], dic['flat_field_order'], dic['stop_condition'],
            dic['max_iterations'], verbose=dic['verbose'])
    
    print(dic)
