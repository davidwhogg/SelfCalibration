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
import true_functions
import save_out
import survey


def run_sim(dic):
    
    if dic['data_dir']:
        if os.path.isdir(dic['data_dir']):  # Create output directories
            os.system('rm -r {0}'.format(dic['data_dir']))
        os.mkdir(dic['data_dir'])
        os.mkdir((dic['data_dir'] + '/FF'))
    
    np.random.seed(dic['seed'])

    # Only find best fit to the true flat-field if not already done
    if 'best_fit_params' in dic:
        if dic['verbose']:
            print("Using pre-existing best fit to the true flat-field.")
    else:
        dic['best_fit_params'] = analysis.best_fit_ff(dic['FoV'],
            dic['ff_samples'], dic['flat_field_order'], dic['stop_condition'],
            dic['max_iterations'], verbose=dic['verbose'])
    
    # Only create a sky catalog if not already generated
    if 'sky_catalog' in dic:
        if dic['verbose']:
            print("Using pre-existing sky catalog.")
        sky_catalog = dic['sky_catalog']
    else:
        if dic['verbose']:
            print("Generating sky catalog...")
        sky_catalog = true_functions.SourceCatalog(dic['sky_limits'],
                        dic['density_of_stars'], dic['m_min'], dic['m_max'],
                        dic['powerlaw_constants'])
        if dic['verbose']:        
            print('...{0} sources generated!'.format(sky_catalog.size))
    if dic['data_dir']:
        save_out.source_catalog(dic['data_dir'], sky_catalog, dic['verbose'])
        save_out.parameters(dic['data_dir'], dic, dic['verbose'])
    
    # Perform survey
    survey_catalog = survey.survey(dic['survey_file'], sky_catalog, dic['FoV'], dic['eta'], dic['delta'], dic['epsilon_max'], data_dir=False, verbose=dic['verbose'])
    if dic['data_dir']:
        save_out.survey(dic['data_dir'], survey_catalog, dic['verbose'])
