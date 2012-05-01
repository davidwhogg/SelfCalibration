# Rory Holmes (MPIA) / David Hogg (NYU)
# 2011 - 2012

# This file contains the high-level self-calibration function that performs
# a full simulation run (sky generation, survey and self-calibration) and
# returns its performance.

# Make Python 3 compatible
from __future__ import division, print_function

# Standard Modules
import numpy as np
import os

# Custom self-calibration modules
import analysis
import true_functions
import save_out
import survey
import self_calibration


def run_sim(dic):
    ''' This is the high level self-calibration function that actually runs
    one self-calibration run (sky generation, survey and self-calibration) and
    returns its performance. This function is included to allow simple Python
    multiprocessing, by calling this function with different parameter
    dictionaries in parallel.

    Input
    -----
    dic         :       dictionary
        A Python dictionary containing all the simulation parameters. See
        default_parameters.py for all the required parameters.

    Return
    ------
    out         :       numpy array
        A list containing the results of the self-calibration simulations
        list[0] = the number of self-calibration iterations required to
        converge to the final fitted solution (if 0 returned, then
        self-calibration did not converge).
        list[1] = the root-mean-squared error between the fitted source fluxes
        and the true source fluxes.
        list[2] = the "badness" between the fitted instrument response and
        the true instrument responses.
        list[3] = the "best-in-badness" between the fitted instrument
        response and the true instrument responses.
        list[4] = the chi2 of the final fitted solution
    '''

    # Create output directories
    if dic['data_dir']:
        if dic['data_dir'][-1] == '/':
            dic['data_dir'] = dic['data_dir'][0:-1]
        if os.path.isdir(dic['data_dir']):
            os.system('rm -r {0}'.format(dic['data_dir']))
        os.mkdir(dic['data_dir'])
        os.mkdir((dic['data_dir'] + '/FF'))

    # repeated runs should return consistent results
    np.random.seed(dic['seed'])

    # Only find best fit to the true flat-field if not already done
    if 'best_fit_params' in dic:
        if dic['verbose']:
            print("Using pre-existing best fit to the true flat-field.")
    else:
        dic['best_fit_params'] = analysis.best_fit_ff(
                                    dic['FoV'],
                                    dic['ff_samples'],
                                    dic['flat_field_order'],
                                    dic['stop_condition'],
                                    dic['max_iterations'],
                                    verbose=dic['verbose'])

    # Only create a sky catalog if one not provided in parameter dictionary
    if 'sky_catalog' in dic:
        if dic['verbose']:
            print("Using pre-existing sky catalog.")
        sky_catalog = dic['sky_catalog']
    else:
        if dic['verbose']:
            print("Generating sky catalog...")
        sky_catalog = true_functions.SourceCatalog(dic['sky_limits'],
                                                    dic['density_of_stars'],
                                                    dic['m_min'],
                                                    dic['m_max'],
                                                    dic['powerlaw_constants'])
        if dic['verbose']:
            print('...{0} sources generated!'.format(sky_catalog.size))
    if dic['data_dir']:
        save_out.source_catalog(dic['data_dir'], sky_catalog, dic['verbose'])
        save_out.parameters(dic['data_dir'], dic, dic['verbose'])

    # Perform sky survey
    measurement_catalog = survey.survey(dic['survey_file'],
                                                sky_catalog,
                                                dic['FoV'],
                                                dic['eta'],
                                                dic['delta'],
                                                dic['epsilon_max'],
                                                data_dir=dic['data_dir'],
                                                verbose=dic['verbose'])
    if dic['data_dir']:
        save_out.measurement_catalog(dic['data_dir'], measurement_catalog,
                                                                dic['verbose'])

    performance = self_calibration.self_calibration(measurement_catalog,
                                                    sky_catalog,
                                                    dic['flat_field_order'],
                                                    dic['FoV'],
                                                    dic['ff_samples'],
                                                    dic['stop_condition'],
                                                    dic['max_iterations'],
                                                    dic['best_fit_params'],
                                                    data_dir=dic['data_dir'],
                                                    verbose=dic['verbose'])

    return performance
