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
import self_calibration
import true
import survey
import save_out
import analysis


def run_sim(dic):
    if os.path.isdir(dic['data_dir']):  # Create output directories
        os.system('rm -r {0}'.format(dic['data_dir']))
    os.mkdir(dic['data_dir'])
    os.mkdir((dic['data_dir'] + '/FF'))

    save_out.parameters(dic)  # Dump parameters
    
    # Only find best fit to the true flat-field if not already done
    if 'best_fit_params' in dic:
        print("Using pre-existing best fit to the true flat-field.")
    else:
        dic['best_fit_params'] = analysis.best_fit_ff(dic)

    # Only create a sky catalog if not already generated
    if 'sky_catalog' in dic:
        print("Using pre-existing sky catalog.")
        sky_catalog = dic['sky_catalog']
    else:
        sky_catalog = true.create_catalog(dic)
    if dic["plotdata"]:
        save_out.source_catalog(dic, sky_catalog)
    
    obs_catalog = survey.survey(dic, sky_catalog,
                                    dic['survey_file'], dic['data_dir'])
    # observation_catalog: *.size, *.pointing_ID, *.star_ID, *.flux,
    # *.invvar, *.x, *.y

    if dic['plotdata']:
        save_out.invvar(dic, obs_catalog, dic['data_dir'])

    sln = self_calibration.self_calibration(dic, obs_catalog,
                                            sky_catalog, dic['data_dir'])
    np.savetxt((dic['data_dir'] + "/result"), sln)
