# Rory Holmes
# Feb 2012

# This file contains the master function that controls a single end-to-end
# self_calibration simulation.

# Make Python 3 compatible
from __future__ import division, print_function

# Standard Modules
import numpy as np
import os
import pickle as pickle

# Custom self-cal modules
import self_calibration
import god
import survey
import save_out


def run_sim(mdic):
    # Create output directories
    os.system('mkdir -p {0}'.format(mdic['data_dir']))
    os.system('rm -r {0}/*'.format(mdic['data_dir']))
    os.mkdir((mdic['data_dir'] + '/FF'))
    # Dump parameters
    pickle.dump(mdic['params'], open((mdic['data_dir'] + '/parameters.p'),
                                                                         "wb"))
    
    if 'bestfit_ff' in mdic:
        print("Using pre-existing best fit to God's flat-field.")
    else:
        save_out.bestfit_ff(mdic['params'], mdic['data_dir']) 

    # Only create a sky catalog if not already generated
    if 'sky_catalog' in mdic:
        print("Using pre-existing sky catalog.")
        sky_catalog = mdic['sky_catalog']
    else:
        sky_catalog = god.create_catalog(mdic['params'], mdic['data_dir'],
                                                plots=mdic['plotdata'])

    obs_catalog = survey.survey(mdic['params'], sky_catalog,
                                    mdic['survey_file'], mdic['data_dir'],
                                    plots=mdic['plotdata'],
                                    verbose=mdic['verbosemode'])
    # observation_catalog: *.size, *.pointing_ID, *.star_ID, *.flux,
    # *.invvar, *.x, *.y

    if mdic['plotdata']:
        save_out.invvar(obs_catalog, mdic['data_dir'],
                                         verbose=mdic['verbosemode'])

    sln = self_calibration.self_calibration(mdic['params'], obs_catalog,
                                            sky_catalog, mdic['data_dir'],
                                            plots=mdic['plotdata'],
                                            verbose=mdic['verbosemode'])
    np.savetxt((mdic['data_dir'] + "/result"), sln)
