# Rory Holmes
# Feb 2012

# This file contains the functions to perform a survey over the synthetic
# sky in the self-calibration simulations.

# Make Python 3 compatible
from __future__ import division, print_function

# Standard Modules
import numpy as np

# Custom self-cal modules
import single_image


def load_survey(filename, verbose=False):

    if verbose:
        print("Loading survey from {0}...".format(filename))
    pointings = np.loadtxt(filename)
    if verbose:
        print("...done!")

    return pointings


def survey(params, sky_catalog, survey_file, data_dir, plots=False, \
                verbose=False):

    pointing = load_survey(survey_file)
    number_pointings = len(pointing[:, 0])

    if verbose:
        print("Surveying sky...")

    obs_cat = None
    for i in range(number_pointings):
        si = single_image.single_image(params, sky_catalog, \
                    np.array([pointing[i, 1], pointing[i, 2]]), \
                    pointing[i, 3], data_dir, plots=plots, verbose=verbose)
        if obs_cat is None:
            obs_cat = si
        else:
            obs_cat.append(si)

    if verbose:
        print("...done!")

    return obs_cat
