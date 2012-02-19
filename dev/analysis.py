# Rory Holmes
# Feb 2012

# This file contains the functions to analyze the fitted flat-field in the
# self-calibration simulations.

# Make Python 3 compatible
from __future__ import division, print_function

# Standard Modules
import numpy as np
import god
import pickle

# Custom Modules
import self_calibration

def badness(p, q, verbose=False):

    if verbose:
        print("Calculating best-in-basis badness of fitted flat-field...")

    nx, ny = p['ff_samples']
    dx = p['FoV'][0] / (nx - 1)
    dy = p['FoV'][1] / (ny - 1)
    x = np.arange(- p['FoV'][0] / 2 + dx / 2, p['FoV'][0] / 2, dx)
    y = np.arange(- p['FoV'][1] / 2 + dy / 2, p['FoV'][1] / 2, dy)
    X, Y = np.meshgrid(x, y)
    temp_x = np.reshape(X, -1)
    temp_y = np.reshape(Y, -1)
    our_ff = self_calibration.evaluate_flat_field(p, temp_x, temp_y, q)
    god_ff = god.flat_field(p, temp_x, temp_y)
    badness = 100. * np.sqrt(np.mean(((our_ff - god_ff) / god_ff) ** 2))

    if verbose:
        print("...done!")

    return badness


def best_in_basis(p, q, data_dir, verbose=False):

    if verbose:
        print("Calculating best-in-basis badness of fitted flat-field...")

    nx, ny = p['ff_samples']
    dx = p['FoV'][0] / (nx - 1)
    dy = p['FoV'][1] / (ny - 1)
    x = np.arange(- p['FoV'][0] / 2 + dx / 2, p['FoV'][0] / 2, dx)
    y = np.arange(- p['FoV'][1] / 2 + dy / 2, p['FoV'][1] / 2, dy)
    X, Y = np.meshgrid(x, y)
    temp_x = np.reshape(X, -1)
    temp_y = np.reshape(Y, -1)
    our_ff = self_calibration.evaluate_flat_field(p, temp_x, temp_y, q)
    pickle_dic = pickle.load(open(data_dir + '/bestfit_ff.p'))
    bestfit_ff_parameters = pickle_dic['fit_parameters']
    bestfit_ff = self_calibration.evaluate_flat_field(p, temp_x, temp_y, \
                                                bestfit_ff_parameters)
    badness = 100. * np.sqrt(
                        np.mean(((our_ff - bestfit_ff) / bestfit_ff) ** 2))

    if verbose:
        print("...done!")

    return badness


def rms_error(flux_estimate, true_flux, verbose=False):
    if verbose:
        print("Calculating RMS Error of fitted sources...")
    rms = 100 * np.sqrt(np.mean(((flux_estimate - true_flux) \
                                                        / true_flux) ** 2))
    print("...done!")
    return rms
