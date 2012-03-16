# Rory Holmes
# Feb 2012

# This file contains the functions to analyze the fitted flat-field in the
# self-calibration simulations.

# Make Python 3 compatible
from __future__ import division, print_function

# Standard Modules
import numpy as np
import god
import scipy.optimize as opt

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
    bestfit_ff_parameters = p['best_fit_params']
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


def best_fit_ff(p, verbose=False):
    ''' Calculates the best fit possible to God's flat-field with the basis
    used to model it during the self-calibration procedure.

    Returns a vector of the best fit parameters

    p                   :   dictionary
        simulation parameters
    '''

    if verbose:
        print("Fitting god's flat-field with basis...")

    order = p['flat_field_order']
    ff_samples = p['ff_samples']
    temp_x = np.linspace(-0.5 * p['FoV'][0], 0.5 * p['FoV'][0], ff_samples[0])
    temp_y = np.linspace(-0.5 * p['FoV'][1], 0.5 * p['FoV'][1], ff_samples[1])
    X, Y = np.meshgrid(temp_x, temp_y)
    x = np.reshape(X, -1)
    y = np.reshape(Y, -1)
    g = self_calibration.evaluate_flat_field_functions(x, y, order)
    god_ff = god.flat_field(p, x, y)
    a = np.zeros((order + 1) * (order + 2) / 2)
    a[0] = 1  # Start roughly near the correct solution
    a[3] = -0.2
    a[5] = 0.5
    fitted_parameters = opt.fmin_bfgs(self_calibration.compare_flats, a,
                        args=(god_ff, g, x, y), \
                        gtol=p['stop_condition'], maxiter=p['max_iterations'])
    fitted_parameters = self_calibration.normalize_flat_field(p,
                                                        fitted_parameters)

    if verbose:
        print("...done!")
        print("Fitted Parameters: ", fitted_parameters)

    return fitted_parameters
