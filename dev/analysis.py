# Rory Holmes
# Feb 2012

# This file contains the functions to analyze the fitted flat-field in the
# self-calibration simulations.

# Make Python 3 compatible
from __future__ import division, print_function

# Standard Modules
import numpy as np
import scipy.optimize as opt

# Custom Modules
import self_calibration as self_cal
import god


def badness(p, q):
    ''' Calculates the 'badness' of the instrument response fit. The
    'badness' is defined as the root-mean-square difference between
    the true instrument response and the fitted instrument response
    measured on a regular grid of sample points across the focal plane.

    Parameters
    ----------
    p          :   dictionary
        the dictionary containing all the simulation parameters
    q          :   numpy array
        the fitted instrument response parameters

    Returns
    -------
    out        :   float
        the calculated badness
    '''

    if p['verbose']:
        print("Calculating best-in-basis badness of fitted flat-field...")

    x = np.linspace(- p['FoV'][0] / 2, p['FoV'][0] / 2, p['ff_samples'][0])
    y = np.linspace(- p['FoV'][1] / 2, p['FoV'][1] / 2, p['ff_samples'][1])
    X, Y = np.meshgrid(x, y)

    our_ff = self_cal.evaluate_flat_field(p, x.flatten(), y.flatten(), q)
    god_ff = god.flat_field(p, x.flatten(), y.flatten())

    badness = 100. * np.sqrt(np.mean(((our_ff - god_ff) / god_ff) ** 2))

    if p['verbose']:
        print("...done!")

    return badness


def best_in_basis(p, q):
    ''' Calculates the 'best-in-basis badness' of the instrument response fit.
    This is defined as the root-mean-square difference between the fitted
    instrument response and the *best fit possible* with the model basis used
    for the fit (on a regular grid of sample points across the focal plane).

    Parameters
    ----------
    p          :   dictionary
        the dictionary containing all the simulation parameters
    q          :   numpy array
        the fitted instrument response parameters

    Returns
    -------
    out        :   float
        the calculated best-in-basis badness
    '''
    if p['verbose']:
        print("Calculating best-in-basis badness of fitted flat-field...")

    x = np.linspace(- p['FoV'][0] / 2, p['FoV'][0] / 2, p['ff_samples'][0])
    y = np.linspace(- p['FoV'][1] / 2, p['FoV'][1] / 2, p['ff_samples'][1])
    X, Y = np.meshgrid(x, y)

    our_ff = self_cal.evaluate_flat_field(p, x.flatten(), y.flatten(), q)
    bf_ff = self_cal.evaluate_flat_field(
                            p, x.flatten(), y.flatten(), p['best_fit_params'])

    best_in_basis_badness = \
                    100. * np.sqrt(np.mean(((our_ff - bf_ff) / bf_ff) ** 2))

    if p['verbose']:
        print("...done!")

    return best_in_basis_badness


def rms_error(p, flux, true_flux):
    ''' Calculates the root-mean-square error between the estimated source
    fluxes and the true source fluxes. This only applies to the bright sources
    selected for the self-calibration procedure and not necessarily all the
    sources on the sky.

    Parameters
    ----------
    p           :   dictionary
        the dictionary containing all the simulation parameters
    flux        :   numpy array
        containing all the source flux estimates
    true_flux   :   numpy array
        containing all the true fluxes for the sources

    Returns
    -------
    out        :   float
        the calculated root-mean-square error
    '''

    if p['verbose']:
        print("Calculating RMS Error of fitted sources...")

    rms = 100 * np.sqrt(np.mean(((flux - true_flux) / true_flux) ** 2))

    if p['verbose']:
        print("...done!")

    return rms


def best_fit_ff(p):
    ''' Calculates the best fit possible to true flat-field with the basis
    used to model it during the self-calibration procedure.

    Parameters
    ----------
    p       :   dictionary
        simulation parameters

    Returns
    -------
    out     :   numpy array
        the fitted parameter
    '''

    if p['verbose']:
        print("Fitting god's flat-field with basis...")

    order = p['flat_field_order']
    x = np.linspace(- p['FoV'][0] / 2, p['FoV'][0] / 2, p['ff_samples'][0])
    y = np.linspace(- p['FoV'][1] / 2, p['FoV'][1] / 2, p['ff_samples'][1])
    X, Y = np.meshgrid(x, y)

    g = self_cal.evaluate_flat_field_functions(x.flatten(), y.flatten(), order)
    god_ff = god.flat_field(p, x.flatten(), y.flatten())
    a = np.zeros((order + 1) * (order + 2) / 2)
    a[0] = 1  # Start roughly near the correct solution
    a[3] = -0.2
    a[5] = 0.5
    fitted_parameters = opt.fmin_bfgs(self_cal.compare_flats, a,
                        args=(god_ff, g, x.flatten(), y.flatten()), \
                        gtol=p['stop_condition'], maxiter=p['max_iterations'])
    fitted_parameters = self_cal.normalize_flat_field(p, fitted_parameters)

    if p['verbose']:
        print("...done!")
        print("Fitted Parameters: ", fitted_parameters)

    return fitted_parameters
