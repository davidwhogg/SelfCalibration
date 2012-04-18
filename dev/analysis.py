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
import true_functions


def best_fit_ff(FoV, ff_samples, order, stop_condition, max_iterations,
                            verbose=False):
    ''' Calculates the best fit possible to true flat-field with the basis
    used to model it during the self-calibration procedure.

    Parameters
    ----------
    FoV             :   Float Array
        The simulate imager's field-of-view in degrees [alpha, beta]
    ff_samples      :   Integer Array
        The number of sample points on the focal plane for the best-in-basis
        fitting and all of the badness measures
    order           :   Integer
        The order of the polynomial flat-field used to fit the instrument
        response in the self-calibration procedure
    stop_condition  :   Float
        The stop condition for the self-calibration procedure and the
        best-in-basis fitting (stop when difference is less than 2 times this)
    max_iterations  :   Integer
        The maximum number of iterations in the self-calibration procedure and
        the best-in-basis fitting
    verbose         :   Boolean
        Set to true to run the simulations in verbose mode

    Returns
    -------
    out     :   numpy array
        the fitted parameter
    '''

    if verbose:
        print("Fitting the true flat-field with basis...")

    x = np.linspace(- FoV[0] / 2, FoV[0] / 2, ff_samples[0])
    y = np.linspace(- FoV[1] / 2, FoV[1] / 2, ff_samples[1])
    X, Y = np.meshgrid(x, y)

    g = self_cal.evaluate_flat_field_functions(X.flatten(), Y.flatten(), order)
    true_ff = true_functions.flat_field(X.flatten(), Y.flatten(), FoV)
    a = np.zeros((order + 1) * (order + 2) / 2)
    fitted_parameters = opt.fmin_bfgs(compare_flats, a,
                        args=(true_ff, g),
                        gtol=stop_condition,
                        maxiter=max_iterations,
                        disp=verbose)
    fitted_parameters = self_cal.normalize_flat_field(fitted_parameters, FoV,
                            ff_samples)
    if verbose:
        print("...done!")
        print("Fitted Parameters: ", fitted_parameters)

    return fitted_parameters


def compare_flats(a, true_ff, g):
    ''' This is the error function used in the best-in-basis flat-field
    fitting function.

    Parameters
    ----------
    a             :     np.array
        The simulate imager's field-of-view in degrees [alpha, beta]
    ff_samples      :   Float Array
        The number of sample points on the focal plane
    func            :   function
        The function describing the flat-field model
    ff_params       :   np.array
        The parameters describing the flat-field model

    Returns
    -------
    out     :   numpy.array
        the mean of the flat-field at the given sample points
    '''
    error = np.sum((true_ff - np.dot(g, a)) ** 2)
    return error
