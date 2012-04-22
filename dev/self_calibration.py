# Rory Holmes
# Feb 2012

# This file contains the functions to self-calibrate the photometric catalog
# in the self-calibration simulations.

# Make Python 3 compatible
from __future__ import division, print_function

# Standard Modules
import numpy as np

# Custom Modules
import true_functions


def evaluate_flat_field_functions(x, y, order):
    ''' This functions returns the g values for the flat-field calculation.

    Parameters
    ----------
    x               :   numpy.array
        The x-coordinates of the sample points
    y               :   numpy.array
        The y-coordinates of the sample points
    order           :   Integer
        The order of the polynomial flat-field used to fit the instrument
        response in the self-calibration procedure

    Returns
    -------
    out             :   numpy array
        the g parameter
    '''

    L = (order + 1) * (order + 2) / 2
    g = np.zeros((len(x), L))
    l = 0
    for n in range(order + 1):
        for m in range(n + 1):
            g[:, l] = (x ** (n - m)) * (y ** m)
            l += 1
    return g


def evaluate_flat_field(x, y, q):
    ''' This functions returns the flat-field values

    Parameters
    ----------
    x               :   numpy.array
        The x-coordinates of the sample points
    y               :   numpy.array
        The y-coordinates of the sample points
    q               :   numpy.array
        The flat-field parameters

    Returns
    -------
    out             :   numpy array
        the flat-field values at the given sample points
    '''

    assert x.shape == y.shape
    order = int(np.around(np.sqrt(0.25 + 2 * len(q)) - 1.5))
    assert(len(q) == ((order + 1) * (order + 2) / 2))
    g = evaluate_flat_field_functions(x, y, order)
    return np.dot(g, q)


def normalize_flat_field(q, FoV, ff_samples):
    ''' This functions modifies the flat-field parameters (q) in such a way
    to normalize the fitted flat-field compared to the *true* flat-field. This
    is necessary as the self-calibration procedure only constrains the
    relative calibration, not the absolute (it is not possible to know if
    all the sources are 50% fainter, or if the instrument response is 50%
    lower). In reality this degeneracy will be solved with absolute standards.

    Parameters
    ----------
    q               :   numpy.array
        The flat-field parameters

    Returns
    -------
    out             :   numpy array
        the modified flat-field parameters, which describe a normalized
        flat-field
    '''
    fit_mean = average_over_ff(FoV, ff_samples, evaluate_flat_field, (q))
    true_mean = average_over_ff(FoV, ff_samples, true_functions.flat_field,
                                    (true_functions.flat_field_parameters()))
    return (q * true_mean / fit_mean)


def average_over_ff(FoV, ff_samples, func, ff_params):
    ''' This functions returns the average of the flat-field values at the
    sample points given.

    Parameters
    ----------
    FoV             :   float array
        The simulate imager's field-of-view in degrees [alpha, beta]
    ff_samples      :   integer array
        The number of sample points on the focal plane for the best-in-basis
        fitting and all of the badness measures
    func            :   function
        The function describing the instrument flat-field
    ff_params       :   int array
        The flat-field parameters describing flat-field

    Returns
    -------
    out     :   numpy array
        the flat-field values at the given sample points
    '''

    x = np.linspace(-0.5 * FoV[0], 0.5 * FoV[0], ff_samples[0])
    y = np.linspace(-0.5 * FoV[1], 0.5 * FoV[1], ff_samples[1])
    X, Y = np.meshgrid(x, y)
    return np.mean(func(X.flatten(), Y.flatten(), ff_params))
