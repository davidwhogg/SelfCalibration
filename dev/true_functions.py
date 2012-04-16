# Rory Holmes
# Feb 2012

# This file contains the functions to create true things:
# the sky and the true flat-field.

# Make Python 3 compatible
from __future__ import division, print_function

# Standard Modules
import numpy as np


def flat_field_parameters():
    ''' This function just contains the parameters used to generate the *true*
    flat-field

    Return
    ------
    out     :   numpy array
        The parameters for the *true* flat-field model
    '''

    return np.append(np.array([1, -0.02, 0.03, -0.2, 0.01, -0.5]),
                                (1e-4) * np.random.uniform(size=256))


def flat_field(x, y, FoV, par=flat_field_parameters()):
    ''' This function returns the *true* flat-field model values at the given
    focal plane sample points

    Input
    -----
    x       :   numpy.array
        A numpy array with the x-coordinates of the sample points
    y       :   numpy.array
        A numpy array with the y-coordinates of the sample points
    FoV     :   Float array
        The simulate imager's field-of-view in degrees [alpha, beta]
    par     :   numpy.array
        The parameters for the *true* flat-field model

    Return
    ------
    out     :   numpy array
        The values of the *true* flat-field at the sample points provided
    '''

    assert x.shape == y.shape
    #  The large-scale flat-field, modeled as second order polynomial
    ff = (par[0] + par[1] * x + par[2] * y \
                    + par[3] * x ** 2 + par[4] * x * y + par[5] * y ** 2)
    # The small-scale flat-field, modeled with sine and cosine contributions
    k = 6
    no_loops = int(np.sqrt((len(par) - 6) / 4))
    for nx in range(no_loops):
        kx = nx * np.pi / FoV[0]
        ckx = np.cos(kx * x)
        skx = np.sin(kx * x)
        for ny in range(no_loops):
            ky = ny * np.pi / FoV[1]
            cky = np.cos(ky * y)
            sky = np.sin(ky * y)
            ff += par[k] * ckx * cky
            ff += par[k + 1] * ckx * sky
            ff += par[k + 2] * skx * cky
            ff += par[k + 3] * skx * sky
            k += 4
    return ff
