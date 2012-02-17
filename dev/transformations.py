# Rory Holmes
# Feb 2012

# This file contains the transformation functions for the
# self-calibration simulations. All data is saved out as pickles.

# Make Python 3 compatible
from __future__ import division, print_function

# Standard Modules
import numpy as np


def mag2flux(mag):
    ''' Converts magnitudes to fluxes

    Returns an numpy array in the same format containing fluxes

    Parameters
    ----------
    mag         :   np.array
        Array of magnitudes to convert to fluxes
    '''

    return np.power(10, (- 0.4 * (mag - 22.5)))


def flux2mag(flux):
    ''' Converts fluxes to magnitudes

    Returns an numpy array in the same format containing magnitudes

    Parameters
    ----------
    flux         :   np.array
        Array of fluxes to convert to magnitudes
    '''

    return 22.5 - 2.5 * np.log10(flux)


def test_mag2flux_flux2mag():
    mag = np.linspace(-50., 50., 1000)
    assert np.sum((mag - flux2mag(mag2flux(mag))) ** 2) < 1e-10
    flux = np.linspace(1e-5, 1e5, 1000)
    assert np.sum((flux - mag2flux(flux2mag(flux))) ** 2) < 1e-10


def sky2fp(alpha, beta, pointing, orientation):
    ''' Transforms the sky coordinates into focal plane coordinates

    Returns focal plane coordinates (x, y)

    Parameters
    ----------
    alpha         :   np.array
        Array of source alpha coordinates
    beta          :   np.array
        Array of source beta coordinates
    pointing      :   np.array()
        Telescope pointing (alpha, beta)
    orientation   :   float
        Telescope orientation
    '''

    theta = - orientation * np.pi / 180.  # telescope rotates NOT sky
    x = (alpha - pointing[0]) * np.cos(theta) \
            - (beta - pointing[1]) * np.sin(theta)
    y = (alpha - pointing[0]) * np.sin(theta) \
            + (beta - pointing[1]) * np.cos(theta)
    return x, y


def fp2sky(x, y, pointing, orientation):
    ''' Transforms the focal plane coordinates into sky coordinates

    Returns sky coordinates coordinates (alpha, beta)

    Parameters
    ----------
    x             :   np.array
        Array of source x focal plane coordinates
    y             :   np.array
        Array of source y focal plane coordinates
    pointing      :   np.array
        Telescope pointing (alpha, beta)
    orientation   :   float
        Telescope orientation
    '''
    theta = - orientation * np.pi / 180.  # telescope rotates NOT sky
    alpha = x * np.cos(theta) + y * np.sin(theta) + pointing[0]
    beta = -x * np.sin(theta) + y * np.cos(theta) + pointing[1]
    return alpha, beta
