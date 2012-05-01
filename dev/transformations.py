# Rory Holmes (MPIA) / David Hogg (NYU)
# 2011 - 2012

# This file contains the transformation functions for the
# self-calibration simulations.

# Make Python 3 compatible
from __future__ import division, print_function

# Standard Modules
import numpy as np


def mag2flux(mag):
    ''' This function converts AB magnitudes to fluxes

    Input
    -----
    mag         :   numpy array
        AB magnitudes

    Returns
    -------
    out         :   numpy array
        Fluxes
    '''

    return np.power(10, (- 0.4 * (mag - 22.5)))


def flux2mag(flux):
    ''' This functions converts fluxes to AB magnitudes

    Input
    ----------
    flux        :   numpy array
        Fluxes

    Returns
    -------
    out         :   numpy array
        AB magnitudes
    '''

    return 22.5 - 2.5 * np.log10(flux)


def test_mag2flux_flux2mag():
    mag = np.linspace(-50., 50., 1000)
    assert np.sum((mag - flux2mag(mag2flux(mag))) ** 2) < 1e-10
    flux = np.linspace(1e-5, 1e5, 1000)
    assert np.sum((flux - mag2flux(flux2mag(flux))) ** 2) < 1e-10


def sky2fp(alpha, beta, pointing, orientation):
    ''' This function transforms sky coordinates into focal plane coordinates.
    alpha and beta are angles from the telescope normal and are not RA, dec
    etc. The focal plane coordinates are measured in degrees.

    Input
    ----------
    alpha           :      numpy array
        Array of alpha angles in degrees
    beta            :      numpy array
        Array of beta angles in degrees
    pointing        :      numpy array
        Telescope pointing (alpha, beta) in degrees
    orientation     :      float
        Telescope orientation in degrees

    Returns
    -------
    x               :       numpy array
        focal plane x positions in degrees
    y               :       numpy array
        focal plane y positions in degrees
    '''

    theta = - orientation * np.pi / 180.  # telescope rotates NOT sky
    x = (alpha - pointing[0]) * np.cos(theta) \
            - (beta - pointing[1]) * np.sin(theta)
    y = (alpha - pointing[0]) * np.sin(theta) \
            + (beta - pointing[1]) * np.cos(theta)
    return x, y


def fp2sky(x, y, pointing, orientation):
    ''' This function transforms focal plane coordinates into sky coordinates.
    alpha and beta are angles from the telescope normal and are not RA, dec
    etc. The focal plane coordinates are measured in degrees.

    Input
    ----------
    x               :       numpy array
        focal plane x positions in degrees
    y               :       numpy array
        focal plane y positions in degrees
    pointing        :      numpy array
        Telescope pointing (alpha, beta) in degrees
    orientation     :      float
        Telescope orientation in degrees

    Returns
    -------
    alpha           :      numpy array
        Array of alpha angles in degrees
    beta            :      numpy array
        Array of beta angles in degrees
    '''

    theta = - orientation * np.pi / 180.  # telescope rotates NOT sky
    alpha = x * np.cos(theta) + y * np.sin(theta) + pointing[0]
    beta = -x * np.sin(theta) + y * np.cos(theta) + pointing[1]
    return alpha, beta
