# Rory Holmes
# Feb 2012

# This file contains the functions to take a single image of the synthetic
# sky in the self-calibration simulations.

# Make Python 3 compatible
from __future__ import division, print_function

# Standard Modules
import numpy as np
import os

# Custom self-cal modules
import transformations as tran
import god
import save_out


class CameraCatalog:
    def __init__(self, sky_catalog, pointing, orientation):
        self.k = sky_catalog.k.astype(int)
        self.mag = sky_catalog.mag
        self.x, self.y = tran.sky2fp(sky_catalog.alpha, sky_catalog.beta,\
                                    pointing, orientation)
        self.size = sky_catalog.size
        self.flux = tran.mag2flux(self.mag)
        self.epsilon = sky_catalog.epsilon


class MeasuredCatalog:
    def __init__(self, params, camera_catalog, inside_FoV):
        self.size = len(inside_FoV[0])
        self.k = camera_catalog.k[inside_FoV].astype(int)
        self.x = camera_catalog.x[inside_FoV]
        self.y = camera_catalog.y[inside_FoV]
        flat = god.flat_field(params, self.x, self.y)
        self.gods_invvar = self.true_invvar(params, camera_catalog,\
                                                inside_FoV, flat)
        self.counts = camera_catalog.flux[inside_FoV] * flat \
                        + np.random.normal(size=self.size) \
                        / np.sqrt(self.gods_invvar)
        self.invvar = self.reported_invvar(params)

    def append(self, other):
        self.size = self.size + other.size
        self.k = np.append(self.k, other.k)
        self.x = np.append(self.x, other.x)
        self.y = np.append(self.y, other.y)
        self.counts = np.append(self.counts, other.counts)
        self.invvar = np.append(self.invvar, other.invvar)
        self.gods_invvar = np.append(self.gods_invvar, other.gods_invvar)

    def mag(self):
        return tran.flux2mag(self.counts)

    def reported_invvar(self, params):
        sky_unc = params['delta']
        var = (sky_unc ** 2 + (params['eta'] ** 2) * self.counts ** 2)
        return 1. / var

    def true_invvar(self, params, camera_catalog, inside_FoV, flat):
        epsilon = camera_catalog.epsilon[inside_FoV]
        sky_unc = params['delta']
        flux = camera_catalog.flux[inside_FoV]
        var = (sky_unc ** 2 * (1. + epsilon ** 2) + \
                (params['eta'] ** 2) * flat ** 2 * flux ** 2)
        return 1. / var


def measure(p, sky_catalog, pointing, orientation, \
                data_dir, plots=None, verbose=None):
    if verbose:
        print("Converting sky catalog to focal plane coordinates...")
    camera_catalog = CameraCatalog(sky_catalog, pointing, orientation)
    if verbose:
        print("...done!")

    if verbose:
        print("Finding stars within camera FoV...")
    x_min = - 0.5 * p['FoV'][0]
    y_min = - 0.5 * p['FoV'][1]
    x_max = 0.5 * p['FoV'][0]
    y_max = 0.5 * p['FoV'][1]
    inside_FoV = np.where((x_min < camera_catalog.x) \
                    & (camera_catalog.x < x_max) \
                    & (y_min < camera_catalog.y) \
                    & (camera_catalog.y < y_max))
    if verbose:
        print("...done!")

    if verbose:
        print("Measuring stars within FoV...")
    measured_catalog = MeasuredCatalog(p, camera_catalog, inside_FoV)
    if verbose:
        print("...done!")

    one_camera_file = os.path.exists((data_dir + '/camera_image.p'))
    if plots and (one_camera_file != True) and (orientation > 30) \
    and (pointing[0] > -1) and (pointing[0] < 1) and (pointing[1] > -1) \
    and (pointing[1] < 1):
        save_out.camera(p, sky_catalog, measured_catalog, inside_FoV,
                            pointing, orientation, data_dir, verbose=verbose)
    return measured_catalog
    # measured_sources  *.size, *.k, *.flux, *.invvar, *.x, *.y
