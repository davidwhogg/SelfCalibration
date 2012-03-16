# Rory Holmes
# Feb 2012

# This file contains the functions to create the things that God only knows:
# the sky and the true flat-field.

# Make Python 3 compatible
from __future__ import division, print_function

# Standard Modules
import numpy as np
import scipy.optimize as opt
import sys
import pickle

# Custom self-cal modules
import transformations as tran


def flat_field_parameters():
    #np.random.seed(1)  # same seed => same FF
    return np.append(np.array([1, -0.02, 0.03, -0.2, 0.01, -0.5]),
                                (1e-4) * np.random.uniform(size=256))


def flat_field(p, x, y, par=flat_field_parameters()):
    # Large Scale
    ff = (par[0] + par[1] * x + par[2] * y \
                    + par[3] * x ** 2 + par[4] * x * y + par[5] * y ** 2)
    # Small Scale
    FoV = p['FoV']
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


def power_law(A, m):
    temp = np.zeros(len(m))
    for indx in range(len(A)):
        temp += A[indx] * m ** indx
    return 10 ** temp


def error(a, d, m):
    error = np.sum(((d - power_law(a, m)) / d) ** 2)
    return error


def prob_dist(m_min, m_max, a, size):
    # Randomly generate source magnitudes according to probability distribution
    r = np.random.uniform(size=size)
    m = (1.0 / a) * np.log10((10. ** (a * m_max) - 10. ** (a * m_min)) \
                        * r + 10. ** (a * m_min))
    return m


def generate_magnitudes(p, number_sources, data_dir, plots=False):
    A = p['powerlaw_constants']
    area = (p['sky_limits'][1] - p['sky_limits'][0]) \
                        * (p['sky_limits'][3] - p['sky_limits'][2])
    # fit for dN/dm = B[0] + B[1] * m
    mag_range = np.arange(p['m_min'] + 0.25, p['m_max'], 0.5)
    B = opt.fmin(error, np.array([0, 0]), \
                    args=(power_law(A, mag_range), mag_range), \
                    xtol=1e-14, maxiter=1e16, maxfun=1e16)
    B = np.append(B, [0])
    # Shift above dN/dm = a + b*m + c*m**2 line
    chck = 1.
    while chck > 0.1:
        bad = power_law(A, mag_range) - power_law(B, mag_range)
        ii = np.where(bad > 0.)
        B[0] += 0.001
        chck = np.sum(bad[ii[0]])

        total_sources = np.sum(power_law(A, mag_range)) * area
    if total_sources * p['useful_fraction'] < number_sources:
        print("Error! You want more sources than there are in the sky...")
        sys.exit()

    mag = prob_dist(p['m_min'], p['m_max'], B[1], p['useful_fraction'] * \
                                                            total_sources)

    # Use rejection method to get correct probability distribution
    c = np.max(power_law(A, mag) / power_law(B, mag))
    difference = power_law(A, mag) / (power_law(B, mag) * c)
    randoms = np.random.uniform(size=len(difference))
    ii = np.zeros((2, 2))
    while len(ii[0]) > 1:
        ii = np.where(randoms > difference)
        mag[ii] = prob_dist(p['m_min'], p['m_max'], B[1], len(ii[0]))
        randoms = np.random.uniform(size=len(ii[0]))
        difference = power_law(A, mag[ii]) / (power_law(B, mag[ii]) * c)

    selected_sources = np.sort(mag)[0:int(number_sources)]
    return selected_sources


class SourceCatalog:
    def __init__(self, params, data_dir, number_sources, plots=False):
        np.random.seed(params['seed'])
        self.k = np.arange(number_sources).astype(int)
        self.alpha = np.random.uniform(low=params['sky_limits'][0], \
                        high=params['sky_limits'][1], size=number_sources)
        self.beta = np.random.uniform(low=params['sky_limits'][2], \
                        high=params['sky_limits'][3], size=number_sources)
        self.mag = generate_magnitudes(params, number_sources, \
                                                    data_dir, plots=plots)
        self.size = number_sources
        self.flux = tran.mag2flux(self.mag)
        self.epsilon = np.random.uniform(0, params['epsilon_max'], \
                                                        size=number_sources)


def create_catalog(p):
    if p["verbose"]:
        print("Generating God's Catalog...")
    # Calculate total number of stars in catalog
    number_sources = p['density_of_stars'] \
                    * (p['sky_limits'][1] - p['sky_limits'][0]) \
                    * (p['sky_limits'][3] - p['sky_limits'][2])
    # Create catalog
    catalog = SourceCatalog(p, p["data_dir"], number_sources, plots=p["plotdata"])
    if p["verbose"]:
        print("...done!")
    return catalog
    # catalog.ID, catalog.mag, catalog.alpha, catalog.beta, catalog.size
