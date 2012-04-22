# Rory Holmes
# Feb 2012

# This file contains the functions to create true things:
# the sky and the true flat-field.

# Make Python 3 compatible
from __future__ import division, print_function

# Standard Modules
import numpy as np
from scipy.integrate import quad
from scipy.optimize import fmin
import sys

# Custom self-calibration modules
import transformations


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
    x       :   numpy array
        A numpy array with the x-coordinates of the sample points
    y       :   numpy array
        A numpy array with the y-coordinates of the sample points
    FoV     :   float array
        The simulate imager's field-of-view in degrees [alpha, beta]
    par     :   numpy array
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


class SourceCatalog:
    def __init__(self, sky_limits, density_of_stars, m_min, m_max, A):
        ''' This class generates the Source Catalog object

        Input
        -----
        sky_limits          :   float array
            The area of sky to generate sources in
            [alpha_min, alpha_max, beta_min, beta_max]
        density_of_stars    :   int
            The maximum number of sources (all magnitude) per unit area
            to generate for the self-calibration simulations
        m_min               :   float
            The saturation limit of the simulated imager
        m_max               :   float
            The 10-sigma detection limit of the simulated imager
        A                   :   numpy array
            The parameters describing the magnitude distribution of the sources
            in the sky, according to: log10(dN/dm) = A + B * mag + C * mag ** 2

        '''

        area = (sky_limits[1] - sky_limits[0]) \
                                    * (sky_limits[3] - sky_limits[2])
        self.size = int(density_of_stars * area)
        self.k = np.arange(self.size).astype(int)
        self.alpha = np.random.uniform(low=sky_limits[0], \
                        high=sky_limits[1], size=self.size)
        self.beta = np.random.uniform(low=sky_limits[2], \
                        high=sky_limits[3], size=self.size)
        self.mag = self.generate_magnitudes(m_min, m_max, A, area, self.size)
        self.flux = transformations.mag2flux(self.mag)

    def power_law(self, A, m):
        ''' This function calculates the value of the source magnitude
        distribution power law.

        Input
        -----
        A           :   numpy.array / float
            The power law parameter(s)
        m           :   numpy.array / float
            The magnitude(s) to calculate the power law at

        Return
        ------
        out         :   numpy array / float
            The value of the power law
        '''
        power = 0 * m
        for indx in range(len(A)):
            power += A[indx] * m ** indx
        return 10 ** power

    def error(self, a, d, m):
        ''' This function calculates the error between the 10 ** [2nd order]
        and the 10 ** [1st order] function. It is used to fit the latter
        to the given power law constants.

        Input
        -----
        a           :   numpy.array
            The 10 ** [1st order] power law functions being fitted
        d           :   numpy.array
            The values of the actual power law at the magnitude sample points
        m           :   numpy.array
            The sample magnitude points

        Return
        ------
        out         :   numpy array / float
            The error to minimize (not RMS to reduce computational load)
        '''
        error = np.sum(((d - self.power_law(a, m)) / d) ** 2)
        return error

    def prob_dist(self, m_min, m_max, a, size):
        ''' This function calculates the 10 ** [1st order] probability
        distribution and is used in the rejection method.

        Input
        -----
        m_min               :   float
            The saturation limit of the simulated imager
        m_max               :   float
            The 10-sigma detection limit of the simulated imager
        a                   :   numpy.array
            The 10 ** [1st order] power law functions being fitted
        size                :   int
            The number of magnitudes to generate

        Return
        ------
        out         :   numpy array
            The random magnitudes generated
        '''
        r = np.random.uniform(size=size)
        m = (1.0 / a) * np.log10((10. ** (a * m_max) - 10. ** (a * m_min)) \
                            * r + 10. ** (a * m_min))
        return m

    def generate_magnitudes(self, m_min, m_max, A, area, size):
        ''' This function uses the rejection method to generate sources with
        a magnitude distribution defined by the given powerlaw parameters. The
        function returns only the brightest sources with the sky, up to the
        source density selected.

        Input
        -----
        m_min               :   float
            The saturation limit of the simulated imager
        m_max               :   float
            The 10-sigma detection limit of the simulated imager
        A                   :   numpy array
            The parameters describing the magnitude distribution of the sources
            in the sky, according to: log10(dN/dm) = A + B * mag + C * mag ** 2
        area                :   float
            The area of the sky patch being investigated
        size                :   int
            The number of sources to generate random magnitudes for

        Return
        ------
        out     :   numpy array
            The magnitude of the sources
        '''

        # For the rejection method we need to have a linear line, just above
        # the true line => fit for dN/dm = B[0] + B[1] * m
        mag_range = np.linspace(m_min, m_max, 20)
        B = fmin(self.error, np.zeros(2), \
                        args=(self.power_law(A, mag_range), mag_range), \
                        xtol=1e-14, maxiter=1e16, maxfun=1e16, disp=False)
        B = np.append(B, [0])
        # Shift above dN/dm = A[0] + A[1] * m + A[2] * m ** 2 line
        chck = 1.
        while chck > 0.1:
            bad = self.power_law(A, mag_range) - self.power_law(B, mag_range)
            ii = np.where(bad > 0.)
            B[0] += 0.01
            chck = np.sum(bad[ii[0]])

        # Calculate the total number of sources per area from the source's
        # magnitude distribution
        func = lambda x, a: 10 ** (a[0] + a[1] * x + a[2] * x ** 2)
        total_sources = int(area * quad(func, m_min, m_max, args=(A,))[0])
        if total_sources < size:
            print("Error! You want more sources than there are in the sky...")
            sys.exit()

        # Rejection method
        mag = self.prob_dist(m_min, m_max, B[1], total_sources)
        c = np.max(self.power_law(A, mag) / self.power_law(B, mag))
        difference = self.power_law(A, mag) / (self.power_law(B, mag) * c)
        randoms = np.random.uniform(size=len(difference))
        ii = np.zeros((2, 2))
        while len(ii[0]) > 1:
            ii = np.where(randoms > difference)
            mag[ii] = self.prob_dist(m_min, m_max, B[1], len(ii[0]))
            randoms = np.random.uniform(size=len(ii[0]))
            difference = self.power_law(A, mag[ii]) \
                                            / (self.power_law(B, mag[ii]) * c)

        selected_sources = np.sort(mag)[0:size]
        return selected_sources
