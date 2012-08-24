# Rory Holmes (MPIA) / David Hogg (NYU)
# 2011 - 2012

# This file contains all the functions and classes used to perform the imaging
# survey in the self-calibration simulations.

# Standard Modules
import numpy as np
import os
import sys

# Custom Modules
import transformations
import true_functions
import save_out


def survey(survey_file, sky_catalog, FoV, eta, delta, epsilon_max,
                                                data_dir=False, verbose=False):
    ''' This functions performs the survey on the synthetic sky and
    returns the observation catalog.

    Input
    -----
    sky_file                :   string
        The path to the survey file, which must be in the format (angles
        in degrees): observation_number alpha beta orientation
    sky_catalog             :   object
        The sky catalog object to survey_file
    FoV                     :   float array
        The images field-of-view (dalpha, dbeta)
    eta, delta, epsilon_max :   all floats
        The parameters describing the noise model according to:
    sigma ** 2 = (1 + epsilon) * delta ** 2 + eta ** 2 * count_rate ** 2
    data_dir                :   boolean / string
        If set to False, the self-calibration simulations do not save out
        any data (only returning the results to the calling program). If a
        string is give, all data is saved out to this directory.
    verbose                 :   boolean
        Set to True to run all functions in verbose mode

    Return
    ------
    out:                    :   object
        The observation catalog
    '''
    if verbose:
        print("Loading survey from {0}...".format(survey_file))
    if os.path.exists(survey_file):
        pointings = np.loadtxt(survey_file)
        number_pointings = len(pointings[:, 0])
    else:
        print('Error - no survey file!')
        sys.exit()
    if verbose:
        print("...{0} pointings loaded!".format(number_pointings))
    if data_dir:
        save_out.survey_strategy(data_dir, pointings, verbose=verbose)

    if verbose:
        print("Surveying sky...")
    obs_cat = None
    for indx in range(pointings.shape[0]):
        si = single_image(sky_catalog, pointings[indx, 1:3],
                            pointings[indx, 3], FoV, eta,
                            delta, epsilon_max, data_dir,
                            verbose)
        if obs_cat is None:
            obs_cat = si
        else:
            obs_cat.append(si)

    if verbose:
        print("{0} individual source measurements!".format(obs_cat.size))
    return obs_cat


class CameraCatalog:

    def __init__(self, sky_catalog, pointing, orientation):
        ''' This class generates a Camera Catalog object, which is just the
        Survey Catalog transfered into focal plane coordinates with the given
        telescope pointing and orientation information.

        Input
        -----
        sky_catalog             :       object
            The sky catalog
        pointing                :       numpy array [alpha, beta]
            The telescope pointing coordinates in degrees
        orientation             :       float
            The telescope orientation in degrees
        '''

        self.k = sky_catalog.k
        self.mag = sky_catalog.mag
        self.alpha = sky_catalog.alpha
        self.beta = sky_catalog.beta
        self.is_in_analysis_region = sky_catalog.is_in_analysis_region
        self.x, self.y = \
            transformations.sky2fp(sky_catalog.alpha, sky_catalog.beta,\
                                                        pointing, orientation)
        self.size = sky_catalog.size
        self.flux = transformations.mag2flux(self.mag)


class MeasuredCatalog:
    def __init__(self, camera_catalog, FoV, inside_FoV, eta, delta,
                                                                epsilon_max):
        ''' This class generates the measured catalog object. Sources within
        in the imagers field-of-view are identified and then "measurements"
        are performed using the *true* instrument response model
        Input
        -----
        camera_catalog          :       object
            The camera catalog
        FoV                     :       float array
            The imagers field-of-view (dalpha, dbeta)
        inside_FoV              :       boolean array
            Array identifying sources within the instruments field of view
        eta, delta, epsilon_max :   all floats
        The parameters describing the noise model according to:
        sigma ** 2 = (1 + epsilon) * delta ** 2 + eta ** 2 * count_rate ** 2
        '''

        self.size = len(inside_FoV[0])
        self.k = camera_catalog.k[inside_FoV]
        self.alpha = camera_catalog.alpha[inside_FoV]
        self.beta = camera_catalog.beta[inside_FoV]
        self.is_in_analysis_region = camera_catalog.is_in_analysis_region[inside_FoV]
        self.x = camera_catalog.x[inside_FoV]
        self.y = camera_catalog.y[inside_FoV]
        flat = true_functions.flat_field(self.x, self.y, FoV)
        true_counts = camera_catalog.flux[inside_FoV] * flat
        self.true_invvar = self.true_invvar_func(true_counts,
                                                    eta, delta, epsilon_max)
        self.invvar = self.reported_invvar(true_counts, eta, delta)
        self.counts = true_counts + np.random.normal(size=self.size) \
                        / np.sqrt(self.true_invvar)

    def reported_invvar(self, true_counts, eta, delta):
        ''' This function calculates the *assumed* inverse variance for each
        measurement that is assumed during the self-calibration procedure
        '''
        var = (delta ** 2) + (eta ** 2) * (true_counts ** 2)
        return 1. / var

    def true_invvar_func(self, true_counts, eta, delta, epsilon_max):
        ''' This function calculates the *true*  inverse invariance for each
        measurement
        '''
        epsilon = np.random.uniform(low=0., high=epsilon_max, size=self.size)
        var = (delta ** 2) * (1 + epsilon) + (eta ** 2) * (true_counts ** 2)
        return 1. / var

    def append(self, other):
        self.size = self.size + other.size
        self.k = np.append(self.k, other.k)
        self.alpha = np.append(self.alpha, other.alpha)
        self.beta = np.append(self.beta, other.beta)
        self.is_in_analysis_region = np.append(self.is_in_analysis_region,
                                                other.is_in_analysis_region)
        self.x = np.append(self.x, other.x)
        self.y = np.append(self.y, other.y)
        self.counts = np.append(self.counts, other.counts)
        self.invvar = np.append(self.invvar, other.invvar)
        self.true_invvar = np.append(self.true_invvar, other.true_invvar)


def single_image(sky_catalog, pointing, orientation, FoV, eta, delta,
                                epsilon_max, data_dir=False, verbose=False):
    ''' This functions performs a single image of the sky at the given pointing
    coordinates.

    Input
    -----
    sky_catalog             :   object
        The sky catalog object to survey_file
    pointing                :   numpy array
        The coordinates of the telecsope pointing in degrees (alpha, beta)
    orientation:            :   float
        The orientation of the telescope in degrees
    FoV                     :   numpy array
        The images field-of-view (dalpha, dbeta)
    eta, delta, epsilon_max :   all floats
        The parameters describing the noise model according to:
    sigma ** 2 = (1 + epsilon) * delta ** 2 + eta ** 2 * count_rate ** 2
    data_dir                :   boolean / string
        If set to False, the self-calibration simulations do not save out
        any data (only returning the results to the calling program). If a
        string is give, all data is saved out to this directory.

    Return
    ------
    out:                    :   object
        The single image object
    '''
    camera_catalog = CameraCatalog(sky_catalog, pointing, orientation)
    inside_FoV = np.where((-0.5 * FoV[0] < camera_catalog.x) \
                    & (camera_catalog.x < 0.5 * FoV[0]) \
                    & (-0.5 * FoV[1] < camera_catalog.y) \
                    & (camera_catalog.y < 0.5 * FoV[1]))
    measured_catalog = MeasuredCatalog(camera_catalog, FoV, inside_FoV, eta,
                                                            delta, epsilon_max)
    if data_dir:
        if os.path.exists((data_dir + '/camera_image.p')) is False:
            if (orientation > 30) and (pointing[0] > -1) and (pointing[0] < 1)\
                            and (pointing[1] > -1) and (pointing[1] < 1):
                save_out.camera(data_dir, sky_catalog, measured_catalog,
                                   FoV, pointing, orientation, verbose=verbose)

    return measured_catalog
