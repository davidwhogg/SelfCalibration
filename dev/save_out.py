# Rory Holmes
# Feb 2012

# This file contains the functions to save out data from the
# self-calibration simulations. All data is saved out as pickles.

# Make Python 3 compatible
from __future__ import division, print_function

# Standard Modules
import pickle
import numpy as np

# Custom Modules
import self_calibration
import transformations


def parameters(data_dir, params, verbose=False):
    ''' Saves out the source catalog object

    Input
    -----
    data_dir            :   string
        The output directory path for the simulation run
    params              :   dictionary
        The parameter dictionary
    '''

    filename = '{0}/parameters.p'.format(data_dir)
    if verbose:
        print("Saving out parameters to {0}".format(filename))
    pickle.dump(params, open(filename, "wb"))
    if verbose:
        print("...done!")


def source_catalog(data_dir, source_catalog, verbose=False):
    ''' Saves out the source catalog object

    Input
    -----
    data_dir            :   string
        The output directory path for the simulation run
    source_catalog      :   obj
        The catalog object # *.k, *.alpha, *.beta, *.mag, *.size, *.flux,
    '''

    filename = '{0}/source_catalog.p'.format(data_dir)
    if verbose:
        print("Saving out source catalog to {0}".format(filename))
    pickle.dump(source_catalog, open(filename, "wb"))
    if verbose:
        print("...done!")


def survey_strategy(data_dir, survey, verbose=False):
    ''' Saves out the survey strategy

    Input
    -----
    data_dir            :   string
        The output directory path for the simulation run
    survey_catalog      :   numpy array
        The survey strategy, each row represents a new pointing in format:
        ID, RA, Dec, Orientation
    '''
    filename = '{0}/survey.txt'.format(data_dir)
    if verbose:
        print("Saving out source catalog to {0}".format(filename))
    np.savetxt(filename, survey)
    if verbose:
        print("...done!")


def camera(data_dir, sky_catalog, measured_catalog, FoV,
            pointing, orientation, verbose=False):
    ''' Saves out a single camera exposure, along with the total sky catalog
    and the meas

    Parameters
    ----------
    data_dir            :   string
        The output directory
    sky_catalog         :   object
        Object of the sky catalog (*.ID, *.mag, *.alpha, *.beta, *.size)
    measured_catalog    :   object
        Object of measured sources on focal plane in a single pointing
        (*.size, *.ID, *.x, *.y, *true_invvar, *.counts, *.invvar)
    FoV                 :   float array
        The imagers field-of-view in degrees (dalpha, dbeta)
    pointing            :   np.array
        Array of telescope pointing (alpha, beta)
    orientation         :   float
        Rotation of the telescope
    verbose             :   boolean
        True to run function in verbose mode
    '''

    filename = '{0}/camera_image.p'.format(data_dir)

    if verbose:
        print("Saving out camera image to {0}...".format(filename))

    x = 0.5 * np.array([-FoV[0], -FoV[0], FoV[0], FoV[0], -FoV[0]])
    y = 0.5 * np.array([-FoV[1], FoV[1], FoV[1], -FoV[1], -FoV[1]])
    alpha, beta = transformations.fp2sky(x, y, pointing, orientation)

    dic = {'sources_x': measured_catalog.x,
            'sources_y': measured_catalog.y,
            'sky_catalog': sky_catalog,
            'pointing': pointing,
            'orientation': orientation,
            'fp_x': x,
            'fp_y': y,
            'fp_alpha': alpha,
            'fp_beta': beta}
    pickle.dump(dic, open(filename, "wb"))

    if verbose:
        print("...done!")


def survey(data_dir, survey_catalog, verbose=False):
    ''' Saves out the survey catalog

    Input
    -----
    data_dir            :   string
        The output directory path for the simulation run
    survey_catalog      :   obj
        The survey object # *.k, *.alpha, *.beta, *.mag, , *.x, *.y *.size,
        *.flux,
    '''

    filename = '{0}/survey_catalog.p'.format(data_dir)
    if verbose:
        print("Saving out source catalog to {0}".format(filename))
    pickle.dump(survey_catalog, open(filename, "wb"))
    if verbose:
        print("...done!")


def fitted_flat_field(q, FoV, ff_samples, it_no, data_dir, verbose=False):
    ''' Saves out data for the fitted flat-field at sample points across the
    focal plane.

    Parameters
    ----------
    q       :   numpy array
        The fitted flat-field parameters
    FoV     :   float array
        The imagers field of view in degrees (dalpha, dbeta)
    it_no   :   int
        The self-calibration iteration number
    data_dir:   string
        The directory name in which the output data should be stored
    verbose :   boolean
        True to run function in verbose mode
    '''

    filename = '{0}/FF/{1}_ff.p'.format(data_dir, it_no)
    if verbose:
        print("Saving out fitted flat-field data to {0}...".format(filename))

    x = np.linspace(-FoV[0] / 2, FoV[0] / 2, ff_samples[0])
    y = np.linspace(-FoV[1] / 2, FoV[1] / 2, ff_samples[1])
    X, Y = np.meshgrid(x, y)
    z = self_calibration.evaluate_flat_field(X.flatten(), Y.flatten(), q)
    fitted_ff = np.reshape(z, (len(X[:, 0]), len(X[0, :])))

    dic = {'x': X, 'y': Y, 'fitted_ff': fitted_ff, 'iteration_number': it_no}
    pickle.dump(dic, open(filename, "wb"))

    if verbose:
        print("...done!")
