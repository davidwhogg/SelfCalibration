# Rory Holmes
# Feb 2012

# This file contains the functions to save out data from the
# self-calibration simulations. All data is saved out as pickles.

# Make Python 3 compatible
from __future__ import division, print_function

# Standard Modules
import numpy as np
import pickle as pickle

# Custom self-cal modules
import self_calibration
import transformations as tran
import god


def parameters(p):
    ''' Saves out the parameters used in the self-calibration simulation.

    Parameters
    ----------
    p          :   dictionary
        the dictionary containing all the simulation parameters
    '''

    filename = "{0}/parameters.p".format(p["data_dir"])
    if p["verbose"]:
        print("Saving out fitted flat-field data to {0}...".format(filename))

    pickle.dump(p, open(filename, "wb"))

    if p["verbose"]:
        print("...done!")


def fitted_flat_field(p, q, it_no, data_dir):
    ''' Saves out data for the fitted flat-field at sample points across the
    focal plane.

    Parameters
    ----------
    p       :   dictionary
        Parameters used in the self-calibration simulation
    q       :   np.array
        The fitted flat-field parameters
    it_no   :   int
        The self-calibration iteration number
    data_dir:   str
        The directory name in which the output data should be stored
    verbose :   boolean
        True to run function in verbose mode
    '''

    filename = '{0}/FF/{1}_ff.p'.format(data_dir, it_no)
    if p["verbose"]:
        print("Saving out fitted flat-field data to {0}...".format(filename))

    x = np.linspace(-p['FoV'][0] / 2, p['FoV'][0] / 2, p['ff_samples'][0])
    y = np.linspace(-p['FoV'][1] / 2, p['FoV'][1] / 2, p['ff_samples'][1])
    X, Y = np.meshgrid(x, y)
    z = self_calibration.evaluate_flat_field(p, X.flatten(), \
                        y.flatten(), q, verbose=False)
    fitted_ff = np.reshape(z, (len(X[:, 0]), len(X[0, :])))

    dic = {'x': X, 'y': Y, 'fitted_ff': fitted_ff, 'iteration_number': it_no}
    pickle.dump(dic, open(filename, "wb"))

    if p["verbose"]:
        print("...done!")


def camera(p, sky_catalog, measured_catalog, inside_FoV,
            pointing, orientation, data_dir,):
    ''' Saves out a single camera exposure, along with the total sky catalog
    and the meas

    Parameters
    ----------
    p                   :   dictionary
        Parameters used in the self-calibration simulation
    sky_catalog         :   object
        Object of the sky catalog (*.ID, *.mag, *.alpha, *.beta, *.size)
    measured_catalog    :   object
        Object of measured sources on focal plane in a single pointing
        (*.size, *.ID, *.x, *.y, *gods_invvar, *.counts, *.invvar)
    inside_FoV          :   np.array
        Boolean array stating if source is inside the camera field-of-view
    pointing            :   np.array
        Array of telescope pointing (alpha, beta)
    orientation         :   float
        Rotation of the telescope
    data_dir            :   str
        The directory name in which the output data should be stored
    verbose             :   boolean
        True to run function in verbose mode
    '''

    filename = data_dir + '/camera_image.p'
    if p["verbose"]:
        print("Saving out camera image to {0}...".format(filename))

    x = 0.5 * np.array([p['FoV'][0], -p['FoV'][0], p['FoV'][0],
                            p['FoV'][0], -p['FoV'][0]])
    y = 0.5 * np.array([-p['FoV'][1], p['FoV'][1], p['FoV'][1],
                            -p['FoV'][1], -p['FoV'[1]]])
    alpha, beta = tran.fp2sky(x, y, pointing, orientation)

    dic = {'measured_catalog': measured_catalog,
            'sky_catalog': sky_catalog,
            'pointing': pointing,
            'orientation': orientation,
            'sky': p['sky_limits'],
            'fp_x': x,
            'fp_y': y,
            'fp_alpha': alpha,
            'fp_beta': beta,
            'inside_FoV': inside_FoV}
    pickle.dump(dic, open(filename, "wb"))

    if p["verbose"]:
        print("...done!")


def invvar(p, obs_cat, data_dir):
    ''' Saves out the true inverse variance, the reported inverse variance and
    the counts from each observation.

    Parameters
    ----------
    obs_cat         :   object
        Object of the sky catalog (*.size, *.pointing_ID, *.star_ID, *.flux,
        *.invvar, *.x, *.y)
    data_dir        :   str
        The directory name in which the output data should be stored
    verbose         :   boolean
        True to run function in verbose mode
    '''

    filename = data_dir + '/invvar.p'
    if p["verbose"]:
        print("Saving out inverse variance to {0}...".format(filename))

    dic = {'counts': obs_cat.counts, 'true_invvar': obs_cat.gods_invvar,
            'reported_invvar': obs_cat.invvar}
    pickle.dump(dic, open(filename, "wb"))

    if p["verbose"]:
        print("...done!")


def best_fit_ff(p, data_dir):
    ''' Calculates the best fit possible to God's flat-field with the basis
    used to model it during the self-calibration procedure. Saves out the
    data at sample points across the focal plane

    p                   :   dictionary
        Parameters used in the self-calibration simulation
    data_dir            :   str
        The directory name in which the output data should be stored
    verbose             :   boolean
        True to run function in verbose mode
    '''

    filename = data_dir + '/bestfit_ff.p'
    if p["verbose"]:
        print("Saving out the best fit to God's flat-field to {0}"
                                                        .format(filename))

    ff_samples = p['ff_samples']
    x = np.linspace(-0.5 * p['FoV'][0], 0.5 * p['FoV'][0], ff_samples[0])
    y = np.linspace(-0.5 * p['FoV'][1], 0.5 * p['FoV'][1], ff_samples[1])
    X, Y = np.meshgrid(x, y)
    z = self_calibration.evaluate_flat_field(p,
                                x.flatten(), y.flatten(),
                                p['best_fit_params'])
    best_fit_ff = np.reshape(z, (len(X[:, 0]), len(X[0, :])))

    dic = {'x': X, 'y': Y, 'best_fit_ff': best_fit_ff}
    pickle.dump(dic, open(filename, "wb"))

    if p["verbose"]:
        print("...done!")


def god_ff(p, data_dir):
    ''' Saves out the God's flat-field at sample points across the focal plane

    p                   :   dictionary
        Parameters used in the self-calibration simulation
    data_dir            :   str
        The directory name in which the output data should be stored
    verbose             :   boolean
        True to run function in verbose mode
    '''

    filename = '{0}/god_ff.p'.format(data_dir)
    if p["verbose"]:
        print("Saving out God's flat-field to {0}".format(filename))

    ff_samples = p['ff_samples']
    x = np.linspace(-0.5 * p['FoV'][0], 0.5 * p['FoV'][0], ff_samples[0])
    y = np.linspace(-0.5 * p['FoV'][1], 0.5 * p['FoV'][1], ff_samples[1])
    X, Y = np.meshgrid(x, y)
    god_ff = god.flat_field(p, X, Y)
    dic = {'x': X, 'y': Y, 'god_ff': god_ff}
    pickle.dump(dic, open(filename, "wb"))

    if p["verbose"]:
        print("...done!")


def source_catalog(p, source_catalog):
    ''' Saves out the source catalog object

    p                   :   dictionary
        Parameters used in the self-calibration simulation
    source_catalog      :   obj
        The catalog object # *.k, *.alpha, *.beta, *.mag, *.size, *.flux,
        *.epsilon
    '''

    filename = '{0}/source_catalog.p'.format(p["data_dir"])
    if p["verbose"]:
        print("Saving out source catalog to {0}".format(filename))
    pickle.dump(source_catalog, open(filename, "wb"))
    if p["verbose"]:
        print("...done!")
