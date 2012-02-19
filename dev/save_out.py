# Rory Holmes
# Feb 2012

# This file contains the functions to save out data from the
# self-calibration simulations. All data is saved out as pickles.

# Make Python 3 compatible
from __future__ import division, print_function

# Standard Modules
import numpy as np
import pickle as pickle
import scipy.optimize as opt

# Custom self-cal modules
import self_calibration
import transformations as tran
import god


def fitted_flat_field(p, q, it_no, data_dir, verbose=False):
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

    filename = '%s/FF/%0*d_ff.p' % (data_dir, 3, it_no)
    if verbose:
        print("Saving out fitted flat-field data to {0}...".format(filename))

    x = np.linspace(-p['FoV'][0] / 2, p['FoV'][0] / 2, p['ff_samples'][0])
    y = np.linspace(-p['FoV'][1] / 2, p['FoV'][1] / 2, p['ff_samples'][1])
    X, Y = np.meshgrid(x, y)
    # Have to reshape so that evaluate_flat_field() works.
    reshape_x = np.reshape(X, -1)
    reshape_y = np.reshape(Y, -1)
    reshape_z = self_calibration.evaluate_flat_field(p, reshape_x,
                                                            reshape_y, q)
    fitted_ff = np.reshape(reshape_z, (len(X[:, 0]), len(X[0, :])))
    dic = {'x': X, 'y': Y, 'fitted_ff': fitted_ff, 'iteration_number': it_no}
    pickle.dump(dic, open(filename, "wb"))

    if verbose:
        print("...done!")


def camera(p, sky_catalog, measured_catalog, inside_FoV,
            pointing, orientation, data_dir, verbose=False):
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
    if verbose:
        print("Saving out camera image to {0}...".format(filename))

    x_min = - p['FoV'][0] / 2
    y_min = - p['FoV'][1] / 2
    x_max = p['FoV'][0] / 2
    y_max = p['FoV'][1] / 2
    x = np.array([x_min, x_min, x_max, x_max, x_min])
    y = np.array([y_min, y_max, y_max, y_min, y_min])
    alpha, beta = tran.fp2sky(x, y, pointing, orientation)
    dic = {'measured_catalog.x': measured_catalog.x, 'measured_catalog.y':
            measured_catalog.y, 'sky_catalog.alpha': sky_catalog.alpha,
            'sky_catalog.beta': sky_catalog.beta, 'pointing': pointing,
            'orientation': orientation, 'sky': p['sky_limits'], 'fp_x': x,
            'fp_y': y, 'fp_alpha': alpha, 'fp_beta': beta,
            'inside_FoV': inside_FoV}
    pickle.dump(dic, open(filename, "wb"))

    if verbose:
        print("...done!")


def invvar(obs_cat, data_dir, verbose=False):
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
    if verbose:
        print("Saving out inverse variance to {0}...".format(filename))

    dic = {'counts': obs_cat.counts, 'true_invvar': obs_cat.gods_invvar,
            'reported_invvar': obs_cat.invvar}
    pickle.dump(dic, open(filename, "wb"))

    if verbose:
        print("...done!")


def bestfit_ff(p, data_dir, verbose=False):
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
    if verbose:
        print("Saving out the best fit to God's flat-field to {0}"
                                                        .format(filename))

    order = p['flat_field_order']
    ff_samples = p['ff_samples']
    temp_x = np.linspace(-0.5 * p['FoV'][0], 0.5 * p['FoV'][0], ff_samples[0])
    temp_y = np.linspace(-0.5 * p['FoV'][1], 0.5 * p['FoV'][1], ff_samples[1])
    X, Y = np.meshgrid(temp_x, temp_y)
    x = np.reshape(X, -1)
    y = np.reshape(Y, -1)
    g = self_calibration.evaluate_flat_field_functions(x, y, order)
    god_ff = god.flat_field(p, x, y)
    a = np.zeros((order + 1) * (order + 2) / 2)
    a[0] = 1
    a[3] = -0.2
    a[5] = 0.5
    print("Fitting god's flat-field with basis...")
    fitted_parameters = opt.fmin_bfgs(self_calibration.compare_flats, a,
                        args=(god_ff, g, x, y), gtol=1e-4, maxiter=1e6)
    print("Done: Fitted Parameters: ", fitted_parameters)
    fitted_parameters = self_calibration.normalize_flat_field(p,
                                                        fitted_parameters)
    bestfit_ff = np.reshape(self_calibration.evaluate_flat_field(p, x, y, \
                            fitted_parameters), (len(X[:, 0]), len(X[0, :])))
    dic = {'x': X, 'y': Y, 'bestfit_ff': bestfit_ff,
                'fit_parameters': fitted_parameters}
    pickle.dump(dic, open(filename, "wb"))

    if verbose:
        print("...done!")


def god_ff(p, data_dir, verbose=False):
    ''' Saves out the God's flat-field at sample points across the focal plane

    p                   :   dictionary
        Parameters used in the self-calibration simulation
    data_dir            :   str
        The directory name in which the output data should be stored
    verbose             :   boolean
        True to run function in verbose mode
    '''

    filename = '{0}/god_ff.p'.format(data_dir)
    if verbose:
        print("Saving out God's flat-field to {0}".format(filename))

    ff_samples = p['ff_samples']
    x = np.linspace(-0.5 * p['FoV'][0], 0.5 * p['FoV'][0], ff_samples[0])
    y = np.linspace(-0.5 * p['FoV'][1], 0.5 * p['FoV'][1], ff_samples[1])
    X, Y = np.meshgrid(x, y)
    god_ff = god.flat_field(p, X, Y)
    dic = {'x': X, 'y': Y, 'god_ff': god_ff}
    pickle.dump(dic, open(filename, "wb"))

    if verbose:
        print("...done!")
