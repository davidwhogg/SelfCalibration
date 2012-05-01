# Rory Holmes (MPIA) / David Hogg (NYU)
# 2011 - 2012

# This file contains the functions to self-calibrate the photometric catalog
# in the self-calibration simulations. These functions both fit for the source
# fluxes and also the instrument response.

# Make Python 3 compatible
from __future__ import division, print_function

# Standard Modules
import numpy as np

# Custom Modules
import true_functions
import save_out
import analysis


def self_calibration(obs_cat, sky_catalog, order, FoV,
                        ff_samples, stop_condition, max_iterations,
                        best_fit_parameters, data_dir=False, verbose=False):
    ''' This functions performs the self-calibration procedure on the survey.

    Parameters
    ----------
    obs_cat                 :       object
        The observation catalog
    sky_catalog             :       object
        The *true* sky catalog
    order                   :       Integer
        The order of the polynomial flat-field used to fit the instrument
        response in the self-calibration procedure
    FoV                     :       float array
        The imagers field-of-view in degrees (dalpha, dbeta)
    ff_samples              :       float array
        The number of points to sample the focal plane with
        (alpha-axis, beta-axis)
    stop_condition          :       float
        The stop condition for the self-calibration procedure and the
        best-in-basis fitting (stop when difference is less than 2 times this)
    max_iterations          :       float
        The maximum number of iterations in the self-calibration procedure and
        the best-in-basis fitting
    best_fit_parameters     :       numpy array
        The best-in-basis fit to the true flat-field
    data_dir                :       boolean / string
        If set to False, the self-calibration simulations do not save out
        any data (only returning the results to the calling program). If a
        string is give, all data is saved out to this directory.
    verbose                 :       boolean
        Set to True to run the simulations in verbose mode
    Returns
    -------
    out                     :       numpy array
        The results of the self-calibration procedure [the number of iteration
        steps required, the rms error of the final fitted source fluxes, their
        badness of the fitted flat-field, the best-in-basis badness of the
        fitted flat-field, the chi2 of the fit].
    '''
    q = np.array([1])
    chi2 = 1e9
    old_chi2 = 1e10
    count = 0
    next_plot_iteration = 1

    if data_dir:
        save_out.fitted_flat_field(q, FoV, ff_samples, count, data_dir,
                                                                    verbose)

    while ((abs(chi2 - old_chi2) > stop_condition) and \
                                            (count < max_iterations)):
        count += 1
        temp_chi2 = chi2
        s, s_invvar = s_step(obs_cat, q)
        q, q_invvar, chi2 = q_step(obs_cat, s, order, FoV, ff_samples)
        old_chi2 = temp_chi2

        indx = [s != 0]
        rms = analysis.rms_error(s[indx], sky_catalog.flux[indx], verbose)
        bdness = analysis.badness(q, FoV, ff_samples, verbose)
        bdness_bestfitff = analysis.best_in_basis(q, FoV, ff_samples,
                                        best_fit_parameters, verbose)

        if verbose:
            print('Fitted parameters: {0}'.format(q))
            print(("{0}: RMS = {1:.6f} %, Badness = {2:.6f} %, "
                    "BestInBasis_Badness = {3:.6f} %, chi2 = {4:.2f} ({5})")
                    .format(count, rms, bdness, bdness_bestfitff, chi2,
                    obs_cat.size))

        if (data_dir and (count == next_plot_iteration)) or \
                                (abs(chi2 - old_chi2) < stop_condition):
            save_out.fitted_flat_field(q, FoV, ff_samples, count,
                                                            data_dir, verbose)
            next_plot_iteration *= 2
    
    if count = max_iterations:
        count = 0  # did not converge error flag

    return np.array([count, rms, bdness, bdness_bestfitff, chi2])


def s_step(obs_cat, q):
    ''' This functions returns a weighted mean of the sources observed multiple
    times in the survey.

    Parameters
    ----------
    obs_cat               :   object
        The observation catalog
    q               :   numpy.array
        The instrument response parameters

    Returns
    -------
    out             :   tuple
        The weighted means and their uncertainties
    '''
    flat_field = evaluate_flat_field(obs_cat.x, obs_cat.y, q)
    fcss = flat_field * obs_cat.counts * obs_cat.invvar
    ffss = flat_field * flat_field * obs_cat.invvar
    s = np.zeros(obs_cat.k.max() + 1)
    s_invvar = np.zeros(obs_cat.k.max())
    for ID in range(obs_cat.k.max()):
        indx = (obs_cat.k == ID)
        denominator = np.sum(ffss[indx])
        s_invvar[ID] = denominator
        if denominator > 0.:
            s[ID] = np.sum(fcss[indx]) / denominator
    return (s, s_invvar)


def q_step(obs_cat, s, order, FoV, ff_samples):
    ''' This functions fits for the flat-field parameters based on the given
    model source fluxes.

    Parameters
    ----------
    obs_cat               :   object
        The observation catalog
    s               :   numpy.array
        The source fluxes
    order           :   Integer
        The order of the polynomial flat-field used to fit the instrument
        response in the self-calibration procedure
    FoV    :    float array
        The imagers field-of-view in degrees (dalpha, dbeta)
    ff_samples  :   float array
        The number of points to sample the focal plane with
        (alpha-axis, beta-axis)

    Returns
    -------
    out             :   tuple
        The results of the flat-field fitting procedure (numpy array of the
        fitted flat-field paramters, their uncertainties and the chi2 of the
        fit)
    '''
    g = evaluate_flat_field_functions(obs_cat.x, obs_cat.y, order)
    ss = s[obs_cat.k]
    q_invvar = np.dot(np.transpose(g) * ss * ss * obs_cat.invvar, g)
    numerator = np.sum(np.transpose(g) * ss * obs_cat.counts * obs_cat.invvar,
                                                        axis=1)
    q = np.dot(np.linalg.inv(q_invvar), numerator)
    q = normalize_flat_field(q, FoV, ff_samples)
    chi2 = np.sum((np.dot(g, q) * ss - obs_cat.counts) ** 2 * obs_cat.invvar)
    return (q, q_invvar, chi2)


def evaluate_flat_field_functions(x, y, order):
    ''' This functions returns the g values for the flat-field calculation.

    Parameters
    ----------
    x               :   numpy.array
        The x-coordinates of the sample points
    y               :   numpy.array
        The y-coordinates of the sample points
    order           :   Integer
        The order of the polynomial flat-field used to fit the instrument
        response in the self-calibration procedure

    Returns
    -------
    out             :   numpy array
        the g parameter
    '''

    L = (order + 1) * (order + 2) / 2
    g = np.zeros((len(x), L))
    l = 0
    for n in range(order + 1):
        for m in range(n + 1):
            g[:, l] = (x ** (n - m)) * (y ** m)
            l += 1
    return g


def evaluate_flat_field(x, y, q):
    ''' This functions returns the flat-field values

    Parameters
    ----------
    x               :   numpy.array
        The x-coordinates of the sample points
    y               :   numpy.array
        The y-coordinates of the sample points
    q               :   numpy.array
        The flat-field parameters

    Returns
    -------
    out             :   numpy array
        the flat-field values at the given sample points
    '''

    assert x.shape == y.shape
    order = int(np.around(np.sqrt(0.25 + 2 * len(q)) - 1.5))
    assert(len(q) == ((order + 1) * (order + 2) / 2))
    g = evaluate_flat_field_functions(x, y, order)
    return np.dot(g, q)


def normalize_flat_field(q, FoV, ff_samples):
    ''' This functions modifies the flat-field parameters (q) in such a way
    to normalize the fitted flat-field compared to the *true* flat-field. This
    is necessary as the self-calibration procedure only constrains the
    relative calibration, not the absolute (it is not possible to know if
    all the sources are 50% fainter, or if the instrument response is 50%
    lower). In reality this degeneracy will be solved with absolute standards.

    Parameters
    ----------
    q               :   numpy.array
        The flat-field parameters

    Returns
    -------
    out             :   numpy array
        the modified flat-field parameters, which describe a normalized
        flat-field
    '''
    fit_mean = average_over_ff(FoV, ff_samples, evaluate_flat_field, (q))
    true_mean = average_over_ff(FoV, ff_samples, true_functions.flat_field,
                                    (true_functions.flat_field_parameters()))
    return (q * true_mean / fit_mean)


def average_over_ff(FoV, ff_samples, func, ff_params):
    ''' This functions returns the average of the flat-field values at the
    sample points given.

    Parameters
    ----------
    FoV             :   float array
        The simulate imager's field-of-view in degrees [alpha, beta]
    ff_samples      :   integer array
        The number of sample points on the focal plane for the best-in-basis
        fitting and all of the badness measures
    func            :   function
        The function describing the instrument flat-field
    ff_params       :   int array
        The flat-field parameters describing flat-field

    Returns
    -------
    out     :   numpy array
        the flat-field values at the given sample points
    '''

    x = np.linspace(-0.5 * FoV[0], 0.5 * FoV[0], ff_samples[0])
    y = np.linspace(-0.5 * FoV[1], 0.5 * FoV[1], ff_samples[1])
    X, Y = np.meshgrid(x, y)
    return np.mean(func(X.flatten(), Y.flatten(), ff_params))
