# Rory Holmes
# Feb 2012

# This file contains the functions to self-calibrate the photometric catalog
# in the self-calibration simulations.

# Make Python 3 compatible
from __future__ import division, print_function

# Standard Modules
import numpy as np

# Custom self-cal modules
import save_out
import analysis
import god


def self_calibration(params, observation_catalog, sky_catalog,
                                data_dir, plots=False, verbose=False):
    order = params['flat_field_order']
    q = np.array([1])
    stop_condition = params['stop_condition']
    max_iterations = params['max_iterations']
    chi2 = 1e9
    old_chi2 = 1e10
    count = 0
    next_plot_iteration = 1

    if plots:
        save_out.fitted_flat_field(params, q, count, data_dir)

    while ((abs(chi2 - old_chi2) > stop_condition) and \
                                            (count < max_iterations)):
        count += 1
        temp_chi2 = chi2
        s, s_invvar = s_step(params, observation_catalog, q)
        q, q_invvar, chi2 = q_step(params, observation_catalog, s,
                                    order, count, plots=plots, verbose=verbose)
        old_chi2 = temp_chi2
        # Calculate rms error in stars
        indx = [s != 0]
        rms = analysis.rms_error(s[indx], sky_catalog.flux[indx])
        bdness = analysis.badness(params, q)
        bdness_bestfitff = analysis.best_in_basis(params, q, data_dir)
        print(("{0}: RMS = {1:.6f} %, Badness = {2:.6f} %, " 
                "BestInBasis_Badness = {3:.6f} %, chi2 = {4:.2f} ({5})")
                .format(count, rms, bdness, bdness_bestfitff, chi2,
                observation_catalog.size))
        print(q)
        if (plots and (count == next_plot_iteration)) or \
                                (abs(chi2 - old_chi2) < stop_condition):
            save_out.fitted_flat_field(params, q, count, data_dir)
            next_plot_iteration *= 2
    return np.array([count, rms, bdness, bdness_bestfitff, chi2])


def evaluate_flat_field_functions(x, y, order, verbose=False):
    L = (order + 1) * (order + 2) / 2
    g = np.zeros((len(x), L))
    l = 0
    for n in range(order + 1):
        for m in range(n + 1):
            g[:, l] = (x ** (n - m)) * (y ** m)
            l += 1
    return g


def evaluate_flat_field(params, x, y, q, verbose=False):
    # Calculate required order
    order = int(np.around(np.sqrt(0.25 + 2 * len(q)) - 1.5))
    assert(len(q) == ((order + 1) * (order + 2) / 2))
    g = evaluate_flat_field_functions(x, y, order)
    return np.dot(g, q)


def normalize_flat_field(params, q, verbose=False):
    fit_mean = average_over_ff(params, evaluate_flat_field, (q))
    god_mean = average_over_ff(params, god.flat_field,
                                            (god.flat_field_parameters()))
    return (q * god_mean / fit_mean)


def average_over_ff(p, func, args, verbose=False):
    nalpha, nbeta = p['ff_samples']
    dalpha = p['FoV'][0] / (nalpha - 1)
    dbeta = p['FoV'][1] / (nbeta - 1)
    x = np.arange(-p['FoV'][0] / 2 + dalpha / 2, p['FoV'][0] / 2, dalpha)
    y = np.arange(-p['FoV'][1] / 2 + dbeta / 2, p['FoV'][1] / 2, dbeta)
    X, Y = np.meshgrid(x, y)
    temp_x = np.reshape(X, -1)
    temp_y = np.reshape(Y, -1)
    temp_ff = func(p, temp_x, temp_y, args)
    ff = np.reshape(temp_ff, (len(X[:, 0]), len(X[0, :])))
    return np.mean(ff)


def s_step(params, obs_cat, q, verbose=False):
    ff = evaluate_flat_field(params, obs_cat.x, obs_cat.y, q)
    fcss = ff * obs_cat.counts * obs_cat.invvar
    ffss = ff * ff * obs_cat.invvar
    max_star = np.max(obs_cat.k)
    s = np.zeros(max_star + 1)
    s_invvar = np.zeros(max_star)
    for ID in range(max_star):
        indx = (obs_cat.k == ID)
        denominator = np.sum(ffss[indx])
        s_invvar[ID] = denominator
        if denominator > 0.:
            s[ID] = np.sum(fcss[indx]) / denominator
    return (s, s_invvar)


def q_step(params, obs_cat, s, order, iteration_number,
                                            plots=False, verbose=False):
    g = evaluate_flat_field_functions(obs_cat.x, obs_cat.y, order)
    ss = s[obs_cat.k]
    q_invvar = np.dot(np.transpose(g) * ss * ss * obs_cat.invvar, g)
    numerator = np.sum(np.transpose(g) * ss * obs_cat.counts * obs_cat.invvar,
                                                        axis=1)
    q = np.dot(np.linalg.inv(q_invvar), numerator)
    q = normalize_flat_field(params, q)
    np.set_printoptions(precision=2)
    chi2 = np.sum((np.dot(g, q) * ss - obs_cat.counts) ** 2 * obs_cat.invvar)
    return (q, q_invvar, chi2)


def compare_flats(a, Z_true, g, X, Y, verbose=False):
    error = np.sum((Z_true - np.dot(g, a)) ** 2)
    return error
