#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Rory Holmes
# Feb 2012

# This file contains the functions to create the figures from the
# self-calibration simulations.

# Make Python 3 compatible
from __future__ import division, print_function

import matplotlib.pyplot as plt
import glob
import pickle
import numpy as np
import os
import sys
import string

import true_functions
import self_calibration
import transformations

GOLDEN_RATIO = 1.61803399


def get_dir_path():
    ''' The function returns the path to the data directory given on the
    command line as this script was called

    Input
    -----

    Return
    ------
    out         :       string
        The data directory
    '''

    if len(sys.argv) == 2:
        dir_path = sys.argv[1]
        if dir_path[-1] == '/':
            dir_path = dir_path[:-1]
        print("Running plotting routine on {0} ...".format(dir_path))
    else:
        print("Error - no plotting directory given!")
        sys.exit()
    return dir_path


def source_catalog(filename, sky_limits, density, m_min, m_max, A,
                            fig_width=8.3, suffix='.png', verbose=False):
    ''' This function plots a sky catalog, both spatially and as a function of
    magnitude

    Input
    -----
    filename            :   string
        The path to the survey files
    sky_limits          :   numpy array
        The area of sky to generate sources in
        [alpha_min, alpha_max, beta_min, beta_max]
    density             :   int
        The maximum number of sources (all magnitude) per unit area
        to generate for the self-calibration simulations
    m_min               :   float
        The saturation limit of the simulated imager
    m_max               :   float
        The 10-sigma detection limit of the simulated imager
    A                   :   numpy array
        The parameters describing the magnitude distribution of the sources
        in the sky, according to: log10(dN/dm) = A + B * mag + C * mag ** 2
    figure_width        :   float
        The width of the figure in inches, default it 8.3
    suffix              :   string
        The format of the saved figure, either '.pdf', '.eps' or '.png'.
        Default is '.png'
    verbose             :   Boolean
        Set to true to run function in verbose mode
    '''

    if verbose:
        print("Plotting Source Catalog from {0}...".format(filename))

    s = pickle.load(open(filename))

    fig = plt.figure(figsize=(figure_width, 0.5 * fig_width))

    ax1 = fig.add_axes([0.075, 0.15, 0.4, 0.8])
    alpha = 10 / density
    if alpha < 0.05:
        alpha = 0.05
    ax1.plot(s.alpha, s.beta, 'k.', alpha=alpha, markersize=1)
    ax1.set_xlabel(r'$\alpha$ (deg$^2$)')
    ax1.set_ylabel(r'$\beta$ (deg$^2$)')
    ax1.set_xlim(1.1 * sky_limits[0], 1.1 * sky_limits[1])
    ax1.set_ylim(1.1 * sky_limits[2], 1.1 * sky_limits[3])

    ax2 = fig.add_axes([0.575, 0.15, 0.4, 0.8])
    mag = np.linspace(m_min, m_max, 20)
    N = 10 ** (A[0] + A[1] * mag + A[2] * mag ** 2)
    ax2.semilogy(mag, 0.5 * N, 'k:')
    area = (sky_limits[1] - sky_limits[0]) * (sky_limits[1] - sky_limits[0])
    mag = np.arange(m_min, m_max + 0.5, 0.5)
    ax2.bar(mag[:-1], np.histogram(s.mag, bins=mag)[0] / area,
                                            width=0.5, color='k', alpha=0.3)
    ax2.set_xlabel(r'Source Magnitude (mag$_{\mathrm{AB}}$)')
    ax2.set_ylabel(r'Density of Sources')

    filename = string.replace(filename, '.p', plot_suffix)
    plt.savefig(filename)
    plt.clf()
    if verbose:
        print("...done!")


def camera_image(filename, sky_limits, density, fig_width=8.3,
                                                suffix='.png', verbose=False):
    ''' This function plots the footprint of the imagers field-of-view on the
    sky (left) and a camera image with the *true* instrument response (right)

    Input
    -----
    filename            :   string
        The path to the survey files
    sky_limits          :   numpy array
        The area of sky to generate sources in
        [alpha_min, alpha_max, beta_min, beta_max]
    density             :   int
        The maximum number of sources (all magnitude) per unit area
        to generate for the self-calibration simulations
    figure_width        :   float
        The width of the figure in inches, default it 8.3
    suffix              :   string
        The format of the saved figure, either '.pdf', '.eps' or '.png'.
        Default is '.png'
    verbose             :   Boolean
        Set to true to run function in verbose mode
    '''

    if verbose:
        print("Plotting Camera Image from {0}...".format(filename))

    dic = pickle.load(open(filename))

    fig = plt.figure(figsize=(fig_width, 0.5 * fig_width))
    fig.clf()

    x = dic['sources_x']
    y = dic['sources_y']
    alpha = dic['sky_catalog'].alpha
    beta = dic['sky_catalog'].beta
    pointing = dic['pointing']
    orientation = dic['orientation']
    fp_x = dic['fp_x']
    fp_y = dic['fp_y']
    fp_alpha = dic['fp_alpha']
    fp_beta = dic['fp_beta']

    ax1 = fig.add_axes([0.075, 0.15, 0.4, 0.8])
    ax1.plot(alpha, beta, 'k.', alpha=(10 / density))
    ax1.plot(fp_alpha, fp_beta, 'k', linewidth=2)
    ax1.set_xlabel(r'Sky Position $\alpha$ (deg)')
    ax1.set_ylabel(r'Sky Position $\beta$ (deg)')
    ax1.set_xlim(1.1 * sky_limits[0], 1.1 * sky_limits[1])
    ax1.set_ylim(1.1 * sky_limits[2], 1.1 * sky_limits[3])

    ax2 = fig.add_axes([0.575, 0.15, 0.4, 0.8])
    ax2.plot(x, y, 'k.', markersize=6)
    ax2.plot(fp_x, fp_y, 'k', linewidth=3)
    ax2.set_xlabel(r'Focal Plane Position $x$ (deg)')
    ax2.set_ylabel(r'Focal Plane Position $y$ (deg)')
    ax2.set_xlim(1.3 * np.min(fp_x), 1.3 * np.max(fp_x))
    ax2.set_ylim(1.3 * np.min(fp_y), 1.3 * np.max(fp_y))

    x = np.linspace(np.min(fp_x), np.max(fp_x), 100)
    y = np.linspace(np.min(fp_y), np.max(fp_y), 100)
    X, Y = np.meshgrid(x, y)
    FoV = [np.max(fp_x) - np.min(fp_x), np.max(fp_y) - np.min(fp_y)]
    true_ff = true_functions.flat_field(X, Y, FoV)
    levels = np.linspace(np.min(true_ff), np.max(true_ff), 25)
    labels = levels[::3]
    cs = ax2.contour(X, Y, true_ff, colors='0.75', levels=levels)
    cs.clabel(labels, fontsize=8)

    filename = string.replace(filename, '.p', plot_suffix)
    fig.savefig(filename)
    fig.clf()
    if verbose:
        print('...done!')


def survey(filename, survey_file, FoV, sky_limits, density, fig_width=8.3,
                                                suffix='.png', verbose=False):
    ''' This function plots the survey strategy, the sky coverage and a
    histogram of the number of times each source is observed.

    Input
    -----
    filename            :   string
        The path to the survey dataset pickle
    survey_file         :   string
        The path to the survey file
    FoV                 :   float array
        The field-of-view of the imager in degrees (dalpha, dbeta)
    sky_limits          :   float_array
        The area of sky to generate sources in
        [alpha_min, alpha_max, beta_min, beta_max]
    density             :   int
        The maximum number of sources (all magnitude) per unit area
        to generate for the self-calibration simulations
    figure_width        :   float
        The width of the figure in inches, default it 8.3
    suffix              :   string
        The format of the saved figure, either '.pdf', '.eps' or '.png'.
        Default is '.png'
    verbose             :   Boolean
        Set to true to run function in verbose mode
    '''
    if verbose:
        print("Plotting Survey Strategy and Coverage from {0} and {1}..."
                    .format(filename, survey_file))

    s = pickle.load(open(filename))
    survey = np.loadtxt(survey_file)

    fig = plt.figure(figsize=(fig_width, fig_width))

    ax1 = fig.add_axes([0.12, 0.55, 0.4, 0.4])

    x_min = -FoV[0] / 2
    y_min = -FoV[1] / 2
    x_max = FoV[0] / 2
    y_max = FoV[1] / 2
    x = np.array([x_min, x_min, x_max, x_max, x_min])
    y = np.array([y_min, y_max, y_max, y_min, y_min])
    for image in survey:
        alpha, beta = transformations.fp2sky(x, y, image[1:3], image[3])
        ax1.plot(alpha, beta, 'k-', alpha=0.5)
    min_sky = np.min(sky_limits) \
                            - 0.1 * (np.max(sky_limits) - np.min(sky_limits))
    max_sky = np.max(sky_limits) \
                            + 0.1 * (np.max(sky_limits) - np.min(sky_limits))
    ax1.set_xlim(min_sky, max_sky)
    ax1.set_ylim(min_sky, max_sky)

    ax2 = fig.add_axes([0.52, 0.55, 0.4, 0.4])
    ax2.set_yticklabels([])
    ax2.get_xticklabels()[0].set_visible(False)

    no_measurements = np.histogram(s.k, bins=range(np.max(s.k)))[0]

    source_ID = np.unique(s.k)
    plot_array = np.zeros((4, source_ID.size))
    plot_array[0, :] = source_ID
    for indx in source_ID:
        found = s.k == indx
        no_obs = sum(found)
        if no_obs > 0:
            plot_array[1, :] = s.alpha[found][0]
            plot_array[2, :] = s.beta[found][0]
            plot_array[3, :] = no_obs
    max_obs = np.max(no_obs)
    for indx in source_ID:
        ax2.plot(plot_array[1, :], plot_array[2, :], '.',
                                                    color=(no_obs / max_obs))

    ax2.set_xlim(min_sky, max_sky)
    ax2.set_ylim(min_sky, max_sky)

    fig.text(0.52, 0.5, r'Sky Position $\alpha$ (deg)', ha='center')
    ax1.set_ylabel(r'Sky Position $\beta$ (deg)')

    ax3 = fig.add_axes([0.12, 0.075, 0.4, 0.4])
    ax3.set_yscale('log')

    ax3.hist(no_measurements, color='k',
                                    bins=range(np.max(no_measurements) + 3))

    ax3.set_xlabel(r'Number of Source Observations')
    ax3.set_ylabel(r'Number of Sources')

    filename = string.replace(filename, '.p', plot_suffix)
    plt.savefig(filename)
    plt.clf()
    if verbose:
        print("...done!")


def variance(filename, fig_width=8.3, suffix='.png', verbose=False):
    ''' This function plots the both the assumed and the *true* measurement
    variances

    Input
    -----
    filename            :   string
        The path to the survey file
    figure_width        :   float
        The width of the figure in inches, default it 8.3
    suffix              :   string
        The format of the saved figure, either '.pdf', '.eps' or '.png'.
        Default is '.png'
    verbose             :   Boolean
        Set to true to run function in verbose mode
    '''
    if verbose:
        print("Plotting Inverse-Variances from {0}...".format(filename))

    s = pickle.load(open(filename))

    fig = plt.figure(figsize=(fig_width, 0.5 * fig_width))
    ax1 = fig.add_axes([0.1, 0.1, 0.43, 0.8])
    ax2 = fig.add_axes([0.53, 0.1, 0.43, 0.8])

    ax1.plot(s.counts, 1 / s.invvar, 'k.', alpha=0.5, markersize=1)
    ax1.text(0.9, 0.1, 'Assumed', bbox=dict(facecolor='w', alpha=0.5),
                            va='center', ha='right', transform=ax1.transAxes)

    ax2.plot(s.counts, 1 / s.true_invvar, 'k.', alpha=0.5, markersize=1)
    ax2.text(0.9, 0.1, 'True', bbox=dict(facecolor='w', alpha=0.5),
                            va='center', ha='right', transform=ax2.transAxes)

    ax1.set_ylim(np.min(0.8 / s.invvar), np.max(1.1 / s.true_invvar))
    ax2.set_ylim(np.min(0.8 / s.invvar), np.max(1.1 / s.true_invvar))

    ax1.set_ylabel(r'$\sigma^2$')
    fig.text(0.53, 0.025, r'Count Rate (s$^{-1}$)', va='center', ha='center')
    ax2.get_xticklabels()[0].set_visible(False)
    ax2.set_yticklabels([])

    variance_filename = string.replace(filename, '.p', '_variance.png')
    fig.savefig(variance_filename)
    fig.clf()
    if verbose:
        print('...done!')


def s2n(filename, fig_width=8.3, suffix='.png', verbose=False):
    ''' This function plots the the assumed and the *true* measurement
    signal-to-noises

    Input
    -----
    filename            :   string
        The path to the survey file
    figure_width        :   float
        The width of the figure in inches, default it 8.3
    suffix              :   string
        The format of the saved figure, either '.pdf', '.eps' or '.png'.
        Default is '.png'
    verbose             :   Boolean
        Set to true to run function in verbose mode
    '''

    if verbose:
        print("Plotting Signal-to-Noise from {0}...".format(filename))

    s = pickle.load(open(filename))

    fig = plt.figure(figsize=(fig_width, 0.5 * fig_width))
    ax1 = fig.add_axes([0.1, 0.1, 0.43, 0.8])
    ax2 = fig.add_axes([0.53, 0.1, 0.43, 0.8])

    s2n = s.counts * np.sqrt(s.invvar)
    true_s2n = s.counts * np.sqrt(s.true_invvar)

    ax1.plot(s.counts, s2n, 'k.', alpha=0.5, markersize=1)
    ax1.text(0.9, 0.1, 'Assumed', bbox=dict(facecolor='w', alpha=0.5),
                            va='center', ha='right', transform=ax1.transAxes)

    ax2.plot(s.counts, true_s2n, 'k.', alpha=0.5, markersize=1)
    ax2.text(0.9, 0.1, 'True', bbox=dict(facecolor='w', alpha=0.5),
                            va='center', ha='right', transform=ax2.transAxes)

    ax1.set_ylim(np.min(0.6 * true_s2n), np.max(1.1 * true_s2n))
    ax2.set_ylim(np.min(0.6 * true_s2n), np.max(1.1 * true_s2n))

    ax1.set_ylabel(r'Signal-to-Noise')
    fig.text(0.53, 0.025, r'Count Rate (s$^{-1}$)', va='center', ha='center')
    ax2.get_xticklabels()[0].set_visible(False)
    ax2.set_yticklabels([])

    s2n_filename = string.replace(filename, '.p', '_s2n.png')
    fig.savefig(s2n_filename)
    fig.clf()

    if verbose:
        print('...done!')


def flat_fields(filename, FoV, ff_samples, best_fit_params, fig_width=8.3,
                                                suffix='.png', verbose=False):
    ''' This function plots the fitted flat-field against the *true* (top) and
    the best-in-basis (bottom) flat-fields. Left: contour plots to compare
    flat-fields, Right: residuals between two flat-fields.

    Input
    -----
    filename            :   string
        The path to the survey files
    FoV                 :   float array
        The size of the imagers field-of-view in degrees (dalpha, dbeta)
    ff_samples          :   float array
        The number of points at which the instrument responses are sampled in
    the (x-direction, y-direction)
        The imager's field-of-view in degrees (dalpha, dbeta)
    best_fit_params     :   numpy array
        Array of the best-in-basis parameters
    figure_width        :   float
        The width of the figure in inches, default it 8.3
    suffix              :   string
        The format of the saved figure, either '.pdf', '.eps' or '.png'.
        Default is '.png'
    verbose             :   Boolean
        Set to true to run function in verbose mode
    '''

    if verbose:
        print("Plotting flat-field from {0}...".format(filename))

    scale = 0.88

    dic = pickle.load(open(filename))

    fig = plt.figure(figsize=(fig_width, scale * fig_width))
    fig.clf()

    ax_cb = fig.add_axes([0.85, 0.075 / scale, 0.05, 0.75 / scale])
    ax1 = fig.add_axes([0.075, 0.075 / scale, 0.375, 0.375 / scale])
    ax2 = fig.add_axes([0.075, 0.45 / scale, 0.375, 0.375 / scale])
    ax3 = fig.add_axes([0.45, 0.075 / scale, 0.375, 0.375 / scale])
    ax4 = fig.add_axes([0.45, 0.45 / scale, 0.375, 0.375 / scale])

    ax2.set_xticks([])
    ax4.set_xticks([])
    ax4.set_yticks([])
    ax3.set_yticks([])
    ax_cb.set_xticks([])

    fig.text(0.45, 0.025 / scale, r'Focal Plane Position $x$ (deg)',
                                                    va='center', ha='center')
    fig.text(0.015, 0.45 / scale, r'Focal Plane Position $y$ (deg)',
                                va='center', ha='center', rotation='vertical')

    true_ff = true_functions.flat_field(dic['x'], dic['y'], FoV)
    best_ff = self_calibration.evaluate_flat_field(dic['x'].flatten(),
                                        dic['y'].flatten(),
                                        best_fit_params).reshape(ff_samples)
    levels = np.linspace(np.min(true_ff), np.max(true_ff), 25)
    labels = levels[::2]
    cs = ax2.contour(dic['x'], dic['y'], true_ff, colors='0.75', levels=levels)
    cs = ax1.contour(dic['x'], dic['y'], best_ff, colors='0.75', levels=levels)
    cs = ax1.contour(dic['x'], dic['y'], dic['fitted_ff'], colors='k',
                                                                levels=levels)
    cs.clabel(labels, fontsize=8)
    cs = ax2.contour(dic['x'], dic['y'], dic['fitted_ff'], colors='k',
                                                                levels=levels)
    cs.clabel(labels, fontsize=8)

    true_residual = 100 * (dic['fitted_ff'] - true_ff) / true_ff
    best_residual = 100 * (dic['fitted_ff'] - best_ff) / best_ff
    fp = (-FoV[0] / 2, FoV[0] / 2, -FoV[1] / 2, FoV[1] / 2)
    a = ax4.imshow(true_residual, extent=fp, vmin=-0.5, vmax=0.5,
                                                    cmap='gray', aspect='auto')
    ax3.imshow(best_residual, extent=fp, vmin=-0.5, vmax=0.5,
                                                    cmap='gray', aspect='auto')

    cbar = fig.colorbar(a, ax_cb, orientation='vertical')
    cbar.set_label(r'Residuals  (%)')

    filename = string.replace(filename, '.p', plot_suffix)
    fig.savefig(filename)
    fig.clf()
    if verbose:
        print('...done!')


if __name__ == "__main__":

    verbose = True
    mult_proc = True
    plot_suffix = ".png"
    figure_width = 8.3

    # Change process nice level to 10 (important if multiprocessing)
    os.nice(10)

    dir_path = get_dir_path()

    old_figures = glob.glob('{0}/*.png'.format(dir_path))
    old_figures += glob.glob('{0}/*.pdf'.format(dir_path))
    for path in old_figures:
        os.remove(path)
    params = pickle.load(open('{0}/parameters.p'.format(dir_path)))

    source_catalog_files = glob.glob('{0}/source_catalog.p'.format(dir_path))
    for path in source_catalog_files:
        source_catalog(path, params['sky_limits'],
                                params['density_of_stars'],
                                params['m_min'], params['m_max'],
                                params['powerlaw_constants'],
                                fig_width=8.3,
                                suffix=plot_suffix,
                                verbose=verbose)

    files = glob.glob('{0}/measurement_catalog.p'.format(dir_path))
    for path in files:
        survey(path, params['survey_file'], params['FoV'],
                        params['sky_limits'], params['density_of_stars'],
                        fig_width=8.3, suffix=plot_suffix, verbose=verbose)
        variance(path, fig_width=8.3, suffix=plot_suffix, verbose=verbose)
        s2n(path, fig_width=8.3, suffix=plot_suffix, verbose=verbose)

    files = glob.glob('{0}/camera_image.p'.format(dir_path))
    for path in files:
        camera_image(path,  params['sky_limits'], params['density_of_stars'],
        fig_width=8.3, suffix=plot_suffix, verbose=verbose)

    files = glob.glob('{0}/FF/*_ff.p'.format(dir_path))
    for path in files:
        flat_fields(path, params['FoV'], params['ff_samples'],
                            params['best_fit_params'], fig_width=8.3,
                            suffix=plot_suffix, verbose=verbose)
