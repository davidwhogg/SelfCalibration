#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Rory Holmes
# Feb 2012

# This file contains the functions to create the figures from the
# self-calibration simulations.

# Make Python 3 compatible
from __future__ import division, print_function

import matplotlib.pyplot as plt
import matplotlib
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
    ax1.set_xlim(sky_limits[0], sky_limits[1])
    ax1.set_ylim(sky_limits[2], sky_limits[3])

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

    filename = string.replace(filename, '.p', suffix)
    plt.savefig(filename)
    plt.clf()
    if verbose:
        print("...done!")


def camera_image(filename, sky_limits, density, output_path=False,
                        fig_width=8.3, suffix='.png', verbose=False):
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
    output_path         :   string
        The path to save the figure to. If no value is given, figure is saved
        in same directory as source pickle
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
    ax1.set_xlim(sky_limits[0], sky_limits[1])
    ax1.set_ylim(sky_limits[2], sky_limits[3])

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

    if output_path:
        filename = output_path + suffix
    else:
        filename = string.replace(filename, '.p', suffix)
    fig.savefig(filename)
    fig.clf()
    if verbose:
        print('...done!')


def survey_coverage(ax, source_filename, measurement_filename,
                                zero_level='0.9', nobs_plot_lim=None, ms=5):
    ''' This function plots the sources on the sky, with a color-coding
    corresponding to the number of times each source was observed.

    Input
    -----
    ax                      :       matplotlib axes instance
        The axes instance to plot the coverage map to
    source_filename         :   string
        The path to the fitted source catalog
    measurement_filename    :   string
        The path to the measurement file
    zero_level              :   string of a float
        The matplotlib color used to plot sources with zero observations
    nobs_plot_lim           :       int
        The maximum number of observations to plot (used homogenize colormaps
        across multiple plots)
    ms                      :       int
        The size of the source markers

    Return
    ------
    out                     :       list
        A two element list. list[0] contains the unique colors used in the plot
        and list[1] contains the corresponding number of observations (for
        colorbar generation)
    '''

    sky_cat = pickle.load(open(source_filename))
    ax.plot(sky_cat.alpha, sky_cat.beta, '.', color=zero_level,
                                    markersize=ms, mec='none', zorder=-2000)

    measurement_catalog = pickle.load(open(measurement_filename))
    source_ID = np.unique(measurement_catalog.k)
    x = np.zeros(source_ID.size)
    y = x.copy()
    nobs = x.copy()
    for i, sid in enumerate(source_ID):
        found = np.where(measurement_catalog.k == sid)[0]
        nobs[i] = found.size
        if nobs[i] > 0:
            x[i] = measurement_catalog.alpha[found[0]]
            y[i] = measurement_catalog.beta[found[0]]

    if nobs_plot_lim == None:
        nobs_plot_lim = np.max(nobs[okay])

    color = np.array(['{0:4.2f}'.format(cn) for cn in
               float(zero_level) * (1 - nobs / nobs_plot_lim)])

    colors = np.unique(color)
    for c in colors:
        found = np.where(color == c)[0]
        if (found.size > 0):
            zorder = -1000 * float(c)
            ax.plot(x[found], y[found], '.', color=c,
                                    mec='none', markersize=ms, zorder=zorder)

    # Provide a colorbar_nobs vector that can be used with the colors vector
    # to define the colorbar
    colorbar_nobs = np.zeros(colors.size)
    for indx in range(colors.size):
        found = np.where(color == colors[indx])[0]
        colorbar_nobs[indx] = nobs[found[0]]

    # Append low and high limits to color range
    colors = np.append('0.00', colors)
    colorbar_nobs = np.append(colorbar_nobs, 0)
    colors = np.append(colors, zero_level)
    colorbar_nobs = np.append(nobs_plot_lim, colorbar_nobs)

    ax.xaxis.set_zorder(3000)
    ax.yaxis.set_zorder(3000)

    return [colors, colorbar_nobs]


def survey_coverage_colorbar(ax_cb, colors, colorbar_nobs):
    ''' This function plots a colorbar for the survey_coverage plot

    Input
    -----
    ax_cb                   :       matplotlib axes instance
        The axes instance to plot the colorbar to
    colors                  :       numpy string array
        Array containing the unique numbers of source observations
    colorbar_nobs           :       numpy array
        The number of observations corresponding to the colors in the colors
        vector
    '''

    assert colors.size == colorbar_nobs.size

    obs_intrp = np.linspace(colorbar_nobs.min(), colorbar_nobs.max(), 100)
    colors_intrp = np.linspace(np.max(colors.astype(float)),
                                            np.min(colors.astype(float)), 100)

    for indx in range(obs_intrp.size - 1):
        ax_cb.fill_between([0, 1], [obs_intrp[indx], obs_intrp[indx]],
                            y2=[obs_intrp[indx + 1], obs_intrp[indx + 1]],
                            color=str(colors_intrp[indx]))

    ax_cb.set_xlim(0, 1)
    ax_cb.set_ylim(0, obs_intrp.max())

    return None


def survey_coverage_histogram(ax, measurement_filename, xlim, ylim=None):
    ''' This function plots a histogram of the number of source observations

    Input
    -----
    ax                      :       matplotlib axes instance
        The axes instance to plot the histogram to
    measurement_filename    :       string
        The path to the measurement file
    xlim                    :       list
        The limits for the x-axis [xmin, xmax].
    xlim                    :       list
        The limits for the x-axis [ymin, ymax].
    '''

    measurement_catalog = pickle.load(open(measurement_filename))
    k = measurement_catalog.k[measurement_catalog.is_in_analysis_region]
    source_ID = np.unique(k)
    nobs = np.zeros(source_ID.size)

    for i, sid in enumerate(source_ID):
        found = np.where(k == sid)[0]
        nobs[i] = found.size

    assert np.sum(nobs).astype(int) == k.size
    hist = ax.hist(nobs, color='k',
                        bins=np.arange(int(xlim[1]) + 5) - 0.5,
                        histtype='step')

    
    ax.set_xlim(xlim)
    if ylim is None:
        ax.set_ylim(-0.1 * np.max(hist[0]), 1.1 * np.max(hist[0]))
    else:
        ax.set_ylim(ylim)

    return None


def survey_source_error(ax, source_filename, fitted_filename, ms=10,
                            lower_limit=None, upper_limit=None,
                            lower_color=0.9, upper_color=0.):
    ''' This function plots a source map on the given matplotlib axes instance,
    with the sources color-coded according to their error in the final fitted
    solution.

    Input
    -----
    ax                      :       matplotlib axes instance
        The axes instance to plot the error map to
    source_filename         :       string
        The path to the source catalog, which contains the *true* source fluxes
    fitted_filename         :       string
        The path to the fitted source catalog
    ms                      :       int
        The markersize for each source
    lims                    :       list
        The lower [0] and upper [1] error limit for the plot
    lower_color             :       float
        The matplotlib color value for the sources with the lowest errors
    upper_color             :       float
        The matplotlib color value for the sources with the highest errors

    Returns
    -------
    out                     :       list
        A two element list. list[0] contains the unique colors used in the plot
        and list[1] contains the corresponding errors (for colorbar generation)
    '''

    true = pickle.load(open(source_filename))
    fitted = pickle.load(open(fitted_filename))

    x = np.zeros(len(fitted.k))
    y = x.copy()
    error = x.copy()

    for i, sid in enumerate(fitted.k):
        found = np.where(fitted.k == sid)[0]
        if len(found) > 0:
            x[i] = fitted.alpha[found[0]]
            y[i] = fitted.beta[found[0]]
            original = true.flux[np.where(true.k == sid)[0]]
            error[i] = (100 * (np.average(fitted.flux[found]) - original)
                                                                    / original)

    if lower_limit is None:
        lower_limit = np.min(error)
    if upper_limit is None:
        upper_limit = np.max(error)

    error = np.clip(error, lower_limit, upper_limit)

    color = np.array(['{0:4.2f}'.format(cn) for cn in
                        0.5 * (upper_color + lower_color)
                        + (upper_color - lower_color)
                        * error / (upper_limit - lower_limit)])

    colors = np.unique(color)
    for c in colors:
        found = np.where(color == c)[0]
        if (found.size > 0):
            zorder = 1000 * np.abs(float(c) - 0.5)
            ax.plot(x[found], y[found], '.', color=c, mec='none',
                                                markersize=ms, zorder=zorder)

    colorbar_errors = np.zeros(colors.size)
    # Provide a error_check vector that can be used to assert that the colormap
    # in the colorbar is correct
    for indx in range(colors.size):
        found = np.where(color == colors[indx])[0]
        colorbar_errors[indx] = error[found[0]]

    # Add end points for colorbar
    colorbar_errors = np.append(colorbar_errors, lower_limit)
    colorbar_errors = np.append(upper_limit, colorbar_errors)
    colors = np.append(colors, lower_color)
    colors = np.append(upper_color, colors)

    return [colors, colorbar_errors]


def survey_source_error_colorbar(ax_cb, colors, errors):
    ''' This function plots a colorbar for the survey_source_error plot

    Input
    -----
    ax_cb                   :       matplotlib axes instance
        The axes instance to plot the colorbar to
    colors                  :       numpy string array
        Array containing the unique numbers of source observations
    errors                  :       numpy array
        The errors corresponding to the colors in the colors
        vector
    '''

    assert colors.size == errors.size

    errors_intrp = np.linspace(errors.min(), errors.max(), 100)
    colors_intrp = np.linspace(np.max(colors.astype(float)),
                                            np.min(colors.astype(float)), 100)

    for indx in range(errors_intrp.size - 1):
        ax_cb.fill_between([0, 1], [errors_intrp[indx], errors_intrp[indx]],
                        y2=[errors_intrp[indx + 1], errors_intrp[indx + 1]],
                        color=str(colors_intrp[indx]))

    ax_cb.set_xlim(0, 1)
    ax_cb.set_ylim(errors_intrp.min(), errors_intrp.max())

    return None


def survey_footprint(ax, survey_filename, FoV):
    ''' This function plots the survey footprint map on the given matplotlib
    axes instance

    Input
    -----
    ax                      :       matplotlib axes instance
        The axes instance to plot the footprint map to
    survey_filename         :       string
        The path to the survey file

    '''

    survey = np.loadtxt(survey_filename)

    x_min = -FoV[0] / 2
    y_min = -FoV[1] / 2
    x_max = FoV[0] / 2
    y_max = FoV[1] / 2
    x = np.array([x_min, x_min, x_max, x_max, x_min])
    y = np.array([y_min, y_max, y_max, y_min, y_min])

    for image in survey:
        alpha, beta = transformations.fp2sky(x, y, image[1:3], image[3])
        ax.plot(alpha, beta, 'k-', alpha=0.5)

    return None


def survey(source_filename, measurement_filename,
                        fitted_filename, survey_filename,
                        out_filename, FoV, sky_limits, density,
                        nobs_plot_lim=50,
                        error_lim=[-3., 3.],
                        fig_width=8.3, suffix='.png', verbose=False):
    ''' This function plots the survey footprint, the sky coverage, a
    histogram of the number of times each source is observed and a map of the
    errors in the fitted sources.

    Input
    -----
    source_filename             :   string
        The path to the source catalog
    measurement_filename        :   string
        The path to the measurement catalog
    fitted_filename             :   string
        The path to the fitted source catalog
    survey_filename             :   string
        The path to the survey file
    out_filename                :   string
        The path to save the figure to
    FoV                     :   float array
        The field-of-view of the imager in degrees (dalpha, dbeta)
    sky_limits                  :   float_array
        The area of sky to generate sources in
        [alpha_min, alpha_max, beta_min, beta_max]
    density                     :   int
        The maximum number of sources (all magnitude) per unit area
        to generate for the self-calibration simulations
    nobs_plot_lim               :   int
        The maximum number of source observations to plot
    error_lim                   :   list
        The lower [0] and the upper [1] error limits for the error subplot
    figure_width                :   float
        The width of the figure in inches, default it 8.3
    suffix                      :   string
        The format of the saved figure, either '.pdf', '.eps' or '.png'.
        Default is '.png'
    verbose                     :   Boolean
        Set to true to run function in verbose mode
    '''

    if verbose:
        print("Generating survey plot...")

    middle = [0.5, 0.5]
    size = [0.40, 0.40]
    fig = plt.figure(figsize=(fig_width, fig_width))

    ax1 = fig.add_axes([middle[0] - size[0], middle[1], size[0], size[1]])
    ax1.set_xlim(sky_limits[0], sky_limits[1])
    ax1.set_ylim(sky_limits[2], sky_limits[3])
    ax1.tick_params(labelleft=True, labelright=False,
                                            labeltop=True, labelbottom=False)

    ax2 = fig.add_axes([middle[0], middle[1], size[0], size[1]])
    ax2.set_xlim(sky_limits[0], sky_limits[1])
    ax2.set_ylim(sky_limits[2], sky_limits[3])
    ax2.tick_params(labelleft=False, labelright=False,
                                            labeltop=True, labelbottom=False)
    ax_cb2 = fig.add_axes([middle[0] + 1.03 * size[0],
                                        middle[1] + 0.05 * size[1],
                                        0.07 * size[0],
                                        0.9 * size[1]])
    ax_cb2.tick_params(labelleft=False, labelright=True,
                                            labeltop=False, labelbottom=False)
    ax_cb2.set_xticks([])

    ax3 = fig.add_axes([middle[0] - size[0], middle[1] - size[1],
                                            0.9 * size[0], 0.9 * size[1]])
    ax3.set_ylabel('Number of Sources')
    ax3.set_xlabel('Number of Observations')

    ax4 = fig.add_axes([middle[0], middle[1] - size[1], size[0], size[1]])
    ax4.set_xlim(sky_limits[0], sky_limits[1])
    ax4.set_ylim(sky_limits[2], sky_limits[3])
    ax4.set_xlabel(r'Sky Position $\alpha$ (deg)')

    ax_cb4 = fig.add_axes([middle[0] + 1.03 * size[0],
                                        middle[1] - 0.95 * size[1],
                                        0.07 * size[0],
                                        0.9 * size[1]])

    ax_cb4.tick_params(labelleft=False, labelright=True,
                                            labeltop=False, labelbottom=False)
    ax_cb4.set_xticks([])

    ax1.text(0.96, 0.04, '(a)', va='center', ha='center',
                                        zorder=2000, transform=ax1.transAxes)
    ax2.text(0.04, 0.04, '(b)', va='center', ha='center',
                                        zorder=2000, transform=ax2.transAxes)
    ax3.text(0.95, 0.96, '(c)', va='center', ha='center',
                                        zorder=2000, transform=ax3.transAxes)
    ax4.text(0.04, 0.96, '(d)', va='center', ha='center',
                                        zorder=2000, transform=ax4.transAxes)

    fig.text(middle[0], middle[1] + 1.1 * size[1],
                    r'Sky Position $\alpha$ (deg)',
                    ha='center', va='center')
    fig.text(middle[0] - 1.1 * size[0], middle[1] + 0.5 * size[1],
                    r'Sky Position $\beta$ (deg)',
                    ha='center', va='center', rotation=90)
    fig.text(middle[0] + 1.22 * size[0], middle[1] + 0.5 * size[1],
                    r'Number of Observations',
                    va='center', ha='center', rotation=90)
    fig.text(middle[0] + 1.22 * size[0], middle[1] - 0.5 * size[1],
                    r'Flux Error (percent)',
                    va='center', ha='center', rotation=90)

    if verbose:
        print("...plotting survey footprint map from {0}..."
                                                    .format(survey_filename))
    survey_footprint(ax1, survey_filename, FoV)
    if verbose:
        print('...done...')

    if verbose:
        print('...plotting coverage map from {0}...'
                                                .format(measurement_filename))
    colors, colorbar_nobs = survey_coverage(ax2, source_filename,
                                                measurement_filename,
                                                nobs_plot_lim=nobs_plot_lim,
                                                zero_level='0.9', ms=6)
    survey_coverage_colorbar(ax_cb2, colors, colorbar_nobs)
    if verbose:
        print('...done...')

    if verbose:
        print('Plotting measurement histogram from {0}...'
                                                .format(measurement_filename))
    survey_coverage_histogram(ax3, measurement_filename, (0, nobs_plot_lim))
    if verbose:
        print('...done...')

    if verbose:
        print('Plotting source error map from {0} and {1}...'
                    .format(source_filename, fitted_filename))
    colors, errors = survey_source_error(ax4, source_filename, fitted_filename,
                            lower_limit=error_lim[0], upper_limit=error_lim[1],
                            lower_color=0.9, upper_color=0., ms=6)
    survey_source_error_colorbar(ax_cb4, colors, errors)
    if verbose:
        print('...done...')

    plt.savefig(out_filename + suffix)
    plt.clf()
    if verbose:
        print("...figure done!")

    return None


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


def flat_fields(filename, FoV, ff_samples, best_fit_params, output_path=False,
                            fig_width=8.3, suffix='.png', verbose=False):
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
    output_path         :   string
        The path to save the figure to. If no value is given, figure is saved
        in same directory as source pickle
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

    ax1.text(0.95, 0.96, '(c)', va='center', ha='center',
                        bbox=dict(facecolor='w', edgecolor='w', alpha=0.7),
                        zorder=2000, transform=ax1.transAxes)
    ax2.text(0.95, 0.04, '(a)', va='center', ha='center',
                        bbox=dict(facecolor='w', edgecolor='w', alpha=0.7),
                        zorder=2000, transform=ax2.transAxes)
    ax3.text(0.05, 0.96, '(d)', va='center', ha='center',
                        zorder=2000, transform=ax3.transAxes)
    ax4.text(0.05, 0.04, '(b)', va='center', ha='center',
                        zorder=2000, transform=ax4.transAxes)

    ax2.set_xticklabels([])
    ax4.set_xticklabels([])
    ax4.set_yticklabels([])
    ax3.set_yticklabels([])
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
    a = ax4.imshow(true_residual, extent=fp, vmin=-2., vmax=2.,
                                                    cmap='gray', aspect='auto')
    ax3.imshow(best_residual, extent=fp, vmin=-2., vmax=2.,
                                                    cmap='gray', aspect='auto')

    cbar = fig.colorbar(a, ax_cb, orientation='vertical')
    cbar.set_label(r'Residuals  (percent)')

    if output_path:
        filename = output_path + suffix
    else:
        filename = string.replace(filename, '.p', suffix)
    fig.savefig(filename)
    fig.clf()
    if verbose:
        print('...done!')


if __name__ == "__main__":

    verbose = True
    mult_proc = True
    plot_suffix = ".png"
    figure_width = 7.

    plt.rcParams['font.family'] = 'Computer Modern'
    plt.rcParams['text.usetex'] = True

    # Change process nice level to 10 (important if multiprocessing)
    os.nice(10)

    dir_path = get_dir_path()

    old_figures = glob.glob('{0}/*.png'.format(dir_path))
    old_figures += glob.glob('{0}/*.pdf'.format(dir_path))
    for path in old_figures:
        os.remove(path)
    params = pickle.load(open('{0}/parameters.p'.format(dir_path)))

    source_catalog('{0}/source_catalog.p'.format(dir_path),
                                    params['sky_limits'],
                                    params['density_of_stars'],
                                    params['m_min'],
                                    params['m_max'],
                                    params['powerlaw_constants'],
                                    fig_width=8.3,
                                    suffix=plot_suffix,
                                    verbose=verbose)

    out_filename = '{0}/{1}'.format(dir_path, 'survey')
    survey('{0}/source_catalog.p'.format(dir_path),
                                '{0}/measurement_catalog.p'.format(dir_path),
                                '{0}/fitted_catalog.p'.format(dir_path),
                                '{0}/survey.txt'.format(dir_path),
                                out_filename,
                                params['FoV'],
                                params['sky_limits'],
                                params['density_of_stars'],
                                fig_width=8.3,
                                suffix=plot_suffix,
                                verbose=verbose)

    variance('{0}/measurement_catalog.p'.format(dir_path),
                                    fig_width=8.3,
                                    suffix=plot_suffix,
                                    verbose=verbose)

    s2n('{0}/measurement_catalog.p'.format(dir_path),
                                    fig_width=8.3,
                                    suffix=plot_suffix,
                                    verbose=verbose)

    if os.path.isfile('{0}/camera_image.p'.format(dir_path)):
        camera_image('{0}/camera_image.p'.format(dir_path),
                                        params['sky_limits'],
                                        params['density_of_stars'],
                                        fig_width=8.3,
                                        suffix=plot_suffix,
                                        verbose=verbose)

    files = glob.glob('{0}/FF/*_ff.p'.format(dir_path))
    for path in files:
        flat_fields(path, params['FoV'], params['ff_samples'],
                                    params['best_fit_params'],
                                    fig_width=8.3,
                                    suffix=plot_suffix,
                                    verbose=verbose)
