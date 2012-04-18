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

mult_proc = True
plot_suffix = ".pdf"

fig_width_pt = 469.75499  # get from LaTeX Document with \showthe\columnwidth
inches_per_pt = 1.0 / 72.27  # Convert pt to inches
fig_width = fig_width_pt * inches_per_pt
fig_height = (np.sqrt(5) - 1.0) / 2.0 * fig_width  # golden ratio

os.nice(10)  # Change process nice level to 10 (important if multiprocessing)


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


def plot_survey(filename, sky_limits, density, m_min, m_max, A, verbose=False):
    ''' This function plots a sky catalog, both spatially and as a function of
    magnitude

    Input
    -----
    filename            :   string
        The path to the survey files
    sky_limits          :   float array
        The area of sky to generate sources in
        [alpha_min, alpha_max, beta_min, beta_max]
    density             :   int
        The maximum number of sources (all magnitude) per unit area
        to generate for the self-calibration simulations
    m_min               :   float
        The saturation limit of the simulated imager
    m_max               :   float
        The 10-sigma detection limit of the simulated imager
    A                   :   float array
        The parameters describing the magnitude distribution of the sources
        in the sky, according to: log10(dN/dm) = A + B * mag + C * mag ** 2
    verbose             :   Boolean
        Set to true to run function in verbose mode
    '''

    if verbose:
        print("Plotting Survey from {0}...".format(filename))

    s = pickle.load(open(filename))

    fig = plt.figure(figsize=(8.3, 4.15))

    ax1 = fig.add_axes([0.075, 0.15, 0.4, 0.8])
    ax1.plot(s.alpha, s.beta, 'k.', alpha=(10 / density))
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


if __name__ == "__main__":
    dir_path = get_dir_path()

    params = pickle.load(open('{0}/parameters.p'.format(dir_path)))

    files = glob.glob('{0}/source_catalog.p'.format(dir_path))
    for paths in files:
        plot_survey(paths, params['sky_limits'], params['density_of_stars'],
                        params['m_min'], params['m_max'],
                        params['powerlaw_constants'], params['verbose'])
