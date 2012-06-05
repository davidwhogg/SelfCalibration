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

sys.path.append('./../..')  # add simulator modules to Python path
import transformations
from plot import camera_image, flat_fields

GOLDEN_RATIO = 1.61803399

def plot_grid(xlabel, ylabel):
    fig = plt.figure(figsize=(figure_width, figure_width))
    fig.clf()

    ax1 = fig.add_axes([0.075, 0.5, 0.425, 0.425])
    ax2 = fig.add_axes([0.5, 0.5, 0.425, 0.425])    
    ax3 = fig.add_axes([0.075, 0.075, 0.425, 0.425])
    ax4 = fig.add_axes([0.5, 0.075, 0.425, 0.425])
    
    ax1.set_xticklabels([''])
    ax2.set_xticklabels([''])
    ax2.set_yticklabels([''])
    ax4.set_yticklabels([''])
    
    fig.text(0.5, 0.025, xlabel, va='center', ha='center')
    fig.text(0.025, 0.5, ylabel, va='center', ha='center', rotation=90)
    
    return [fig, [ax1, ax2, ax3, ax4]]

def paper_surveys(dirs):
    xlabel = r'Sky Position $\alpha$ (deg$^2$)'
    ylabel = r'Sky Position $\beta$ (deg$^2$)'
    fig, axes = plot_grid(xlabel, ylabel)
    for indx in range(len(dirs)):
        params = pickle.load(open('{0}/parameters.p'.format(dirs[indx])))
        single_survey(fig, axes[indx], (dirs[indx] + '/survey.txt'), 
                            params['FoV'], params['sky_limits'],
                            fig_width=8.3, suffix='.png', verbose=False)
       
    fig.savefig('simple_surveys.pdf')
    fig.clf() 

def single_survey(fig, ax, survey_file, FoV, sky_limits, fig_width=8.3,
                                                suffix='.png', verbose=False):

    if verbose:
        print("Plotting Survey Strategy from {0}...".format(survey))
    
    x_min = -FoV[0] / 2
    y_min = -FoV[1] / 2
    x_max = FoV[0] / 2
    y_max = FoV[1] / 2
    x = np.array([x_min, x_min, x_max, x_max, x_min])
    y = np.array([y_min, y_max, y_max, y_min, y_min])
    
    survey = np.loadtxt(survey_file)
    for image in survey:
        alpha, beta = transformations.fp2sky(x, y, image[1:3], image[3])
        ax.plot(alpha, beta, 'k-', alpha=0.5)
    min_sky = np.min(sky_limits) \
                            - 0.1 * (np.max(sky_limits) - np.min(sky_limits))
    max_sky = np.max(sky_limits) \
                            + 0.1 * (np.max(sky_limits) - np.min(sky_limits))
    ax.set_xlim(min_sky, max_sky)
    ax.set_ylim(min_sky, max_sky)
    
    if verbose:
        print("...done!")


def paper_camera(dirs, fig_width=8.3, suffix='png', verbose=False):
    files = []
    for dir_path in dirs:
        files += glob.glob('{0}/camera_image.p'.format(dir_path))
    camera_file = files[0]
    parameter_file = camera_file.replace('camera_image.p', 'parameters.p')
    params = pickle.load(open(parameter_file))
    camera_image(camera_file,  params['sky_limits'],
                        params['density_of_stars'],
                        output_path='camera_image',
                        fig_width=fig_width,
                        suffix=suffix,
                        verbose=verbose)


def paper_flat_fields(dirs, fig_width=8.3, suffix='png', verbose=False):
    
    for directory in dirs:
        files = glob.glob('{0}/FF/*_ff.p'.format(directory))
        
        count = 0
        final_solution = ''
        for ff in files:
            print(ff)
            iteration = int(ff.split('/FF/')[1].split('_ff.p')[0])
            if iteration > count:
                final_solution = ff
                count = iteration
                
        filename = directory + '_' + str(count) + '_ff'
        params = pickle.load(open('{0}/parameters.p'.format(directory)))
        flat_fields(final_solution, params['FoV'], params['ff_samples'],
                            params['best_fit_params'], output_path=filename,
                            fig_width=8.3, suffix=plot_suffix, verbose=verbose)

if __name__ == "__main__":

    old_figures = glob.glob('*.png')
    old_figures += glob.glob('*.pdf')
    for path in old_figures:
        os.remove(path)

    dirs = ['A', 'B', 'C', 'D']
    plot_suffix = ".png"
    figure_width = 8.3
    verbose = True

    if len(dirs) != 4:
        print("Error: {0} directories specified, only 4 are allowed!")
        sys.exit(1)

    paper_surveys(dirs)
    paper_camera(dirs)
    paper_flat_fields(dirs)
