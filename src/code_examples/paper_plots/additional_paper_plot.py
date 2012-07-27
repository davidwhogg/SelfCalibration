#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Rory Holmes
# Feb 2012

# This script contains the 2x2 survey comparison figure from the paper. All
# other figures are made with the standard plot.py script

# Make Python 3 compatible
from __future__ import division, print_function

import matplotlib.pyplot as plt
import glob
import pickle
import numpy as np
import sys
import os

sys.path.append('./../..')  # add simulator modules to Python path
import transformations

old_figures = glob.glob('*.png')
old_figures += glob.glob('*.pdf')
for path in old_figures:
    os.remove(path)

dirs = ['A', 'B', 'C', 'D']
plot_suffix = ".pdf"
verbose = True

plt.rcParams['font.family'] = 'Computer Modern'
plt.rcParams['text.usetex'] = True

fig = plt.figure(figsize=(7, 7))

middle = [0.525, 0.525]
size = [0.45, 0.45]
ax = [
    fig.add_axes([middle[0] - size[0], middle[1], size[0], size[1]]),
    fig.add_axes([middle[0], middle[1], size[0], size[1]]),
    fig.add_axes([middle[0] - size[0], middle[1] - size[1], size[0], size[1]]),
    fig.add_axes([middle[0], middle[1] - size[1], size[0], size[1]])
    ]

ax[0].set_xticklabels([''])
ax[1].set_yticklabels([''])
ax[3].set_yticklabels([''])

fig.text(0.5, 0.025, r'Sky Position $\alpha$ (deg$^2$)',
                                    va='center', ha='center')
fig.text(0.025, 0.5, r'Sky Position $\beta$ (deg$^2$)',
                                    va='center', ha='center', rotation=90)

for indx in range(len(dirs)):

    params = pickle.load(open('{0}/parameters.p'.format(dirs[indx])))
    FoV = params['FoV']
    sky_limits = params['sky_limits']

    x = 0.5 * np.array([-FoV[0], -FoV[0], FoV[0], FoV[0], -FoV[0]])
    y = 0.5 * np.array([-FoV[1], FoV[1], FoV[1], -FoV[1], -FoV[1]])

    survey = np.loadtxt((dirs[indx] + '/survey.txt'))

    for image in survey:
        alpha, beta = transformations.fp2sky(x, y, image[1:3], image[3])
        ax[indx].plot(alpha, beta, 'k-', alpha=0.2)

    min_sky = np.min(sky_limits) \
                        - 0.02 * (np.max(sky_limits) - np.min(sky_limits))
    max_sky = np.max(sky_limits) \
                        + 0.02 * (np.max(sky_limits) - np.min(sky_limits))
    ax[indx].set_xlim(min_sky, max_sky)
    ax[indx].set_ylim(min_sky, max_sky)

ax[0].text(0.95, 0.04, 'A', va='center', ha='center',
                                                    transform=ax[0].transAxes)
ax[1].text(0.05, 0.04, 'B', va='center', ha='center',
                                                    transform=ax[1].transAxes)
ax[2].text(0.95, 0.96, 'C', va='center', ha='center',
                                                    transform=ax[2].transAxes)
ax[3].text(0.05, 0.96, 'D', va='center', ha='center',
                                                    transform=ax[3].transAxes)

fig.savefig('simple_surveys' + plot_suffix)
