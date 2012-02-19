#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Make Python 3 compatible
from __future__ import division, print_function

import default_parameters
import simulation
import save_out
import numpy as np
import god

params = default_parameters.dic
files = []
dic = {'params': params, 'plotdata': True, 'verbosemode': False}

                                                
# Create Survey Files
pointings = 1000
samples = 10
bb_badness = np.zeros(samples)
survey_size = np.linspace(0.35, 8, samples)
survey = np.zeros((pointings, 4))
survey[:,0] = np.arange(pointings)

dirs = []
for indx in range(samples):
    dirs.append(("test_" + str(indx)))


for indx in range(samples):
    survey[:,1] = survey_size[indx] * (np.random.random(pointings) - 0.5)
    survey[:,2] = survey_size[indx] * (np.random.random(pointings) - 0.5)
    survey[:,3] = 360. * (np.random.random(pointings))
    np.savetxt(("random_" + str(indx)), survey)
    dic['data_dir'] = dirs[indx]
    dic['survey_file'] = 'random_' + str(indx)
    simulation.run_sim(dic)
    filename = dirs[indx] + '/result'
    bb_badness[indx] = np.loadtxt(filename)[2]

import matplotlib.pylab as plt
plt.plot(survey_size/0.7, bb_badness)
plt.xlabel("Deep Survey Field Size (focal plane size)")
plt.ylabel("Best-in-Basis Badness")
plt.savefig("test.png")
