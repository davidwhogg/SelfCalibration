#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pylab as plt
import pickle
import os

plt.figure()
fontsize = 20
tick_fontsize = 12
guide_frac = 0.9 # what coverage to draw guide line at

# Select current or old pickle file
filename = "coverage.p"
strategies_filename = "./strategies.p"
if os.path.exists(strategies_filename) != True:
  filename += '_old'
  strategies_filename += '_old'
  plt.title('(Old Data)')

sky_surveys = pickle.load(open(strategies_filename))

for surveys in sky_surveys:
  dic = pickle.load(open(('./'+surveys+'_'+filename)))
  number_stars = dic['number_stars']
  k = dic['k']
  strategy = dic['strategy']

  num_obs_of_star = np.bincount(k)
  num_star_with_N_obs = np.bincount(num_obs_of_star)
  # "Integrate" 
  for i in range(1,len(num_star_with_N_obs)):
    num_star_with_N_obs[i] += num_star_with_N_obs[i-1]
    
  fraction = (1.*num_star_with_N_obs)/number_stars

  plt.plot(fraction, label = ('Strategy '+strategy))
  gt90 = np.where (fraction > guide_frac)
  print gt90[0][0]
  guide_line = np.array([[0,float(gt90[0][0]), float(gt90[0][0])], [fraction[gt90[0][0]], fraction[gt90[0][0]],0.]])
  #plt.plot(guide_line[0], guide_line[1],'k--', alpha = 0.25)
  
plt.xlabel(r"Number of Observations")
plt.ylabel(r"Fraction of Sources Covered")
plt.title('Survey Coverage')
plt.legend(loc=4)
filename = "./coverage.png"
ax = plt.gca()
for tick in ax.xaxis.get_major_ticks():
  tick.label1.set_fontsize(tick_fontsize)
for tick in ax.yaxis.get_major_ticks():
  tick.label1.set_fontsize(tick_fontsize)
plt.ylim(0.,1.1)
plt.savefig(filename,bbox_inches='tight')#,pad_inches=0.5)
plt.clf()
