#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pylab as plt
import pickle
import os


plt.figure()
fontsize = 20
tick_fontsize = 12

# Select current or old pickle file
filename = "./invvar.p"
if os.path.exists(filename) != True:
  filename += '_old'
  plt.title('(Old Data)')

# Load data from pickle file
dic = pickle.load(open(filename))
counts = dic['counts']
true_invvar = dic['true_invvar']
reported_invvar = dic['reported_invvar']

# Sort true invvar by counts so we can plot as a line
sort = np.zeros((len(counts), 2))
sort[:,0] = counts[:]
sort[:,1] = reported_invvar[:]
sort = sort[sort[:,0].argsort(),:]
sort_reported_invvar = sort[:,1]
sort_counts = sort[:,0]

plt.plot(counts, np.log10((1/true_invvar)/counts**2),'r.', markersize = 1., label = "Actual Variance")
plt.plot(sort_counts, np.log10((1/sort_reported_invvar)/sort_counts**2),'k', label = "Reported Variance")
plt.xlabel(r'$c_i$', fontsize=fontsize)
plt.ylabel(ur'$\log_{10}(\frac{{\sigma_i}^2}{c_i^2})$',fontsize=fontsize)#, rotation='horizontal')
plt.legend()
ax = plt.gca()
for tick in ax.xaxis.get_major_ticks():
  tick.label1.set_fontsize(tick_fontsize)
for tick in ax.yaxis.get_major_ticks():
  tick.label1.set_fontsize(tick_fontsize)


filename = './invvar_plot.png'
plt.savefig(filename,bbox_inches='tight',pad_inches=0.)  