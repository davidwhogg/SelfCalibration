#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pylab as plt
import pickle
import os


plt.figure(figsize=(13,6))
fontsize = 20
tick_fontsize = 12

# Select current or old pickle file
filename = "./camera_image.p"
if os.path.exists(filename) != True:
  filename += '_old'
  plt.title('(Old Data)')

# Load data from pickle file
dic = pickle.load(open(filename))
x = dic['measured_catalog.x'] 
y = dic['measured_catalog.y']
alpha = dic['sky_catalog.alpha'] 
beta = dic['sky_catalog.beta']
sky = dic['sky'] 
pointing = dic['pointing']
orientation = dic['orientation']
fp_x = dic['fp_x']
fp_y = dic['fp_y']
fp_alpha = dic['fp_alpha']
fp_beta = dic['fp_beta']
inside_FoV = dic['inside_FoV']


title = ur'$\alpha$ = %.1lf, $\beta$ = %.1lf, $\theta$ = %.1lf$^\circ$' % (pointing[0],pointing[1], orientation)
plt.suptitle(title, fontsize=fontsize)
plt.subplot(121)
plt.plot(alpha,beta,'.', markersize=2)
plt.plot(alpha[inside_FoV],beta[inside_FoV],'r.',markersize=2)  
plt.plot(fp_alpha,fp_beta,'k', linewidth=2)
plt.xlabel(ur'$\alpha$', fontsize=fontsize)
plt.ylabel(ur'$\beta$', fontsize=fontsize)
ax = plt.gca()
for tick in ax.xaxis.get_major_ticks():
  tick.label1.set_fontsize(tick_fontsize)
for tick in ax.yaxis.get_major_ticks():
  tick.label1.set_fontsize(tick_fontsize)
dalpha_sky = sky[1] - sky[0]
dbeta_sky = sky[3] - sky[2]
buffer_sky = 0.1
plt.xlim(sky[0] - buffer_sky*dalpha_sky, sky[1] + buffer_sky*dalpha_sky)
plt.ylim(sky[2] - buffer_sky*dbeta_sky, sky[3] + buffer_sky*dbeta_sky)

# Plot sources on focal plane
plt.subplot(122)
plt.plot(x,y,'r.', markersize=10)
plt.plot(fp_x, fp_y, 'k', linewidth=3)
plt.xlabel(ur'$x$', fontsize=fontsize)
plt.ylabel(ur'$y$', fontsize=fontsize)

fp_buffer = 0.1
dx = np.max(fp_x) - np.min(fp_x)
dy = np.max(fp_y) - np.min(fp_y)
plt.xlim(np.min(fp_x) - dx*fp_buffer, np.max(fp_x) + dx*fp_buffer)
plt.ylim(np.min(fp_y) - dy*fp_buffer, np.max(fp_y) + dy*fp_buffer)
plt.savefig("./camera_image.png",bbox_inches='tight',pad_inches=0.5)
plt.clf()
  
'''

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



  
 
'''