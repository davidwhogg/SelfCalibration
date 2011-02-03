#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pylab as plt
import pickle
import os
import string
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']))
rc('text', usetex=True)


# Data directory
dir_path = 'temp'

# General Plotting Parameters
fontsize = 20
tick_fontsize = 12
double_fig_width = 14
single_fig_width = 6
single_fig_height = 6

# Import general parameters
simulation_parameters = pickle.load(open(('./%s/simulation_parameters.p' % dir_path)))
sky_limits = simulation_parameters['sky_limits']
m_min = simulation_parameters['m_min']
m_max = simulation_parameters['m_max']
survey_strategies = simulation_parameters['survey_strategies']
FoV = simulation_parameters['FoV'] 

def plot_sky_catalog(sky_catalog_filename):
  pickle_dic = pickle.load(open((sky_catalog_filename)))
  alpha = pickle_dic['alpha'] 
  beta = pickle_dic['beta'] 
  mag = pickle_dic['mag']

  plt.figure(figsize=(double_fig_width, single_fig_height))
  st = plt.suptitle("Sky Catalog", fontsize=fontsize)
  # Plot portion of sky
  plt.subplot(121)
  plt.plot(alpha,beta,'.k',markersize=1)
  plt.xlabel(ur'$\alpha$', fontsize=fontsize)
  plt.ylabel(ur'$\beta$', fontsize=fontsize)
  plt.xlim(sky_limits[0] - (sky_limits[1]-sky_limits[0])*0.1, sky_limits[1] + (sky_limits[1]-sky_limits[0])*0.1)
  plt.ylim(sky_limits[2] - (sky_limits[3]-sky_limits[2])*0.1, sky_limits[3] + (sky_limits[3]-sky_limits[2])*0.1)
  ax = plt.gca()
  for tick in ax.xaxis.get_major_ticks():
    tick.label1.set_fontsize(fontsize)
  for tick in ax.yaxis.get_major_ticks():
    tick.label1.set_fontsize(fontsize)
  # Histogram of source magnitude
  plt.subplot(122)
  bins=np.arange(m_min,m_max,0.05)
  plt.hist(mag,bins=bins, log=True)
  plt.xlabel("Source Magnitude", fontsize=fontsize)
  plt.ylabel("log(N)", fontsize=fontsize)
  ax = plt.gca()
  for tick in ax.xaxis.get_major_ticks():
    tick.label1.set_fontsize(fontsize)
  for tick in ax.yaxis.get_major_ticks():
    tick.label1.set_fontsize(fontsize)
    
  filename = dir_path + '/sky_catalog.png'
  plt.savefig(filename,bbox_inches='tight',pad_inches=0.)
  plt.clf()

def plot_invvar(invvar_filename):
  plt.figure(figsize=(single_fig_width, single_fig_height))
  
  pickle_dic = pickle.load(open(invvar_filename))
  counts = pickle_dic['counts']
  true_invvar = pickle_dic['true_invvar']
  reported_invvar = pickle_dic['reported_invvar']  
  # Sort true invvar by counts so we can plot as a line
  sort = np.zeros((len(counts), 2))
  sort[:,0] = counts[:]
  sort[:,1] = reported_invvar[:]
  sort = sort[sort[:,0].argsort(),:]
  sort_reported_invvar = sort[:,1]
  sort_counts = sort[:,0]
  plt.plot(counts, np.log10((1/true_invvar)/counts**2),'r.', markersize = 1., label = "Actual Variance")
  plt.plot(sort_counts, np.log10((1/sort_reported_invvar)/sort_counts**2),'k', label = "Assumed Variance")
  plt.xlabel(r'$c_i$', fontsize=fontsize)
  plt.ylabel(ur'$\log_{10}(\frac{{\sigma_i}^2}{c_i^2})$',fontsize=fontsize)#, rotation='horizontal')
  plt.legend()
  ax = plt.gca()
  for tick in ax.xaxis.get_major_ticks():
    tick.label1.set_fontsize(tick_fontsize)
  for tick in ax.yaxis.get_major_ticks():
    tick.label1.set_fontsize(tick_fontsize)


  filename = dir_path + '/invvar.png'
  plt.savefig(filename,bbox_inches='tight',pad_inches=0.)
  plt.clf()

def plot_coverage():  
  plt.figure(figsize=(single_fig_width, single_fig_height))
  for j in range(len(survey_strategies)):
    coverage_filename =('%s/%s_coverage.p' % (dir_path, survey_strategies[j])) 
    if os.path.exists(coverage_filename):
      print ("Plotting sky coverage for Strategy %s..." % survey_strategies[j])
      pickle_dic = pickle.load(open((coverage_filename)))
      number_stars = pickle_dic['number_stars']
      k = pickle_dic['k']
      strategy = pickle_dic['strategy']
      num_obs_of_star = np.bincount(k)
      num_star_with_N_obs = np.bincount(num_obs_of_star)
      # "Integrate" 
      for i in range(1,len(num_star_with_N_obs)):
        num_star_with_N_obs[i] += num_star_with_N_obs[i-1]
      fraction = (1.*num_star_with_N_obs)/number_stars
      plt.plot(fraction, label = ('Strategy '+survey_strategies[j]))  
      print "...done!"
  
    else: print ('... no file for Strategy %s!' % survey_strategies[j])
  
  plt.xlabel(r"Number of Observations")
  plt.ylabel(r"Fraction of Sources Covered")
  plt.title('Survey Coverage')
  plt.legend(loc=4)
  filename = "%s/coverage.png" % dir_path
  ax = plt.gca()
  for tick in ax.xaxis.get_major_ticks():
    tick.label1.set_fontsize(tick_fontsize)
  for tick in ax.yaxis.get_major_ticks():
    tick.label1.set_fontsize(tick_fontsize)
  plt.ylim(0.,1.1)
  filename = dir_path + '/coverage.png'
  plt.savefig(filename,bbox_inches='tight',pad_inches=0.)
  plt.clf()

def camera_image(camera_filename):
  plt.figure(figsize=(double_fig_width, single_fig_height))
  dic = pickle.load(open(camera_filename))
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
  filename = dir_path + '/camera.png'
  plt.savefig(filename,bbox_inches='tight',pad_inches=0.5)
  plt.clf()

def plot_flat_fields(ff_filename, strategy):
  plt.figure(figsize=(double_fig_width, single_fig_height))
  dic = pickle.load(open(ff_filename))
  X = dic['x']
  Y = dic['y']
  our_ff = dic['our_ff']
  god_ff = dic['god_ff']
  iteration_number = dic['iteration_number']
  
  plt.subplot(121)
  plt.suptitle('Survey %s' % strategy, fontsize = fontsize)
  plt.title(r"Flat-Fields (God's = Black; Fitted = Red) Iteration: %i" % (iteration_number))
  plt.xlabel(r"$\alpha$")
  plt.ylabel(r"$\beta$")

  # Find parameters for contour plot
  god_ff_max = np.max(god_ff)
  god_ff_min = np.min(god_ff)
  # XX magic number
  levels = np.arange(0.5,1.5,0.01)
  CS = plt.contour(X,Y,god_ff,levels ,colors='k')
  plt.clabel(CS, fontsize=9, inline=1)
  CS2 = plt.contour(X,Y,our_ff,levels,colors='r',alpha=0.5)

  plt.xlim(-FoV[0]/2, FoV[0]/2)
  plt.ylim(-FoV[1]/2, FoV[1]/2)
  
  # Plot residual in flat-field
  plt.subplot(122)
  plt.title(r"Residual (Fit - God) in Flat-Field (\%)")
  a = plt.imshow((100*(our_ff-god_ff)/god_ff),extent=(-FoV[0]/2,FoV[0]/2,-FoV[1]/2,FoV[1]/2), vmin = -1,vmax = 1, cmap='gray')
  plt.colorbar(a,shrink=0.7)
  plt.xlabel(r"$\alpha$")
  plt.ylabel(r"$\beta$")
  
  filename = string.replace(ff_filename, '.p', '.png')  
  plt.savefig(filename,bbox_inches='tight',pad_inches=0.5)
  plt.clf()  

if __name__ == "__main__":
  os.system(('rm %s/*.png' % dir_path))
  os.system(('rm %s/*.gif' % dir_path))
  for ii in range(len(survey_strategies)):
    os.system(('rm %s/%s/*.png' % (dir_path, survey_strategies[ii])))
  # Plot survey catalog
  sky_catalog_filename = dir_path + '/source_catalog.p'
  if os.path.exists(sky_catalog_filename) == True:
    print "Plotting sky catalog..."
    plot_sky_catalog(sky_catalog_filename)
    print "...done!"
  else: print "No file for sky plots..."
  
  # Plot invvar
  invvar_filename = dir_path + '/invvar.p'
  if os.path.exists(invvar_filename) == True:
    print "Plotting invvar..."
    plot_invvar(invvar_filename)
    print "...done!"
  else: print "No file for invvar plots..."

  # Plot camera image
  camera_filename = dir_path + '/camera_image.p'
  if os.path.exists(camera_filename) == True:
    print "Plotting camera image..."
    camera_image(camera_filename)
    print "...done!"
  else: print "No file for camera image..."
  
  # Plot Flat Fields
  for ii in range(len(survey_strategies)):
    print ("Plotting Survey %s flat fields..." % survey_strategies[ii])
    number_ff = os.listdir(('%s/%s/' % (dir_path, survey_strategies[ii])))
    for jj in range(len(number_ff)):
      ff_filename = '%s/%s/%s' % (dir_path, survey_strategies[ii], number_ff[jj])
      plot_flat_fields(ff_filename, survey_strategies[ii])
    
    # Create Animations
    png_dir = ('%s/%s/' % (dir_path, survey_strategies[ii]))
    out_dir = ('%s/' % (dir_path))
    print "...animating..."
    command = ('convert -delay 20 -loop 0 %s*.png %s%s_animation.gif' % (png_dir, out_dir, 
    survey_strategies[ii]))
    os.system(command)
    print "...done!"
  
  # Plotting Coverage
  plot_coverage()