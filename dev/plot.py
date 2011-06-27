#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Rory Holmes 2011
# Script to plot all the output generated in the cross-calibration simulation
# Must call with the output directory from cross-cal, e.g: "./plot.py dir"

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pylab as plt
import pickle
import os
import sys
import string
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']))
rc('text', usetex=True)

if len(sys.argv) == 2:
  print "Running plotting routine..."
  out_dir = sys.argv[1]
else:
  print "Error - no plotting directory given!"
  sys.exit()

# General Plotting Parameters
fontsize = 20
tick_fontsize = 12
double_fig_width = 14
single_fig_width = 6
single_fig_height = 6


def plot_flat_fields(params, ff_filename, strategy):
  print ff_filename
  FoV = params['FoV']
  plt.figure(figsize=(double_fig_width, single_fig_height))
  dic = pickle.load(open(ff_filename))
  X = dic['x']
  Y = dic['y']
  our_ff = dic['our_ff']
  god_ff = dic['god_ff']
  iteration_number = dic['iteration_number']
  plt.subplot(121)
  plt.suptitle('Survey %s' % strategy, fontsize = fontsize)
  plt.title(r"Flat-Fields (True = Black; Fitted = Red) Iteration: %i" % (iteration_number))
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
  plt.title(r"Residual (Fit - True) in Flat-Field (\%)")
  a = plt.imshow((100*(our_ff-god_ff)/god_ff),extent=(-FoV[0]/2,FoV[0]/2,-FoV[1]/2,FoV[1]/2), vmin = -0.5,vmax = 0.5, cmap='gray')
  plt.colorbar(a,shrink=0.7)
  plt.xlabel(r"$\alpha$")
  plt.ylabel(r"$\beta$")
  
  filename = string.replace(ff_filename, '.p', '.png')  
  plt.savefig(filename,bbox_inches='tight',pad_inches=0.5)
  plt.clf() 
  
def camera_image(params, camera_filename):
  print camera_filename
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
  filename = string.replace(camera_filename, '.p', '.png')
  plt.savefig(filename,bbox_inches='tight',pad_inches=0.5)
  plt.clf()

def plot_invvar(params, invvar_filename):
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
  filename = string.replace(invvar_filename, '.p', '.png')
  plt.savefig(filename,bbox_inches='tight',pad_inches=0.1)
  plt.clf()

def plot_solution(solution_path):
  print solution_path
  plt.figure()
  pickle_dic = pickle.load(open(solution_path)) # [modified_param, iteration number, rms, badness, chi2]
  param_range = pickle_dic['parameter_range']
  mod_param = pickle_dic['modified_parameter']
  mod_param = string.replace(mod_param, '_', ' ') 
  rms = pickle_dic['rms']
  bdnss = pickle_dic['bdnss']
  it_num = pickle_dic['it_num']
  chi2 = pickle_dic['chi2']
  plt.subplot(221)
  plt.plot(param_range, it_num)
  plt.ylabel("Number of Iterations to Converge")
  plt.xlabel(mod_param)
  plt.subplot(222)
  plt.plot(param_range, rms)
  plt.ylabel(u"RMS Source Error (\%)")
  plt.xlabel(mod_param)
  plt.subplot(223)
  plt.plot(param_range, bdnss)
  plt.ylabel(u"Badness of Flat-Field (\%)")
  plt.xlabel(mod_param)
  plt.subplot(224)
  plt.plot(param_range, chi2)
  plt.ylabel("$\chi^2$")
  plt.xlabel(mod_param)
  ax = plt.gca()
  ax.yaxis.major.formatter.set_powerlimits((0,0))
  filename = string.replace(solution_path, '.p', '.png')
  plt.savefig(filename,bbox_inches='tight',pad_inches=0.1)
  plt.clf()

if __name__ == "__main__":
  # Find the surveys done
  surveys = os.listdir(out_dir)
  temp = 1*surveys
  for isfile in temp:
    if os.path.isfile(out_dir+ '/' + isfile):
      surveys.remove(isfile)
  for srvy in surveys:
    # Find Modified Parameter Directories
    modified_value_dir = os.listdir(out_dir + '/' + srvy)
    temp = 1*modified_value_dir
    for isfile in temp:
      if os.path.isfile(out_dir+ '/' + srvy + '/' + isfile):
	modified_value_dir.remove(isfile)
    for mod_val_dir in modified_value_dir:
      dir_path = out_dir + '/' + srvy + '/' + mod_val_dir
      # Import Simulation Parameters
      params = pickle.load(open(('./%s/parameters.p' % dir_path)))
      # Plot Flat-Fields
      ff_path = dir_path + '/FF'
      os.system("rm -r %s/*.png" % ff_path)
      FFs = os.listdir(ff_path)
      for ff in FFs:
	plot_flat_fields(params, (ff_path + '/' + ff), srvy)
      # Create Animations
      png_dir = ff_path
      command = ('convert -delay 50 -loop 0 %s/*.png %s/ff_animation.gif' % (png_dir, dir_path))
      os.system(command)
      # Plot Camera Image
      camera_image(params, (dir_path + '/camera_image.p'))
      # Plot Inverse Invariance 
      plot_invvar(params, (dir_path + '/invvar.p'))
    
    # plot iteration number, badness, rms, chi2 
    solution_path = out_dir + '/' + srvy + '/solution.p'
    if os.path.isfile(solution_path): plot_solution(solution_path)

# Import general parameters
'''
simulation_parameters = pickle.load(open(('./%ssimulation_parameters.p' % dir_path)))
sky_limits = simulation_parameters['sky_limits']
m_min = simulation_parameters['m_min']
m_max = simulation_parameters['m_max']
survey_strategies = simulation_parameters['survey_strategies']
FoV = simulation_parameters['FoV'] 

  if plot_bdness:
    # Plot badness against parameter
    select_files = os.listdir(dir_path)
    temp = 1*select_files
    for isfile in temp:
      if os.path.isfile(dir_path + isfile) != True:
        select_files.remove(isfile)
    select_files.remove('simulation_parameters.p')
    for data_files in select_files:
      data = np.loadtxt((dir_path+data_files))
      xdata = data[:,0]
      badness = data[:,1]
      plt.plot((xdata), (badness), label = data_files[0])
      plt.xlabel((bdness_plot_xlabel))
      plt.ylabel("Badness")
    plt.legend()
    filename = dir_path+'Badness.png'
    plt.savefig(filename,bbox_inches='tight',pad_inches=0.5)
    plt.clf() 
    
  if plot_rms:
    # Plot badness against parameter
    select_files = os.listdir(dir_path)
    temp = 1*select_files
    for isfile in temp:
      if os.path.isfile(dir_path + isfile) != True:
        select_files.remove(isfile)
    select_files.remove('simulation_parameters.p')
    select_files.remove('Badness.png')
    for data_files in select_files:
      data = np.loadtxt((dir_path+data_files))
      xdata = data[:,0]
      rms = data[:,2]
      plt.plot((xdata), (rms), label = data_files[0])
      plt.xlabel((bdness_plot_xlabel))
      plt.ylabel("rms")      
    plt.legend()
    filename = dir_path+'/rms.png'
    plt.savefig(filename,bbox_inches='tight',pad_inches=0.5)
    plt.clf() 
    
'''