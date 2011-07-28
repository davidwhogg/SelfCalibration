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
import functions as f
import copy
from multiprocessing import Pool

mult_proc = False
plot_ff = True
expct_perf = True

alpha = r'Sky Position $\alpha$ (deg$^2$)'
beta = r'Sky Position $\beta$ (deg$^2$)'

scale = 2
fig_width_pt = scale*415.55  # Get this from LaTeX using \showthe\columnwidth
inches_per_pt = 1.0/72.27               # Convert pt to inches
golden_mean = (np.sqrt(5)-1.0)/2.0         # Aesthetic ratio
fig_width = fig_width_pt*inches_per_pt  # width in inches
fig_height =fig_width*golden_mean       # height in inches
params = {'backend': 'pdf',
          'axes.labelsize': scale*10,
          'text.fontsize': scale*10,
          'legend.fontsize': scale*10,
          'xtick.labelsize': scale*9,
          'ytick.labelsize': scale*9,
          'linewidth' : scale*1,
          'text.usetex': True}
          #'font': {'family':'sans-serif','sans-serif':['Helvetica']}}
plt.rcParams.update(params)
plt.rc('font', family = 'serif', serif = 'cmr10')

if len(sys.argv) == 2:
  print "Running plotting routine..."
  out_dir = sys.argv[1]
  if out_dir[-1] == '/': out_dir = out_dir[:-1]
else:
  print "Error - no plotting directory given!"
  sys.exit()

# General Plotting Parameters
fontsize = 40
tick_fontsize = 24
double_fig_width = 28
single_fig_width = 12
single_fig_height = 12


def plot_flat_fields(map_dic):
  params = map_dic['params']
  out_dir = map_dic['out_dir']
  ff_filename = map_dic['ff_filename']
  strategy = map_dic['survey']
  print ff_filename
  FoV = params['FoV']
  plt.rcParams.update({'figure.figsize': [fig_width, fig_width*0.94], 'xtick.labelsize': scale*7, 'ytick.labelsize': scale*7})
  fig = plt.figure()
  # Axis for colorbar
  cax = fig.add_axes([0.93, 0.1, 0.03, 0.8])
  # Load Pickles
  dic = pickle.load(open(ff_filename))
  our_X = dic['x']
  our_Y = dic['y']
  our_ff = dic['our_ff']
  iteration_number = dic['iteration_number']
  god_dic = pickle.load(open(out_dir+'/god_ff.p'))
  god_ff = god_dic['god_ff']
  god_X = god_dic['x']
  god_Y = god_dic['y']
  bestfit_dic = pickle.load(open(out_dir+'/bestfit_ff.p'))
  bestfit_ff = bestfit_dic['bestfit_ff']
  bestfit_X = bestfit_dic['x']
  bestfit_Y = bestfit_dic['y']
  # Plot God and Fitted
  plt.subplot(221)
  # Find parameters for contour plot
  god_ff_max = np.max(god_ff)
  god_ff_min = np.min(god_ff)
  levels = np.arange(0.5,1.5,0.01)
  CS = plt.contour(god_X,god_Y,god_ff,levels ,colors='k',alpha=0.3)
  CS2 = plt.contour(our_X, our_Y, our_ff, levels,colors='k')
  plt.ylabel(beta)
  plt.clabel(CS2, inline=1)
  plt.xlim(-FoV[0]/2, FoV[0]/2)
  plt.ylim(-FoV[1]/2, FoV[1]/2)

  plt.gca().set_xticklabels([])
  # Plot residual in god and fitted
  plt.subplot(222)
  a = plt.imshow((100*(our_ff-god_ff)/god_ff),extent=(-FoV[0]/2,FoV[0]/2,-FoV[1]/2,FoV[1]/2), vmin = -0.5,vmax = 0.5, cmap='gray')
  plt.gca().set_xticklabels([])
  plt.gca().set_yticklabels([])
  # Plot best-in-basis and fitted
  plt.subplot(223)
  plt.xlabel(alpha)
  plt.ylabel(beta)
  levels = np.arange(0.5,1.5,0.01)
  CS = plt.contour(bestfit_X,bestfit_Y,bestfit_ff,levels ,colors='k',alpha=0.3)
  CS2 = plt.contour(our_X, our_Y, our_ff, levels,colors='k')
  plt.clabel(CS2, inline=1)
  plt.xlim(-FoV[0]/2, FoV[0]/2)
  plt.ylim(-FoV[1]/2, FoV[1]/2)

  # Plot residual between fitted and best-in-basis
  plt.subplot(224)
  a = plt.imshow((100*(our_ff-bestfit_ff) /bestfit_ff),extent=(-FoV[0]/2,FoV[0]/2,-FoV[1]/2,FoV[1]/2), vmin = -0.5,vmax = 0.5, cmap='gray')
  plt.xlabel(alpha) 

  plt.gca().set_yticklabels([])
  # Squeeze subplots together
  plt.subplots_adjust(wspace=0.0,hspace=0.0)
  # Add Colorbar
  fig.colorbar(a, cax, orientation = 'vertical')
  filename = string.replace(ff_filename, '.p', '.png')  
  plt.savefig(filename,bbox_inches='tight')
  plt.clf() 
  
def camera_image(params, out_dir, camera_filename):
  print camera_filename
  plt.figure(figsize = (fig_width, 0.5*fig_width))
  plt.clf()
  
  dic = pickle.load(open(camera_filename))
  x = dic['measured_catalog.x'] 
  y = dic['measured_catalog.y']
  Alpha = dic['sky_catalog.alpha'] 
  Beta = dic['sky_catalog.beta']
  sky = dic['sky'] 
  pointing = dic['pointing']
  orientation = dic['orientation']
  fp_x = dic['fp_x']
  fp_y = dic['fp_y']
  fp_alpha = dic['fp_alpha']
  fp_beta = dic['fp_beta']
  inside_FoV = dic['inside_FoV']
  plt.subplot(121)
  plt.plot(Alpha,Beta,'k.', alpha = 0.2, markersize=3)
  #plt.plot(alpha[inside_FoV],beta[inside_FoV],'k.',markersize=5)  
  plt.plot(fp_alpha,fp_beta,'k', linewidth=2)
  plt.xlabel(alpha)
  plt.ylabel(beta)
  dalpha_sky = sky[1] - sky[0]
  dbeta_sky = sky[3] - sky[2]
  buffer_sky = 0.1
  plt.axis('scaled')
  plt.xlim(sky[0] - buffer_sky*dalpha_sky, sky[1] + buffer_sky*dalpha_sky)
  plt.ylim(sky[2] - buffer_sky*dbeta_sky, sky[3] + buffer_sky*dbeta_sky)
  

  # Plot sources on focal plane
  plt.subplot(122)
  plt.plot(x,y,'k.', markersize=6)
  plt.plot(fp_x, fp_y, 'k', linewidth=3)
  plt.xlabel(ur'Focal Plane Position $x$ (deg$^2$)')
  plt.ylabel(ur'Focal Plane Position $y$ (deg$^2$)')
  fp_buffer = 0.1
  dx = np.max(fp_x) - np.min(fp_x)
  dy = np.max(fp_y) - np.min(fp_y)
  plt.axis('scaled')
  plt.xlim(np.min(fp_x) - dx*fp_buffer, np.max(fp_x) + dx*fp_buffer)
  plt.ylim(np.min(fp_y) - dy*fp_buffer, np.max(fp_y) + dy*fp_buffer)
  god_dic = pickle.load(open(out_dir+'/god_ff.p'))
  god_ff = god_dic['god_ff']
  god_X = god_dic['x']
  god_Y = god_dic['y'] 
  levels = np.arange(0.5,1.5,0.01)
  CS = plt.contour(god_X,god_Y,god_ff,levels ,colors='k',alpha=0.3)
  #plt.clabel(CS, inline=1,color = '0.3')
  
  filename = string.replace(camera_filename, '.p', '.pdf')
  plt.savefig(filename,bbox_inches='tight')
  plt.clf()

def SurveyInnvar(params, survey_filemane, invvar_filename):
  print survey_filename
  plt.figure(figsize = (fig_width, 0.5*fig_width))
  plt.clf()
  
  plt.subplot(121)
  survey_dic = pickle.load(open(survey_filename))
  mag = survey_dic['mag']
  all_mag = survey_dic['all_sources']
  fit_mag = survey_dic['fit_mag']
  fit_dens = survey_dic['fit_dens']
  plt.semilogy(fit_mag, fit_dens, 'k:')
  hist_mag = np.arange(params['m_min'], params['m_max']+0.5, 0.5)
  plt.bar(hist_mag[:-1], np.histogram(all_mag, bins = hist_mag)[0], width = 0.5, color = 'w')
  plt.bar(hist_mag[:-1], np.histogram(mag, bins = hist_mag)[0], width = 0.5, hatch = '/', color = 'k', alpha =0.3)
  plt.xlabel(r'Source Magnitude (AB)')
  plt.ylabel(r'Sources within Survey ((0.5 mag)$^{-1})$')
  
  print invvar_filename
  plt.subplot(122)
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
  plt.plot(np.log10(counts), np.log10((1/true_invvar)/counts**2),'k.', markersize = 2., label = "Actual Variance", alpha = 0.1)
  plt.plot(np.log10(sort_counts), np.log10((1/sort_reported_invvar)/sort_counts**2),'k', label = "Assumed Variance")
  plt.xlabel(r'$\log_{10}(c_i)$')
  plt.ylabel(ur'$\log_{10}(\frac{{\sigma_i}^2}{c_i^2})$')
  plt.xlim(np.min(np.log10(counts)), np.max(np.log10(counts)))
  plt.subplots_adjust(wspace=0.3)
  filename = string.replace(invvar_filename, '.p', '.pdf')
  plt.savefig(filename,bbox_inches='tight')
  plt.clf()

def plot_performance(sln, mod_param, solution_path):
  print "Plotting Solution: ", solution_path
  sln = sln[sln[:,0].argsort(),:]
  if mod_param == 'density_of_stars':
    mod_param = r'Density of Sources (deg$^{-2}$)'
  elif mod_param == 'epsilon_max':
    mod_param = r'$\epsilon_{max}$'
  else:
    mod_param = string.replace(mod_param, '_', ' ')
  plt.figure(figsize = (1.1*fig_width, 0.5*fig_width))
  plt.subplot(121)
  plt.loglog(sln[:,0], sln[:,3], 'kx',label = r'True')
  plt.loglog(sln[:,0], sln[:,4], 'k.', label = r'Best-in-Basis')
  if expct_perf: plt.loglog(sln[:,0], (0.04*(sln[:,0])**-0.5), 'k:')
  #plt.legend(loc = 'upper right')
  plt.ylabel(r"Badness (\%)")
  plt.xlabel(mod_param)
  plt.subplot(122)
  plt.loglog(sln[:,0], sln[:,2], 'k.')
  plt.ylabel(r"Source Error (\%)")
  plt.xlabel(mod_param)
  plt.subplots_adjust(wspace=0.3)
  plt.savefig(solution_path + '/performance.pdf',bbox_inches='tight')
  plt.clf()


def plot_solution(sln, mod_param, solution_path):
  sln = sln[sln[:,0].argsort(),:]
  mod_param = string.replace(mod_param, '_', ' ')
  print solution_path
  plt.figure()
  plt.subplot(221)
  plt.plot(sln[:,0], sln[:,1])
  plt.ylabel("Number of Iterations to Converge")
  plt.xlabel(mod_param)
  plt.subplot(222)
  plt.plot(sln[:,0], sln[:,2])
  plt.ylabel(r"RMS Source Error (\%)")
  plt.xlabel(mod_param)
  plt.subplot(223)
  plt.plot(sln[:,0], sln[:,3])
  plt.plot(sln[:,0], sln[:,4])
  plt.ylabel(r"Badness of Flat-Field (\%)")
  plt.xlabel(mod_param)
  plt.subplot(224)
  plt.plot(sln[:,0], sln[:,5])
  plt.ylabel("$\chi^2$")
  plt.xlabel(mod_param)
  ax = plt.gca()
  ax.yaxis.major.formatter.set_powerlimits((0,0))
  plt.savefig(solution_path + '/solution.pdf',bbox_inches='tight',pad_inches=0.1)
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
    sln = np.array([])
    for isfile in temp:
      if os.path.isfile(out_dir+ '/' + srvy + '/' + isfile):
	modified_value_dir.remove(isfile)
    for mod_val_dir in modified_value_dir:
      dir_path = out_dir + '/' + srvy + '/' + mod_val_dir
      # Import Simulation Parameters
      params = pickle.load(open(('./%s/parameters.p' % dir_path)))
      # Plot Flat-Fields
      if plot_ff:
        ff_path = dir_path + '/FF'
        os.system("rm -r %s/*.png" % ff_path)
        os.system("rm -r %s/*.pdf" % ff_path)
        FFs = os.listdir(ff_path)
        map_dictionaries = []
        for ff in FFs:
          map_dic = {}
          map_dic['params'] = params
          map_dic['out_dir'] = out_dir
          map_dic['ff_filename'] = (ff_path + '/' + ff)
          map_dic['survey'] = srvy
          map_dictionaries.append(copy.deepcopy(map_dic))
        if mult_proc:
          p = Pool(8)
          p.map(plot_flat_fields,map_dictionaries)      
        else: map(plot_flat_fields,map_dictionaries)
        # Create Animations
        png_dir = ff_path
        command = ('convert -delay 50 -loop 0 %s/*.png %s/ff_animation.gif' % (png_dir, dir_path))
        os.system(command)
      # Plot Camera Image
      if os.path.isfile((dir_path + '/camera_image.p')): camera_image(params, out_dir, (dir_path + '/camera_image.p'))
      # Plot Inverse Invariance 
      survey_filename = dir_path+'/source_catalog.p'
      invvar_filename = dir_path + '/invvar.p'
      if os.path.isfile(survey_filename) and os.path.isfile(invvar_filename): SurveyInnvar(params, dir_path, invvar_filename)
    
    # plot iteration number, badness, rms, chi2 
      solution_path = dir_path + '/solution.p'
      if os.path.isfile(solution_path):
        sln_dic = pickle.load(open(solution_path))
        temp_sln = np.array([[sln_dic['mod_value'], sln_dic['iter_no'], sln_dic['rms'], sln_dic['badness'], sln_dic['badness_bestfit'], sln_dic['chi2']]])
        if len(sln) == 0: sln = temp_sln
        else: sln = np.append(sln, temp_sln, axis=0)
    if len(sln.shape) > 1:
     if len(sln[:,0]) > 1:
       plot_solution(sln, sln_dic['mod_param'], (out_dir + '/' + srvy))
       plot_performance(sln, sln_dic['mod_param'], (out_dir + '/' + srvy))


