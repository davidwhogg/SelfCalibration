#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pylab as plt
import os
import pickle

# Custom Modules
import parameters
import god



#*****************************************************
#************ Transformation Functions ***************
#*****************************************************
def mag2flux(mag):
  return np.power(10,(-0.4*(mag-22.5)))
  
def flux2mag(flux):
  return 22.5-2.5*np.log10(flux)

def sky2fp(alpha, beta, pointing, orientation): 
  theta = - orientation*np.pi/180 # telescope rotates NOT sky
  x = (alpha - pointing[0]) * np.cos(theta) - (beta - pointing[1]) * np.sin(theta)
  y = (alpha - pointing[0]) * np.sin(theta) + (beta - pointing[1]) * np.cos(theta)
  return x,y 

def fp2sky(x, y, pointing, orientation):
  theta = - orientation*np.pi/180 # telescope rotates NOT sky
  alpha = x * np.cos(theta) + y * np.sin(theta) + pointing[0]
  beta = -x * np.sin(theta) + y * np.cos(theta) + pointing[1]
  return alpha, beta

#*****************************************************
#************* Single Image Functions ****************
#*****************************************************

class CameraCatalog:
    def __init__(self, sky_catalog, pointing,orientation): 
      self.k = sky_catalog.k.astype(int) 
      self.mag = sky_catalog.mag
      self.x, self.y = sky2fp(sky_catalog.alpha, sky_catalog.beta, pointing, orientation)
      self.size = sky_catalog.size
      self.flux = mag2flux(self.mag)
      self.epsilon = sky_catalog.epsilon

class MeasuredCatalog:
    def __init__(self, camera_catalog, inside_FoV): 
      self.size = len(inside_FoV[0])
      self.k = camera_catalog.k[inside_FoV] .astype(int)
      self.x = camera_catalog.x[inside_FoV]
      self.y = camera_catalog.y[inside_FoV]
      flat = god.flat_field(self.x,self.y)
      self.gods_invvar = self.true_invvar(camera_catalog, inside_FoV, flat)
      self.counts = camera_catalog.flux[inside_FoV] * flat + np.random.normal(size=self.size)/np.sqrt(self.gods_invvar)
      self.invvar = self.reported_invvar()

    def append(self, other):
       self.size = self.size + other.size
       self.k = np.append(self.k, other.k)
       self.x = np.append(self.x, other.x)
       self.y = np.append(self.y, other.y)
       self.counts = np.append(self.counts, other.counts)
       self.invvar = np.append(self.invvar, other.invvar)
       self.gods_invvar = np.append(self.gods_invvar, other.gods_invvar)

    def mag(self):
      return flux2mag(self.counts)
      
    def reported_invvar(self):
      sky_unc = 0.1*mag2flux(parameters.mag_at_ten_sigma())
      var = (sky_unc**2 + (parameters.eta()**2) * self.counts**2)
      return 1. / var

    def true_invvar(self, camera_catalog, inside_FoV, flat):
      epsilon = camera_catalog.epsilon[inside_FoV]
      sky_unc = 0.1*mag2flux(parameters.mag_at_ten_sigma())
      flux = camera_catalog.flux[inside_FoV]
      var = (sky_unc**2 * (1. + epsilon**2) + (parameters.eta()**2) * flat**2 * flux**2)
      return 1. / var

def single_image(sky_catalog, pointing, orientation, plots=None, verbose=None):
  if verbose != None: print "Converting sky catalog to focal plane coordinates..."
  camera_catalog = CameraCatalog(sky_catalog, pointing, orientation)
  if verbose != None: print "...done!"

  if verbose != None: print "Finding stars within camera FoV..."
  FoV = parameters.FoV()
  x_min = -FoV[0]/2; y_min = -FoV[1]/2
  x_max = FoV[0]/2; y_max = FoV[1]/2
  inside_FoV = np.where((x_min<camera_catalog.x) & (camera_catalog.x<x_max) & (y_min<camera_catalog.y) & (camera_catalog.y<y_max))
  if verbose != None: print "...done!"

  if verbose != None: print "Measuring stars within FoV..."
  measured_catalog = MeasuredCatalog(camera_catalog, inside_FoV)
  if verbose != None: print "...done!"
    
  if plots != None:
    # Plot sky (highlight extracted sources), with FoV
    plt.figure(2000, figsize=(13,6))
    title = ur'$\alpha$ = %.1lf, $\beta$ = %.1lf, $\theta$ = %.1lf$^\circ$' % (pointing[0],pointing[1], orientation)
    plt.suptitle(title, fontsize=20)
    plt.subplot(121)
    x = np.array([x_min, x_min, x_max, x_max, x_min])
    y = np.array([y_min, y_max, y_max, y_min, y_min])
    alpha, beta = fp2sky(x,y,pointing, orientation)
    plt.plot(sky_catalog.alpha,sky_catalog.beta,'o', markersize=2)
    plt.plot(sky_catalog.alpha[inside_FoV],sky_catalog.beta[inside_FoV],'ro',markersize=2)
    plt.plot(alpha,beta,'k', linewidth=2)
    plt.xlabel(ur'$\alpha$', fontsize=20)
    plt.ylabel(ur'$\beta$', fontsize=20)
    sky = parameters.sky_limits()
    plt.xlim(sky[0],sky[1])
    plt.ylim(sky[2],sky[3])
    plt.axis('equal')
    # Plot sources on focal plane
    plt.subplot(122)
    plt.plot(measured_catalog.x,measured_catalog.y,'o', markersize=2)
    plt.xlabel(ur'$x$', fontsize=20)
    plt.ylabel(ur'$y$', fontsize=20)
    plt.xlim(x_min,x_max)
    plt.ylim(y_min,y_max)
    filename = "Figures/Camera_Images/%s_alpha_%.1lf_beta_%.1lf_rot_%.1lf.png" % (plots, pointing[0],pointing[1], orientation)
    print filename
    plt.savefig(filename,bbox_inches='tight',pad_inches=0.5)
    plt.clf()

  return measured_catalog
  # measured_sources  *.size, *.k, *.flux, *.invvar, *.x, *.y

#*****************************************************
#**************** Survey Functions *******************
#*****************************************************
      
def survey(sky_catalog, survey_file, plots=None, verbose=None):  
  if verbose != None: print "Loading survey..."
  pointing = np.loadtxt(survey_file)
  number_pointings = len(pointing[:,0])  
  if verbose != None: print "...done!"
  
  if verbose != None: print "Surveying sky..."
  obs_cat = None
  for i in range(number_pointings):
    si = single_image(sky_catalog, [pointing[i,1],pointing[i,2]], pointing[i,3], plots=plots, verbose = verbose)
    if obs_cat is None:
      obs_cat = si
    else:
      obs_cat.append(si)

  if verbose != None: print "...done!"
  return obs_cat

#*****************************************************
#**************** Ubercal Functions ******************
#*****************************************************

def ubercalibration(observation_catalog,sky_catalog, strategy,ff_plots = None):
  order = parameters.flat_field_order
  q = np.array([1])
  stop_condition = 1e-4
  chi2 = 1e9
  old_chi2 = 1e10
  iteration_number = 0
  while (abs(chi2 - old_chi2) > stop_condition):
    temp_chi2 = chi2
    s, s_invvar = s_step(observation_catalog,q)
    q, q_invvar, chi2 = q_step(observation_catalog, s, order,iteration_number,plots=ff_plots)
    old_chi2 = temp_chi2
    
     # Calculate rms error in stars
    indx = [s != 0]
    rms = rms_error(s[indx],sky_catalog.flux[indx])
    bdness = badness(q)
    print "%i: RMS = %.6f %%, Badness = %0.6f %%, chi2 = %0.2f (%i)" % (iteration_number+1, rms, 100*bdness,chi2, observation_catalog.size)
    print q

    if (ff_plots == 'all') or (abs(chi2 - old_chi2) < stop_condition): plot_flat_fields(q, (iteration_number), strategy=strategy)
    
    if ff_plots == 'all' and (abs(chi2 - old_chi2) < stop_condition): 
      os.system(("convert -delay 20 -loop 0 ./Figures/Flat_Fields/%s*.png ./Figures/Flat_Fields/%s_00_animation.gif" % (strategy,strategy)))
    iteration_number += 1
  return 

def evaluate_flat_field_functions(x, y, order):
  L = (order + 1)*(order + 2)/2
  g = np.zeros((len(x),L))
  l = 0
  for n in range(order+1):
    for m in range(n+1):
      g[:,l] = (x**(n-m)) * (y**m)
      l += 1 
  return g

def evaluate_flat_field(x, y, q):  
  # Calculate required order
  order = int(np.around(np.sqrt(0.25+2*len(q))-1.5))
  assert(len(q) == ((order + 1) * (order + 2) / 2))
  g = evaluate_flat_field_functions(x, y, order)
  return np.dot(g, q)
  
def normalize_flat_field(q):
  fit_mean = average_over_ff(evaluate_flat_field, (q))
  return (q * god_mean_flat_field/fit_mean)

def diff_flat_field_squared(x, y, q):
  q = normalize_flat_field(q)
  return (evaluate_flat_field(x, y, q) - god.flat_field(x, y))**2

def badness(q):
  return np.sqrt(average_over_ff(diff_flat_field_squared, (q)))

def s_step(obs_cat, q):
  ff = evaluate_flat_field(obs_cat.x, obs_cat.y, q)
  fcss = ff * obs_cat.counts * obs_cat.invvar  
  ffss = ff * ff * obs_cat.invvar
  max_star = np.max(obs_cat.k)
  s = np.zeros(max_star+1)
  s_invvar = np.zeros(max_star)
  for ID in range(max_star):
    indx = (obs_cat.k == ID)
    denominator = np.sum(ffss[indx])
    s_invvar[ID] = denominator
    if denominator > 0.:
      s[ID] = np.sum(fcss[indx])/denominator
  return (s, s_invvar)
  
def q_step(obs_cat, s, order, iteration_number, plots=None):
  g = evaluate_flat_field_functions(obs_cat.x,obs_cat.y, order)
  ss = s[obs_cat.k]
  q_invvar = np.dot(np.transpose(g)*ss*ss*obs_cat.invvar,g)
  numerator = np.sum(np.transpose(g)*ss*obs_cat.counts*obs_cat.invvar,axis=1)
  q = np.dot(np.linalg.inv(q_invvar),numerator)
  q = normalize_flat_field(q)
  np.set_printoptions(precision=2) 
  chi2 = np.sum((np.dot(g,q) * ss - obs_cat.counts)**2 * obs_cat.invvar)
  return (q, q_invvar, chi2)
  
  
#*****************************************************
#*************** Analysis Functions ******************
#*****************************************************

def rms_error(flux_estimate, true_flux):
  return 100*np.sqrt(np.mean((flux_estimate-true_flux)/true_flux)**2)

def average_over_ff(func, args):
  FoV = parameters.FoV()
  nalpha, nbeta = parameters.ff_samples()
  dalpha = FoV[0]/(nalpha-1)
  dbeta = FoV[1]/(nbeta-1)
  x = np.arange(-FoV[0]/2+dalpha/2,FoV[0]/2,dalpha)
  y = np.arange(-FoV[1]/2+dbeta/2,FoV[1]/2,dbeta)
  X, Y = np.meshgrid(x, y)
  temp_x = np.reshape(X,-1)
  temp_y = np.reshape(Y,-1)
  temp_ff = func(temp_x, temp_y, args)
  ff = np.reshape(temp_ff, (len(X[:,0]),len(X[0,:])))
  return np.mean(ff)

#*****************************************************
#*************** Plotting Functions ******************
#*****************************************************

def plot_flat_fields(our_q, iteration_number,strategy=None):
  FoV = parameters.FoV() 
  
  plt.figure(3000,figsize=(13, 6)) # want square figure
  plt.subplot(121)
  plt.suptitle('Survey %s' % strategy, fontsize = 20)

  plt.title(r"Flat-Fields (God's = Black; Fitted = Red) Iteration: %i" % (iteration_number+1))
  plt.xlabel(r"$\alpha$")
  plt.ylabel(r"$\beta$")
  x = np.arange(-FoV[0]/2,FoV[0]/2,0.01)
  y = np.arange(-FoV[1]/2,FoV[1]/2,0.01)
  X, Y = np.meshgrid(x, y)
  
  # Have to reshape so that evaluate_flat_field() works 
  temp_x = np.reshape(X,-1)
  temp_y = np.reshape(Y,-1)
  temp_z = evaluate_flat_field(temp_x,temp_y,our_q)
  our_ff = np.reshape(temp_z, (len(X[:,0]),len(X[0,:])))
  god_q = god.flat_field_parameters()
  temp_z = god.flat_field(temp_x,temp_y)
  god_ff = np.reshape(temp_z, (len(X[:,0]),len(X[0,:])))

  # Find parameters for contour plot
  god_ff_max = np.max(god_ff)
  god_ff_min = np.min(god_ff)
  step = (god_ff_max-god_ff_min)/10
  
  CS = plt.contour(X,Y,god_ff,np.arange(god_ff_min,god_ff_max,step),colors='k')
  plt.clabel(CS, fontsize=9, inline=1)
  CS2 = plt.contour(X,Y,our_ff,np.arange(god_ff_min,god_ff_max,step),colors='r',alpha=0.5)

  FoV = parameters.FoV()
  plt.xlim(-FoV[0]/2, FoV[0]/2)
  plt.ylim(-FoV[1]/2, FoV[1]/2)
  
  # Write formulas on plot
  #our_formula = r"$f(x,y) =%+.2f%+.2fx%+.2fy%+.2fxx%+.2fxy%+.2fyy$" % (our_q[0], our_q[1], our_q[2], our_q[3], our_q[4], our_q[5])
  #plt.text(-.33,-.29,our_formula, color='r',bbox = dict(boxstyle="square",ec='w',fc='w', alpha=0.9), fontsize=9.5)
  #god_formula = r"$f(x,y) =%+.2f%+.2fx%+.2fy%+.2fxx%+.2fxy%+.2fyy$" % (god_q[0], god_q[1], god_q[2], god_q[3], god_q[4], god_q[5])
  #plt.text(-.33,-.33,god_formula, color='k',bbox = dict(boxstyle="square",ec='w',fc='w', alpha=0.9), fontsize=9.5)
  
  # Plot residual in flat-field
  plt.subplot(122)
  plt.title(r"Residual Error in Fitted Flat-Field (\%)")
  a = plt.imshow((100*(our_ff-god_ff)/god_ff),extent=(-FoV[0]/2,FoV[0]/2,-FoV[1]/2,FoV[1]/2), vmin = -1,vmax = 1, cmap='gray')
  plt.colorbar(a,shrink=0.7)
  plt.xlabel(r"$\alpha$")
  plt.ylabel(r"$\beta$")
  
  # Save figure
  if iteration_number < 9:
    filename = 'Figures/Flat_Fields/%s_0%d_ff.png' % (strategy,iteration_number+1)
  else:
    filename = 'Figures/Flat_Fields/%s_%d_ff.png' % (strategy, iteration_number+1)
  plt.savefig(filename,bbox_inches='tight',pad_inches=0.5)
  plt.clf()  

def coverage(obs_cat, strategy):
  limits = parameters.sky_limits()
  num_stars = int(np.around(parameters.density_of_stars()*(limits[1]-limits[0]) * (limits[3]-limits[2])))
  num_obs = np.zeros((num_stars,2))
  for i in range(num_stars):
    num_obs[i,0] = i
    indx = np.where(obs_cat.k == i)
    num_obs[i,1] = len(obs_cat.k[indx])
  max_obs = int(np.around(np.max(num_obs[:,1])))
  int_obs = np.zeros((max_obs+1,2))
  
  for i in range(max_obs+1):
    int_obs[i,0] = i
    indx = np.where(num_obs[:,1] <= i)
    int_obs[i,1] = len(indx[0])
  plt.figure()
  plt.plot(int_obs[:,0],int_obs[:,1]/float(num_stars))
  plt.xlabel(r"Number of Observations")
  plt.ylabel(r"Fraction of Sources Covered")
  plt.title('Survey %s Coverage' % strategy)
  plt.ylim(0.,1.)
  filename = "./Figures/%s_coverage.png" % strategy
  plt.savefig(filename,bbox_inches='tight')#,pad_inches=0.5)
  plt.clf()
  return num_obs

def invvar_saveout(obs_cat):
  dic = {}
  dic['true_invvar'] = obs_cat.gods_invvar
  dic['reported_invvar'] = obs_cat.invvar
  pickle.dump(dic, open( "./Plotting_Data/invvar.p", "wb" ) )
  plt.clf()
  plt.plot(flux2mag(obs_cat.counts), np.log10((1/obs_cat.gods_invvar)/obs_cat.counts**2),'rx')
  plt.plot(flux2mag(obs_cat.counts), np.log10((1/obs_cat.invvar)/obs_cat.counts**2),'k.')
  plt.show()
  return 0
# constants
god_mean_flat_field = average_over_ff(god.flat_field,god.flat_field_parameters())