#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import os
import pickle
import string
import scipy.optimize
# Custom Modules
import god

from master import init_func
pdic, directory_path, temp1, temp2 = init_func() # import parameter database from main module

#*****************************************************
#************ Transformation Functions ***************
#*****************************************************
def mag2flux(mag):
  return np.power(10,(-0.4*(mag-22.5)))
  
def flux2mag(flux):
  return 22.5-2.5*np.log10(flux)

def sky2fp(alpha, beta, pointing, orientation): 
  theta = - orientation*np.pi/180. # telescope rotates NOT sky
  x = (alpha - pointing[0]) * np.cos(theta) - (beta - pointing[1]) * np.sin(theta)
  y = (alpha - pointing[0]) * np.sin(theta) + (beta - pointing[1]) * np.cos(theta)
  return x,y 

def fp2sky(x, y, pointing, orientation):
  theta = - orientation*np.pi/180. # telescope rotates NOT sky
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
      sky_unc = 0.1*mag2flux(pdic['mag_at_ten_sigma'])
      var = (sky_unc**2 + (pdic['eta']**2) * self.counts**2)
      return 1. / var

    def true_invvar(self, camera_catalog, inside_FoV, flat):
      epsilon = camera_catalog.epsilon[inside_FoV]
      sky_unc = 0.1*mag2flux(pdic['mag_at_ten_sigma'])
      flux = camera_catalog.flux[inside_FoV]
      var = (sky_unc**2 * (1. + epsilon**2) + (pdic['eta']**2) * flat**2 * flux**2)
      return 1. / var

def single_image(sky_catalog, pointing, orientation, plots=None, verbose=None):
  if verbose: print "Converting sky catalog to focal plane coordinates..."
  camera_catalog = CameraCatalog(sky_catalog, pointing, orientation)
  if verbose: print "...done!"

  if verbose: print "Finding stars within camera FoV..."
  FoV = pdic['FoV']
  x_min = -0.5*FoV[0]; y_min = -0.5*FoV[1]
  x_max = 0.5*FoV[0]; y_max = 0.5*FoV[1]
  inside_FoV = np.where((x_min<camera_catalog.x) & (camera_catalog.x<x_max) & (y_min<camera_catalog.y) & (camera_catalog.y<y_max))
  if verbose: print "...done!"

  if verbose: print "Measuring stars within FoV..."
  measured_catalog = MeasuredCatalog(camera_catalog, inside_FoV)
  if verbose: print "...done!"
  temp_filename = directory_path + '/camera_image.p'
  one_camera_file = os.path.exists(temp_filename)
  if plots and (one_camera_file != True) and (len(inside_FoV[0]) >= 5): save_camera(sky_catalog, measured_catalog, inside_FoV, pointing, orientation, verbose = verbose)
  return measured_catalog
  # measured_sources  *.size, *.k, *.flux, *.invvar, *.x, *.y

#*****************************************************
#**************** Survey Functions *******************
#*****************************************************
      
def survey(sky_catalog, survey_file, plots=None, verbose=None):  
  if verbose: print "Loading survey..."
  pointing = np.loadtxt(survey_file)
  number_pointings = len(pointing[:,0])  
  if verbose: print "...done!"
  
  if verbose: print "Surveying sky..."
  obs_cat = None
  for i in range(number_pointings):
    si = single_image(sky_catalog, [pointing[i,1],pointing[i,2]], pointing[i,3], plots=plots, verbose = verbose)
    if obs_cat is None:
      obs_cat = si
    else:
      obs_cat.append(si)

  if verbose: print "...done!"
  return obs_cat

#*****************************************************
#**************** Ubercal Functions ******************
#*****************************************************

def ubercalibration(observation_catalog,sky_catalog, strategy,modified_parameter, modified_value, ff_plots = None):
  order = pdic['flat_field_order']
  q = np.array([1])
  stop_condition = 1e-6
  chi2 = 1e9
  old_chi2 = 1e10
  iteration_number = 0
  next_plot_iteration = 1
  while ((abs(chi2 - old_chi2) > stop_condition) and (iteration_number < 258)):
    iteration_number += 1
    temp_chi2 = chi2
    s, s_invvar = s_step(observation_catalog,q)
    q, q_invvar, chi2 = q_step(observation_catalog, s, order,iteration_number,plots=ff_plots)
    old_chi2 = temp_chi2
    
     # Calculate rms error in stars
    indx = [s != 0]
    rms = rms_error(s[indx],sky_catalog.flux[indx])
    bdness = badness(q)
    print "%i: RMS = %.6f %%, Badness = %0.6f %%, chi2 = %0.2f (%i)" % (iteration_number, rms, bdness,chi2, observation_catalog.size)
    print q

    if (ff_plots and (iteration_number == next_plot_iteration)) or (abs(chi2 - old_chi2) < stop_condition): 
      saveout_flat_fields(q, iteration_number, strategy=strategy)
      next_plot_iteration *= 2
    
    if (ff_plots and (abs(chi2 - old_chi2) < stop_condition)):
      if modified_value != None:
        modified_value = float(modified_value)
        out = np.zeros((1,4))
        out[0,0] = modified_value
        out[0,1] = bdness
        out[0,2] = rms
        out[0,3] = chi2
        filename = '%s/%s_%s.txt' % (string.rstrip(directory_path, ('/'+str(modified_value))),  strategy, modified_parameter)
        #print filename
        f_handle = file(filename, 'a')
        np.savetxt(f_handle, out)
        f_handle.close()
    
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

def badness_old(s, s_true):
  f0 = 1.0
  df = 0.01
  x = np.array([-1.,0.,1.])
  f = f0 + df * x
  b = [np.mean(((f1*s - s_true)/s_true)**2) for f1 in f]
  # second-order polynomial fit a0 + a1 x + a2 xx
  a = [b[1], 0.5*(b[2]-b[0]), 0.5*(b[2]-2.*b[1]+b[0])]
  a = np.polyfit(f,b,2)
  xmin = - 0.5 * a[1] / a[2]
  return np.sqrt(a[0] + a[1] * xmin + a[2] * xmin * xmin)

def badness(q):
  FoV = pdic['FoV']
  nalpha, nbeta = pdic['ff_samples']
  dalpha = FoV[0]/(nalpha-1)
  dbeta = FoV[1]/(nbeta-1)
  x = np.arange(-FoV[0]/2+dalpha/2,FoV[0]/2,dalpha)
  y = np.arange(-FoV[1]/2+dbeta/2,FoV[1]/2,dbeta)
  X, Y = np.meshgrid(x, y)
  temp_x = np.reshape(X,-1)
  temp_y = np.reshape(Y,-1)
  our_ff = evaluate_flat_field(temp_x, temp_y, q)
  god_ff = god.flat_field(temp_x, temp_y, god.flat_field_parameters())
  return 100*np.sqrt(np.mean(((our_ff-god_ff)/god_ff)**2))

def rms_error(flux_estimate, true_flux):
  return 100*np.sqrt(np.mean(((flux_estimate-true_flux)/true_flux)**2))

def average_over_ff(func, args):
  FoV = pdic['FoV']
  nalpha, nbeta = pdic['ff_samples']
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
#********** Saving pickles for plotting **************
#*****************************************************

def saveout_flat_fields(q, iteration_number, strategy):
  FoV = pdic['FoV']
  x = np.arange(-FoV[0]/2,FoV[0]/2,0.01)
  y = np.arange(-FoV[1]/2,FoV[1]/2,0.01)
  X, Y = np.meshgrid(x, y)
  # Have to reshape so that evaluate_flat_field() works 
  temp_x = np.reshape(X,-1)
  temp_y = np.reshape(Y,-1)
  temp_z = evaluate_flat_field(temp_x,temp_y,q)
  our_ff = np.reshape(temp_z, (len(X[:,0]),len(X[0,:])))
  god_q = god.flat_field_parameters()
  temp_z = god.flat_field(temp_x,temp_y)
  god_ff = np.reshape(temp_z, (len(X[:,0]),len(X[0,:])))
  dic = {}
  dic['x'] = X
  dic['y'] = Y
  dic['god_ff'] = god_ff
  dic['our_ff'] = our_ff
  dic['iteration_number'] = iteration_number
  filename = '%s/%s/%04d_ff.p' % (directory_path, strategy, iteration_number)
  pickle.dump(dic, open(filename, "wb"))

def coverage(obs_cat, strategy):
  dic = {}
  sky_limits = pdic['sky_limits']
  dic['number_stars'] = int(np.around(pdic['density_of_stars']*(sky_limits[1]-sky_limits[0]) * (sky_limits[3]-sky_limits[2])))
  dic['k'] = obs_cat.k 
  dic['strategy'] = strategy
  filename = "%s/%s_coverage.p" % (directory_path, strategy)
  pickle.dump(dic, open(filename, "wb" )) 
  return 0

def save_camera(sky_catalog, measured_catalog, inside_FoV, pointing, orientation, verbose=None):
  if verbose: print "Writing out camera pickle..."
  FoV = pdic['FoV']
  x_min = -FoV[0]/2; y_min = -FoV[1]/2
  x_max = FoV[0]/2; y_max = FoV[1]/2
  x = np.array([x_min, x_min, x_max, x_max, x_min])
  y = np.array([y_min, y_max, y_max, y_min, y_min])
  alpha, beta = fp2sky(x,y,pointing, orientation)    
  dic = {}
  dic['measured_catalog.x'] = measured_catalog.x
  dic['measured_catalog.y'] = measured_catalog.y
  dic['sky_catalog.alpha'] = sky_catalog.alpha
  dic['sky_catalog.beta'] = sky_catalog.beta
  dic['pointing'] = pointing
  dic['orientation'] = orientation
  dic['sky'] = pdic['sky_limits']
  dic['fp_x'] = x
  dic['fp_y'] = y
  dic['fp_alpha'] = alpha
  dic['fp_beta'] = beta
  dic['inside_FoV'] = inside_FoV
  filename = directory_path + '/camera_image.p'
  pickle.dump(dic, open(filename, "wb"))
  if verbose: print "...done!"

def invvar_saveout(obs_cat):
  dic = {}
  dic['counts'] = obs_cat.counts
  dic['true_invvar'] = obs_cat.gods_invvar
  dic['reported_invvar'] = obs_cat.invvar
  filename = directory_path + '/invvar.p'  
  print filename
  pickle.dump(dic, open(filename, "wb"))
  return 0
# constants
god_mean_flat_field = average_over_ff(god.flat_field,god.flat_field_parameters())