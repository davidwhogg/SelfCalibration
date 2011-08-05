#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Rory Holmes 2011
# Contains all the functions used in the cross-calibration simulation [./cross_cal.py]

import numpy as np
import os
import pickle
import string
import scipy.optimize as sci
# Custom Modules
import god


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
    def __init__(self, params, camera_catalog, inside_FoV): 
      self.size = len(inside_FoV[0])
      self.k = camera_catalog.k[inside_FoV].astype(int)
      self.x = camera_catalog.x[inside_FoV]
      self.y = camera_catalog.y[inside_FoV]
      flat = god.flat_field(params, self.x,self.y)
      self.gods_invvar = self.true_invvar(params, camera_catalog, inside_FoV, flat)
      self.counts = camera_catalog.flux[inside_FoV] * flat + np.random.normal(size=self.size)/np.sqrt(self.gods_invvar)
      self.invvar = self.reported_invvar(params)

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
      
    def reported_invvar(self, params):
      sky_unc = 0.1*mag2flux(params['mag_at_ten_sigma'])
      var = (sky_unc**2 + (params['eta']**2) * self.counts**2)
      return 1. / var

    def true_invvar(self, params, camera_catalog, inside_FoV, flat):
      epsilon = camera_catalog.epsilon[inside_FoV]
      sky_unc = 0.1*mag2flux(params['mag_at_ten_sigma'])
      flux = camera_catalog.flux[inside_FoV]
      var = (sky_unc**2 * (1. + epsilon**2) + (params['eta']**2) * flat**2 * flux**2)
      return 1. / var

def single_image(params, sky_catalog, pointing, orientation, data_dir, plots=None, verbose=None):
  if verbose: print "Converting sky catalog to focal plane coordinates..."
  camera_catalog = CameraCatalog(sky_catalog, pointing, orientation)
  if verbose: print "...done!"
  if verbose: print "Finding stars within camera FoV..."
  FoV = params['FoV']
  x_min = -0.5*FoV[0]; y_min = -0.5*FoV[1]
  x_max = 0.5*FoV[0]; y_max = 0.5*FoV[1]
  inside_FoV = np.where((x_min<camera_catalog.x) & (camera_catalog.x<x_max) & (y_min<camera_catalog.y) & (camera_catalog.y<y_max))
  if verbose: print "...done!"

  if verbose: print "Measuring stars within FoV..."
  measured_catalog = MeasuredCatalog(params, camera_catalog, inside_FoV)
  if verbose: print "...done!"
  one_camera_file = os.path.exists((data_dir + '/camera_image.p'))
  if plots and (one_camera_file != True) and (orientation > 30) and (pointing[0] > -1) and (pointing[0] < 1) and (pointing[1] > -1) and (pointing[1] < 1): save_camera(params, sky_catalog, measured_catalog, inside_FoV, pointing, orientation, data_dir, verbose = verbose)
  return measured_catalog
  # measured_sources  *.size, *.k, *.flux, *.invvar, *.x, *.y

#*****************************************************
#**************** Survey Functions *******************
#*****************************************************

def survey(params, sky_catalog, survey_file, data_dir, plots=None, verbose=None):  
  if verbose: print "Loading survey..."
  pointing = np.loadtxt(survey_file)
  number_pointings = len(pointing[:,0])  
  if verbose: print "...done!"
  if verbose: print "Surveying sky..."
  obs_cat = None
  for i in range(number_pointings):
    si = single_image(params, sky_catalog, [pointing[i,1],pointing[i,2]], pointing[i,3], data_dir, plots=plots, verbose = verbose)
    if obs_cat is None:
      obs_cat = si
    else:
      obs_cat.append(si)
  if verbose: print "...done!"
  return obs_cat

#*****************************************************
#**************** Ubercal Functions ******************
#*****************************************************

def ubercalibration(params, observation_catalog, sky_catalog, strategy, out_dir, data_dir, plots = None):
  order = params['flat_field_order']
  q = np.array([1])
  stop_condition = params['stop_condition']
  max_iterations = params['max_iterations']
  chi2 = 1e9
  old_chi2 = 1e10
  iteration_number = 0
  next_plot_iteration = 1 
  if plots: saveout_flat_fields(params, q, iteration_number, data_dir)
  while ((abs(chi2 - old_chi2) > stop_condition) and (iteration_number < max_iterations)):
    iteration_number += 1
    temp_chi2 = chi2
    s, s_invvar = s_step(params, observation_catalog,q)
    q, q_invvar, chi2 = q_step(params, observation_catalog, s, order,iteration_number,plots=plots)
    old_chi2 = temp_chi2
     # Calculate rms error in stars
    indx = [s != 0]
    rms = rms_error(s[indx],sky_catalog.flux[indx])
    bdness = badness(params, q)
    bdness_bestfitff = badness_bestinbasis(params, q, out_dir) 
    print "%i: RMS = %.6f %%, Badness = %0.6f %%, BestInBasis_Badness = %0.6f %%, chi2 = %0.2f (%i)" % (iteration_number, rms, bdness, bdness_bestfitff, chi2, observation_catalog.size)
    print q
    if (plots and (iteration_number == next_plot_iteration)) or (abs(chi2 - old_chi2) < stop_condition): 
      saveout_flat_fields(params, q, iteration_number, data_dir)
      next_plot_iteration *= 2
  return np.array([iteration_number, rms, bdness, bdness_bestfitff, chi2])

def evaluate_flat_field_functions(x, y, order):
  L = (order + 1)*(order + 2)/2
  g = np.zeros((len(x),L))
  l = 0
  for n in range(order+1):
    for m in range(n+1):
      g[:,l] = (x**(n-m)) * (y**m)
      l += 1 
  return g

def evaluate_flat_field(params, x, y, q):  
  # Calculate required order
  order = int(np.around(np.sqrt(0.25+2*len(q))-1.5))
  assert(len(q) == ((order + 1) * (order + 2) / 2))
  g = evaluate_flat_field_functions(x, y, order)
  return np.dot(g, q)

def normalize_flat_field(params, q):
  fit_mean = average_over_ff(params, evaluate_flat_field, (q))
  god_mean = average_over_ff(params, god.flat_field,(god.flat_field_parameters()))
  return (q * god_mean/fit_mean)
  
def average_over_ff(params, func, args):
  FoV = params['FoV']
  nalpha, nbeta = params['ff_samples']
  dalpha = FoV[0]/(nalpha-1)
  dbeta = FoV[1]/(nbeta-1)
  x = np.arange(-FoV[0]/2+dalpha/2,FoV[0]/2,dalpha)
  y = np.arange(-FoV[1]/2+dbeta/2,FoV[1]/2,dbeta)
  X, Y = np.meshgrid(x, y)
  temp_x = np.reshape(X,-1)
  temp_y = np.reshape(Y,-1)
  temp_ff = func(params, temp_x, temp_y, args)
  ff = np.reshape(temp_ff, (len(X[:,0]),len(X[0,:])))
  return np.mean(ff)

def s_step(params, obs_cat, q):
  ff = evaluate_flat_field(params, obs_cat.x, obs_cat.y, q)
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

def q_step(params, obs_cat, s, order, iteration_number, plots=None):
  g = evaluate_flat_field_functions(obs_cat.x,obs_cat.y, order)
  ss = s[obs_cat.k]
  q_invvar = np.dot(np.transpose(g)*ss*ss*obs_cat.invvar,g)
  numerator = np.sum(np.transpose(g)*ss*obs_cat.counts*obs_cat.invvar,axis=1)
  q = np.dot(np.linalg.inv(q_invvar),numerator)
  q = normalize_flat_field(params, q)
  np.set_printoptions(precision=2) 
  chi2 = np.sum((np.dot(g,q) * ss - obs_cat.counts)**2 * obs_cat.invvar)
  return (q, q_invvar, chi2)

#*****************************************************
#*************** Best Fit of God's *******************
#*****************************************************

def compare_flats(a, Z_true, g, X, Y):
  error = np.sum((Z_true - np.dot(g,a))**2)
  #print error
  return error

def bestfit_ff(params, out_dir):
  order = params['flat_field_order']
  ff_samples = params['ff_samples']
  FoV = params['FoV']
  temp_x = np.linspace(-0.5*FoV[0], 0.5*FoV[0], ff_samples[0])
  temp_y = np.linspace(-0.5*FoV[1], 0.5*FoV[1], ff_samples[1])
  X, Y = np.meshgrid(temp_x,temp_y)
  x = np.reshape(X,-1)
  y = np.reshape(Y,-1)
  g = evaluate_flat_field_functions(x,y, order)
  god_ff = god.flat_field(params,x,y)
  a = np.zeros((order + 1)*(order + 2)/2)
  a[0] = 1
  a[3] = -0.2
  a[5] = 0.5
  print "Fitting god's flat-field with basis..."
  fitted_parameters = sci.fmin_bfgs(compare_flats, a, args = (god_ff, g, x, y), gtol = 1e-8, maxiter = 1e16)
  print "Fitted Parameters: ", fitted_parameters
  "...done!"
  fitted_parameters = normalize_flat_field(params, fitted_parameters)
  bestfit_ff = np.reshape(evaluate_flat_field(params, x, y, fitted_parameters),
(len(X[:,0]),len(X[0,:])))
  dic = {}
  dic['x'] = X
  dic['y'] = Y
  dic['bestfit_ff'] = bestfit_ff
  dic['fit_parameters'] = fitted_parameters
  filename = '%s/bestfit_ff.p' % (out_dir)
  pickle.dump(dic, open(filename, "wb"))

#*****************************************************
#*************** Analysis Functions ******************
#*****************************************************

def badness(params, q):
  FoV = params['FoV']
  nalpha, nbeta = params['ff_samples']
  x = np.linspace(-FoV[0]/2,FoV[0]/2,nalpha)
  y = np.linspace(-FoV[1]/2,FoV[1]/2,nbeta)
  X, Y = np.meshgrid(x, y)
  temp_x = np.reshape(X,-1)
  temp_y = np.reshape(Y,-1)
  our_ff = evaluate_flat_field(params, temp_x, temp_y, q)
  god_ff = god.flat_field(params,temp_x, temp_y)
  return 100*np.sqrt(np.mean(((our_ff-god_ff)/god_ff)**2))

def badness_bestinbasis(params, q, out_dir):
  FoV = params['FoV']
  nalpha, nbeta = params['ff_samples']
  x = np.linspace(-FoV[0]/2,FoV[0]/2,nalpha)
  y = np.linspace(-FoV[1]/2,FoV[1]/2,nbeta)
  X, Y = np.meshgrid(x, y)
  temp_x = np.reshape(X,-1)
  temp_y = np.reshape(Y,-1)
  our_ff = evaluate_flat_field(params, temp_x, temp_y, q)
  pickle_dic = pickle.load(open(out_dir+'/bestfit_ff.p'))
  bestfit_ff_parameters = pickle_dic['fit_parameters']
  bestfit_ff = evaluate_flat_field(params, temp_x, temp_y, bestfit_ff_parameters)
  return 100*np.sqrt(np.mean(((our_ff-bestfit_ff)/bestfit_ff)**2))

def rms_error(flux_estimate, true_flux):
  return 100*np.sqrt(np.mean(((flux_estimate-true_flux)/true_flux)**2))


#*****************************************************
#********** Saving pickles for plotting **************
#*****************************************************

def saveout_flat_fields(params, q, iteration_number, data_dir):
  FoV = params['FoV']
  ff_samples = params['ff_samples']
  x = np.linspace(-FoV[0]/2,FoV[0]/2,ff_samples[0])
  y = np.linspace(-FoV[1]/2,FoV[1]/2,ff_samples[1])
  X, Y = np.meshgrid(x, y)
  # Have to reshape so that evaluate_flat_field() works 
  temp_x = np.reshape(X,-1)
  temp_y = np.reshape(Y,-1)
  temp_z = evaluate_flat_field(params, temp_x,temp_y,q)
  our_ff = np.reshape(temp_z, (len(X[:,0]),len(X[0,:])))
  god_q = god.flat_field_parameters()
  temp_z = god.flat_field(params,temp_x,temp_y)
  god_ff = np.reshape(temp_z, (len(X[:,0]),len(X[0,:])))
  dic = {}
  dic['x'] = X
  dic['y'] = Y
  dic['our_ff'] = our_ff
  dic['iteration_number'] = iteration_number
  filename = '%s/FF/%0*d_ff.p' % (data_dir, 3, iteration_number)
  pickle.dump(dic, open(filename, "wb"))


def coverage(params, obs_cat, strategy, data_dir):
  dic = {}
  sky_limits = params['sky_limits']
  dic['number_stars'] = int(np.around(params['density_of_stars']*(sky_limits[1]-sky_limits[0]) * (sky_limits[3]-sky_limits[2])))
  dic['k'] = obs_cat.k 
  dic['strategy'] = strategy
  filename = data_dir + "/coverage.p"
  pickle.dump(dic, open(filename, "wb" )) 
  return 0

def save_camera(params, sky_catalog, measured_catalog, inside_FoV, pointing, orientation, data_dir, verbose=None):
  if verbose: print "Writing out camera pickle..."
  FoV = params['FoV']
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
  dic['sky'] = params['sky_limits']
  dic['fp_x'] = x
  dic['fp_y'] = y
  dic['fp_alpha'] = alpha
  dic['fp_beta'] = beta
  dic['inside_FoV'] = inside_FoV
  filename = data_dir + '/camera_image.p'
  pickle.dump(dic, open(filename, "wb"))
  if verbose: print "...done!"

def invvar_saveout(obs_cat, data_dir):
  dic = {}
  dic['counts'] = obs_cat.counts
  dic['true_invvar'] = obs_cat.gods_invvar
  dic['reported_invvar'] = obs_cat.invvar
  filename = data_dir + '/invvar.p'  
  print filename
  pickle.dump(dic, open(filename, "wb"))
  return 0

