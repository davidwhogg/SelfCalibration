#!/usr/bin/env python
# -*- coding: utf-8 -*-
  
import numpy as np
import pickle
import functions as f
import default_parameters
import scipy.optimize as opt
import sys

# magic numbers
def flat_field_parameters():
  np.random.seed(1) # same seed => same FF
  return np.append(np.array([1, 0, 0, -0.2, 0, -0.5]),
    (1e-4) * np.random.uniform(size = 256))

def flat_field(params,x,y,par = flat_field_parameters()):
  # Large Scale
  ff = (par[0] + par[1]*x + par[2]*y + par[3]*x**2 + par[4]*x*y + par[5]*y**2)
  # Smaller Scale
  k = 6 
  fov = params['FoV']
  no_loops = int(np.sqrt((len(par) - 6)/4))
  for nx in range(no_loops):
    kx = nx * np.pi / fov[0]
    ckx = np.cos(kx * x)
    skx = np.sin(kx * x)
    for ny in range(no_loops):
      ky = ny * np.pi / fov[1]
      cky = np.cos(ky * y)
      sky = np.sin(ky * y)
      ff += par[k] * ckx * cky
      ff += par[k+1] * ckx * sky
      ff += par[k+2] * skx * cky
      ff += par[k+3] * skx * sky
      k += 4
  return ff

def saveout_god_ff(params, out_dir):
  flat_field_order = params['flat_field_order']
  ff_samples = params['ff_samples']
  FoV = params['FoV']
  x = np.linspace(-0.5*FoV[0], 0.5*FoV[0], ff_samples[0])
  y = np.linspace(-0.5*FoV[1], 0.5*FoV[1], ff_samples[1])
  X, Y = np.meshgrid(x,y)
  god_ff = flat_field(params, X, Y)
  dic = {}
  dic['x'] = X
  dic['y'] = Y
  dic['god_ff'] = god_ff
  filename = '%s/god_ff.p' % (out_dir)
  pickle.dump(dic, open(filename, "wb"))

def power_law(A, m):
 temp = np.zeros(len(m))
 for indx in range(len(A)):
  temp += A[indx]*m**indx
 return 10**temp 

def error(a, d, m):
 error = np.sum(((d-power_law(a, m))/d)**2)
 return error

def prob_dist(m_min, m_max, a, size):
  # Randomly generate source magnitudes according to probability distribution
  r = np.random.uniform(size=size)
  m = (1.0/a) * np.log10((10.0**(a*m_max)-10.0**(a*m_min)) * r+ 10**(a*m_min))
  return m

def generate_magnitudes(params, number_sources, out_dir, plots = False):
  A = params['powerlaw_constants']
  area = (params['sky_limits'][1]-params['sky_limits'][0])*(params['sky_limits'][3]-params['sky_limits'][2])
  # fit for dN/dm = B[0] + B[1]*m
  mag_range = np.arange(params['m_min']+0.25, params['m_max'], 0.5)
  B = opt.fmin(error, np.array([0,0]), args = (power_law(A, mag_range), mag_range), xtol=1e-14, maxiter = 1e16, maxfun = 1e16)
  B = np.append(B, [0])
  # Shift above dN/dm = a + b*m + c*m**2 line
  chck = 1.
  while chck > 0.1:
   bad = power_law(A, mag_range) - power_law(B, mag_range)
   ii = np.where(bad > 0.)
   B[0] += 0.001
   chck = np.sum(bad[ii[0]])
  
  total_sources = np.sum(power_law(A, mag_range))*area
  if total_sources*params['useful_fraction'] < number_sources:
   print "Error! You want more sources than there are in the sky..."
   sys.exit()
  
  mag = prob_dist(params['m_min'], params['m_max'], B[1],params['useful_fraction']*total_sources)
  
  # Use rejection method to get correct probability distribution
  c = np.max(power_law(A, mag)/power_law(B, mag))
  difference = power_law(A, mag)/(power_law(B, mag)*c)
  randoms = np.random.uniform(size = len(difference)) 
  ii = np.zeros((2,2))
  while len(ii[0]) > 1:
   ii = np.where(randoms > difference)
   k = len(ii[0])
   mag[ii] = prob_dist(params['m_min'], params['m_max'], B[1],len(ii[0]))
   randoms = np.random.uniform(size = len(ii[0]))
   difference = power_law(A, mag[ii])/(power_law(B, mag[ii])*c)
  
  selected_sources = np.sort(mag)[0:int(number_sources)]
  # Print out catalog for plotting
  if plots:
    survey_dic = {}
    survey_dic['mag'] = selected_sources
    survey_dic['all_sources'] = mag
    survey_dic['fit_mag'] = np.arange(14, 26, 0.01)
    survey_dic['fit_dens'] = power_law(A, np.arange(14, 26, 0.01))*area
    pickle_path = out_dir+'/source_catalog.p'
    pickle.dump(survey_dic, open(pickle_path, "wb"))
  return selected_sources

class SourceCatalog:
    def __init__(self, params, out_dir, number_sources, plots = False):
      np.random.seed(params['seed'])
      self.k = np.arange(number_sources).astype(int)
      self.alpha = np.random.uniform(low=params['sky_limits'][0],high=params['sky_limits'][1],size=number_sources) # α
      self.beta = np.random.uniform(low=params['sky_limits'][2],high=params['sky_limits'][3],size=number_sources) # β
      self.mag = generate_magnitudes(params, number_sources, out_dir, plots = plots)
      self.size = number_sources
      self.flux = f.mag2flux(self.mag)
      self.epsilon = np.random.uniform(0,params['epsilon_max'], size=number_sources)
      
def create_catalog(params, out_dir, plots=False, verbose=False):
  if verbose: print "Generating God's Catalog..."
  # Calculate total number of stars in catalog
  number_sources = params['density_of_stars'] * (params['sky_limits'][1]-params['sky_limits'][0])*(params['sky_limits'][3]-params['sky_limits'][2])
  # Create catalog
  catalog = SourceCatalog(params, out_dir, number_sources, plots = plots)
  if verbose: print "...done!"
  if plots:
    pickle_path = out_dir+'/source_catalog.p'
    survey_dic = pickle.load(open(pickle_path))
    survey_dic['alpha'] = catalog.alpha
    survey_dic['beta'] = catalog.beta
    pickle_path = out_dir+'/source_catalog.p'
    pickle.dump(survey_dic, open(pickle_path, "wb"))
  return catalog
  # catalog.ID, catalog.mag, catalog.alpha, catalog.beta, catalog.size
  
