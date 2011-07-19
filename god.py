#!/usr/bin/env python
# -*- coding: utf-8 -*-
  
import numpy as np
import pickle
import functions as f
import default_parameters
import scipy.optimize as sci

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

def generate_magnitudes(m_min,m_max,powerlaw,size):
  # Randomly generate source magnitudes according to probability distribution
  r = np.random.uniform(size=size)
  m = (1.0/powerlaw) * np.log10((10.0**(powerlaw*m_max)-10.0**(powerlaw*m_min)) * r+ 10**(powerlaw*m_min))
  return m
  
class SourceCatalog:
    def __init__(self, params, out_dir, size, m_min, m_max, powerlaw, limits, seed):
      np.random.seed(seed)
      self.k = np.arange(size).astype(int)
      self.alpha = np.random.uniform(low=limits[0],high=limits[1],size=size) # α
      self.beta = np.random.uniform(low=limits[2],high=limits[3],size=size) # β
      self.mag = generate_magnitudes(m_min,m_max,powerlaw,size)
      self.size = size
      self.flux = f.mag2flux(self.mag)
      self.epsilon = np.random.uniform(0,params['epsilon_max'], size=size)
      
def create_catalog(params, out_dir, plots=None, verbose=False):
  if verbose: print "Generating God's Catalog..."
  M = params['density_of_stars']
  m_min = params['m_min']
  m_max = params['m_max']
  powerlaw = params['powerlaw']
  limits = params['sky_limits']
  seed = params['seed']
  
  # Calculate total number of stars in catalog
  number_stars = M * (limits[1]-limits[0])*(limits[3]-limits[2])
  # Create catalog
  catalog = SourceCatalog(params, out_dir, number_stars, m_min, m_max,powerlaw, limits, seed)
  if verbose: print "...done!"

  # Print out catalog for plotting
  if plots:
    if verbose: print "Saving source catalog data..."
    pickle_dic = {}
    pickle_dic['alpha'] = catalog.alpha
    pickle_dic['beta'] = catalog.beta
    pickle_dic['mag'] = catalog.mag
    pickle_path = out_dir+'/source_catalog.p'
    pickle.dump(pickle_dic, open(pickle_path, "wb"))
    if verbose: print "...done!"
  
  return catalog
  # catalog.ID, catalog.mag, catalog.alpha, catalog.beta, catalog.size
  
