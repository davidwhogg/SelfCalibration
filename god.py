#!/usr/bin/env python
# -*- coding: utf-8 -*-
  
import numpy as np
import math
import matplotlib.pylab as plt
import pickle

import functions as f
from master import init_func 
pdic, directory_path, temp1, temp2 = init_func() # import parameter database from main module

# magic numbers
def flat_field_parameters():
  return np.append(np.array([1, 0, 0, -0.2, 0, -0.5]),
    (0.003 / 16.) * np.array([ -2.46640144e-01,   7.98451075e-01, 4.15607201e-01,
        -4.70986300e-01,   8.48211144e-01,  -1.45938404e+00,
         1.48463628e-01,   8.17257101e-01,  -2.60560568e-01,
         1.55530761e-02,  -1.09873618e-01,  -6.52229300e-01,
        -8.17172195e-01,  -6.83910165e-01,   2.91446321e-01,
        -3.30006370e-01,   5.41434486e-02,  -9.09158505e-02,
         1.45073016e+00,  -2.62954213e+00,  -1.26893671e+00,
         4.82147110e-01,   4.52229224e-01,  -1.65617743e+00,
         1.45802412e+00,  -1.63648027e+00,  -8.34149101e-01,
        -9.77390857e-02,   2.14295348e-01,   1.36684839e+00,
        -9.33976099e-01,   4.52120550e-01,   3.56813818e-01,
         8.87511851e-01,   7.86754603e-02,  -5.16448573e-01,
         5.67890585e-01,   5.40594950e-01,  -5.58005823e-01,
        -5.43919673e-01,   1.28633646e+00,  -1.14865468e+00,
        -4.65055531e-01,   9.59426420e-01,  -4.67244948e-01,
        -1.33755986e+00,   1.57364497e+00,  -1.86430907e+00,
         1.47085593e+00,  -1.13733344e+00,   1.18986782e-01,
        -3.94974777e-01,  -7.95232892e-01,  -1.38971887e+00,
         1.37771983e+00,   8.33680274e-01,  -1.47274333e+00,
         1.15160848e+00,  -3.64293316e-01,  -3.35002574e-01,
        -2.48736707e-01,   1.76853593e-01,  -5.10398674e-01,
         1.39728196e+00,   1.50372453e-01,  -8.33288140e-01,
         1.20570437e+00,  -6.58470259e-01,  -3.19115018e-01,
        -2.27777233e-01,   9.88676751e-01,   3.96942286e-01,
        -1.14943975e-01,   9.36292074e-02,   2.98004622e-01,
        -1.44645332e+00,  -2.37799423e-04,  -1.94516831e+00,
        -4.71249238e-01,   4.85317021e-01,  -1.06633432e+00,
        -7.29061834e-01,  -7.33823059e-01,   3.91340951e-01,
        -1.49986434e+00,  -3.15765143e-01,  -1.73605656e+00,
         7.30753634e-01,  -2.76776026e-01,  -1.56526642e-01,
         5.46853528e-01,   8.55326967e-02,  -2.62225693e-01,
        -1.11248990e-01,  -3.05262711e-01,   2.01564740e-01,
        -1.50518297e+00,  -1.62010376e+00,  -1.30051218e+00,
         6.46755325e-01,  -1.01408576e+00,   2.08391768e+00,
         7.84008734e-01,  -1.23928888e+00,  -1.15179742e+00,
        -2.01366835e+00,  -1.39564290e-01,   1.08175535e+00,
        -3.64176383e-01,  -2.40227761e-01,  -1.10681091e-01,
        -2.71131092e-01,  -8.67828529e-01,   1.61204842e+00,
        -6.27667494e-01,   1.18222300e+00,  -1.62734327e+00,
        -6.64306227e-01,   1.80972465e+00,  -2.21554881e-02,
         1.56562142e+00,  -1.90106161e+00,   2.78195257e-01,
         3.16278370e-01,   4.86245659e-01,  -1.22536446e+00,
        -8.85969027e-01,   5.98483215e-01,   4.81636317e-01,
        -1.52289276e+00,   1.10097058e+00,   9.19964053e-02,
         4.37198026e-01,   8.86136398e-01,   2.51617981e+00,
         4.71431834e-01,  -2.09274864e-01,  -7.74074526e-01,
         5.36217421e-01,  -1.34475975e+00,  -7.88510125e-01,
        -5.74502581e-01,   3.31090084e-01,  -2.38729072e-01,
        -1.77847087e+00,  -1.61094388e-01,   1.75468712e+00,
        -2.07971838e+00,   4.62452449e-01,   1.59985887e+00,
         1.36801268e+00,   7.07677848e-01,   3.73954642e-02,
         2.40363867e+00,  -1.19501856e-02,   1.65095079e-01,
         8.47791981e-01,  -7.95697277e-01,  -1.05970895e+00,
         9.04195222e-01,   1.74832606e+00,   2.97715939e-01,
         2.03838511e+00,   4.75460629e-01,   1.85219525e-01,
         1.60181763e-01,   1.51497929e+00,   1.66015944e+00,
         8.28551543e-02,  -1.12971888e+00,  -3.38545531e-01,
        -1.50687961e-01,   1.18792895e+00,   2.36463105e-01,
        -5.40825858e-01,  -1.30226928e+00,  -3.72279622e-01,
        -1.40340128e-02,  -1.82565834e+00,  -3.82187025e-01,
         2.67037744e-01,   5.66950838e-01,   6.77644161e-01,
         1.01787440e+00,   8.21296613e-01,   1.42132007e+00,
         8.28210856e-01,   2.16844870e-01,   1.25936836e-01,
         1.65362976e-01,  -6.71116024e-01,  -1.35098450e+00,
         6.33290724e-01,  -2.35895832e-01,   1.25877054e-01,
        -1.32839395e+00,  -1.21716079e+00,  -1.43434583e+00,
         9.77545338e-01,   7.24054063e-01,  -2.30030003e-02,
         1.33408426e+00,  -5.39331399e-01,  -3.23801981e-01,
        -5.03217410e-01,   3.31567328e-01,   1.29523041e+00,
         9.41461039e-01,   4.22569332e-01,  -1.33850832e+00,
        -1.51889310e+00,   1.88904111e-02,  -8.01725336e-01,
        -5.92959289e-01,  -1.65831899e+00,   1.43462485e+00,
         8.73636918e-01,   4.71088802e-01,  -1.18486217e-01,
         4.07890745e-01,   1.06523283e+00,  -1.50455849e+00,
         3.08367670e-01,  -7.42079255e-01,   9.35106561e-03,
         1.31435647e+00,   1.07063199e+00,   1.25015898e+00,
        -1.19967374e+00,  -1.86221629e-01,   1.24486453e+00,
         7.51979378e-01,  -9.84257037e-01,  -6.22882840e-01,
         2.85946452e-01,   1.72643322e-01,  -3.50057065e-01,
        -1.20065938e+00,  -9.74609731e-01,  -1.94489260e+00,
        -1.04200832e-01,   7.05679235e-01,   8.40274731e-01,
         3.82345240e-01,  -4.13907641e-01,  -8.33327666e-01,
        -3.34106966e+00,  -2.55433786e-01,  -2.01373211e-01,
         1.15798517e+00,   1.24301783e+00,   1.05456054e+00,
        -6.13243663e-01,  -5.12256121e-01,  -1.65605714e-01,
         7.52010421e-01]))


def flat_field(x,y,par = flat_field_parameters()):
  # illumination part
  ff = (par[0] + par[1]*x + par[2]*y + par[3]*x**2 + par[4]*x*y + par[5]*y**2)
  # detector part
  k = 6
  fov = pdic['FoV']
  for nx in range(8):
    kx = nx * np.pi / fov[0]
    ckx = np.cos(kx * x)
    skx = np.sin(kx * x)
    for ny in range(8):
      ky = ny * np.pi / fov[1]
      cky = np.cos(ky * y)
      sky = np.sin(ky * y)
      ff += par[k] * ckx * cky
      ff += par[k] * ckx * sky
      ff += par[k] * skx * cky
      ff += par[k] * skx * sky
      k += 4
  return ff

def generate_magnitudes(m_min,m_max,powerlaw,size):
  # Randomly generate source magnitudes according to probability distribution
  r = np.random.uniform(size=size)
  m = (1.0/powerlaw)*np.log10((10.0**(powerlaw*m_max)-10.0**(powerlaw*m_min))*r+10**(powerlaw*m_min))
  return m
  
class SourceCatalog:
    def __init__(self, size, m_min, m_max, powerlaw, limits, seed):
      np.random.seed(seed)
      self.k = np.arange(size).astype(int)
      self.alpha = np.random.uniform(low=limits[0],high=limits[1],size=size) # α
      self.beta = np.random.uniform(low=limits[2],high=limits[3],size=size) # β
      self.mag = generate_magnitudes(m_min,m_max,powerlaw,size)
      self.size = size
      self.flux = f.mag2flux(self.mag)
      self.epsilon = np.random.uniform(0,pdic['epsilon_max'], size=size)
      
def create_catalog(seed, plots=None, verbose=False):
  if verbose: print "Generating God's Catalog..."
  M = pdic['density_of_stars']
  m_min = pdic['m_min']
  m_max = pdic['m_max']
  powerlaw = pdic['powerlaw']
  limits = pdic['sky_limits']
  
  # Calculate total number of stars in catalog
  number_stars = M * (limits[1]-limits[0])*(limits[3]-limits[2])
  # Create catalog
  catalog = SourceCatalog(number_stars, m_min, m_max,powerlaw, limits, seed)
  if verbose: print "...done!"

  #**********************************************************
  #********************* Plot Catalog ***********************
  #**********************************************************
    
  if plots:
    if verbose: print "Saving source catalog data..."
    pickle_dic = {}
    pickle_dic['alpha'] = catalog.alpha
    pickle_dic['beta'] = catalog.beta
    pickle_dic['mag'] = catalog.mag
    pickle_path = directory_path+'/source_catalog.p'
    pickle.dump(pickle_dic, open(pickle_path, "wb"))
    if verbose: print "...done!"
  
  return catalog
  # catalog.ID, catalog.mag, catalog.alpha, catalog.beta, catalog.size
  