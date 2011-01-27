#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import parameters

def mag2flux(mag):
  return np.power(10,(-0.4*(mag-22.5)))
  
def flux2mag(flux):
  return 22.5-2.5*np.log10(flux)
  
def flux_uncertainty_variance(mag, eta):
  sky_error= 0.1*mag2flux(parameters.mag_at_ten_sigma())
  return (np.power(sky_error,2) + (eta**2)*np.power(mag2flux(mag),2))
  
def evaluate_flat_field(x, y, q):
  # Function returns flat field at focal plane position specified. The function also calculates the order of the flat field equation from the number of flat field parameters
  
  # Calculate required order
  order = int(np.around(np.sqrt(0.25+2*len(q))-1.5))
  assert(len(q) == ((order + 1) * (order + 2) / 2)) # assert to break code if condition not met
 
  # Calculate flat-field
  l = 0
  total = 0
  for n in range(order+1):
    for m in range(n+1):
      total += q[l] * (x**(n-m)) * (y**m)
      l += 1 
  
  return total # Return flat-field at specified [x,y] focal plane position

def normalize_flat_field(q):
  return q

def s_step(obs_cat, q):
  # First step in ubercalibration. This function calculates the new flux estimates.
  # obs_cat = observation_catalog = [pointing ID, star ID, observed_flux, observed_invvar, focal plane x, focal plane y]
  
  # Calculate flat-field based on source positions and current flat field parameters 
  
  ff = evaluate_flat_field(obs_cat[:,4], obs_cat[:,5], q)
  fcss = ff * obs_cat[:,2] * obs_cat[:,3]  
  ffss = ff * ff * obs_cat[:,3]
  star_ID = np.around(obs_cat[:,1]).astype('int')
  max_star = np.max(star_ID)
  s = np.zeros(max_star)
  s_invvar = np.zeros(max_star)
  for ID in range(max_star):
    indx = (star_ID == ID)
    denominator = np.sum(ffss[indx])
    s_invvar[ID] = denominator
    if denominator > 0.:
      s[ID] = np.sum(fcss[indx])/denominator
  return (s, s_invvar)
  
#def q_step(obs_cat, s, order):
  
  