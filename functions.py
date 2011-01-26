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
  
def evaluate_flat_field(focal_position, q, plots, verbose):
  return 1

def normalize_flat_field(q):
  return q