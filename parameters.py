#!/usr/bin/env python
# -*- coding: utf-8 -*-

def mag_at_ten_sigma():
  return 22.5

def eta():
  return 0.01

def m_min():
  return 19.5

def m_max():
  return 22.5

def density_of_stars():
  # density of stars (all magnitudes)
  return 50

def powerlaw():
  #B in log10(dN/dm) = A + B*m
  return 0.25
  
def FoV():
  # Camera's field-of-view
  return [0.7,0.7] # [Δα, Δβ]
  
def sky_limits():
  # limits of spacial/angular coordinates 
  return [0.0,10.0,0.0,10.0] #[α_min, α_max, β_min, β_max]
  
def ff_samples():
  return [100,100]
  
flat_field_order = 9

epsilon_max = 1