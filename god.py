#!/usr/bin/env python
# -*- coding: utf-8 -*-
  
import numpy as np
import math
import sys
import matplotlib.pylab as plt

def flat_field(x,y):
  # God's flat_field
  a = 1
  b = 0
  c = 0
  d = 0.5
  e = -1
  return a + b*x + c*y + d*x**2 + e*y**2

def generate_magnitudes(m_min,m_max,powerlaw,size):
  # Randomly generate source magnitudes according to probability distribution
  r = np.random.uniform(size=size)
  m = (1.0/powerlaw)*np.log10((10.0**(powerlaw*m_max)-10.0**(powerlaw*m_min))*r+10**(powerlaw*m_min))
  return m
  
  
def create_catalog(M, m_min, m_max, powerlaw, sky, limits, seed, plots, verbose):
  # M = density of stars 
  # m_min = minimum magnitude (i.e. saturation limit)
  # m_max = maximum magnitude limit per dither (i.e. sensitivity limit 10σ)
  # powerlaw = B in log10(dN/dm) = A + B*m
  # sky = create a catalog on a square (='sq') or sphere (='sp')
  # seed = seed for random number generation, providing same seed generate the same sky
  # limits = limits of spacial/angular coordinates [α_min, α_max, β_min, β_max]

  if sky != 'sq':
    print "Exiting - only square sky implemented!"
    sys.exit() 
  
  # Calculate total number of stars and create catalog
  number_stars = M * (limits[1]-limits[0])*(limits[3]-limits[2])
  catalog = np. zeros((number_stars, 4))
  
  # Seed the random number generator
  np.random.seed(seed)
  
  # Add source ID to catalog
  catalog[:,0] = np.arange(number_stars)
  
  size = (number_stars)
  # Add star magnitudes to catalog with defined probability distribution
  catalog[:,1] = generate_magnitudes(m_min, m_max, powerlaw, size)
  # Add random coordinates
  catalog[:,2] = np.random.uniform(low=limits[0],high=limits[1],size=size) # α
  catalog[:,3] = np.random.uniform(low=limits[2],high=limits[3],size=size) # β
  
  #***********************************************************
  #********************* Plot Catalog ************************
  #***********************************************************
    
  if plots != None:
    # Plot portion of sky
    if verbose != None: print "Plotting portion of sky..."
    plt.figure()
    plt.plot(catalog[0:1000,2],catalog[0:1000,3],'ko',markersize=2)
    plt.xlabel(ur'$\alpha$', fontsize=20)
    plt.ylabel(ur'$\beta$', fontsize=20)
    plt.title("God's Sky")
    plt.savefig("Figures/gods_sky.png")
    if verbose != None: print "...done!"
    
    # Histogram of source magnitude
    plt.figure()
    bin=np.arange(m_min,m_max,0.05)
    if verbose != None: print "Plotting histogram of source magnitude..."
    plt.hist(catalog[:,1],bins=bin, log=True)
    plt.title("%i Sources in Sample" % len(catalog[:,0]))
    plt.xlabel("Source Magnitude")
    plt.ylabel("log(N)")
    plt.savefig("Figures/source_mag_histogram.png")
    if verbose != None: print "...done!"
  
  return catalog
  # catalog = [Star ID, magnitude, α, β]