#!/usr/bin/env python
# -*- coding: utf-8 -*-
  
import numpy as np
import math
import sys
import matplotlib.pylab as plt

def flat_field_parameters():
  return np.array([1, 0, 0, -0.3, 0, -0.5])

def flat_field(x,y):
  # God's flat_field
  par = flat_field_parameters()
  return (par[0] + par[1]*x + par[2]*y + par[3]*x**2 + par[4]*x*y + par[5]*y**2) 

def generate_magnitudes(m_min,m_max,powerlaw,size):
  # Randomly generate source magnitudes according to probability distribution
  r = np.random.uniform(size=size)
  m = (1.0/powerlaw)*np.log10((10.0**(powerlaw*m_max)-10.0**(powerlaw*m_min))*r+10**(powerlaw*m_min))
  return m
  
  
def create_catalog(M, m_min, m_max, powerlaw, limits, seed, plots=None, verbose=None):
  # M = density of stars 
  # m_min = minimum magnitude (i.e. saturation limit)
  # m_max = maximum magnitude limit per dither (i.e. sensitivity limit 10σ)
  # powerlaw = B in log10(dN/dm) = A + B*m
  # sky = create a catalog on a square (='sq') or sphere (='sp')
  # seed = seed for random number generation, providing same seed generate the same sky
  # limits = limits of spacial/angular coordinates [α_min, α_max, β_min, β_max]

  
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
    fontsize = 25
    if verbose != None: print "Plotting portion of sky..."
    plt.figure()
    plt.plot(catalog[0:1000,2],catalog[0:1000,3],'ko',markersize=2)
    plt.xlabel(ur'$\alpha$', fontsize=fontsize)
    plt.ylabel(ur'$\beta$', fontsize=fontsize)
    plt.title("God's Sky", fontsize=fontsize)
    ax = plt.gca()
    for tick in ax.xaxis.get_major_ticks():
      tick.label1.set_fontsize(fontsize)
    for tick in ax.yaxis.get_major_ticks():
      tick.label1.set_fontsize(fontsize)
    plt.savefig("Figures/gods_sky.png",bbox_inches='tight',pad_inches=0.)
    if verbose != None: print "...done!"
    
    # Histogram of source magnitude
    plt.figure()
    fontsize = 25
    bin=np.arange(m_min,m_max,0.05)
    if verbose != None: print "Plotting histogram of source magnitude..."
    plt.hist(catalog[:,1],bins=bin, log=True)
    plt.title("%i Sources in Sample" % len(catalog[:,0]), fontsize=fontsize)
    plt.xlabel("Source Magnitude", fontsize=fontsize)
    plt.ylabel("log(N)", fontsize=fontsize)
    ax = plt.gca()
    for tick in ax.xaxis.get_major_ticks():
      tick.label1.set_fontsize(fontsize)
    for tick in ax.yaxis.get_major_ticks():
      tick.label1.set_fontsize(fontsize)

    
    plt.savefig("Figures/source_mag_histogram.png",bbox_inches='tight',pad_inches=0.)
    if verbose != None: print "...done!"
  
  return catalog
  # catalog = [Star ID, magnitude, α, β]