#!/usr/bin/env python
# -*- coding: utf-8 -*-
  
import numpy as np
import math
import sys
import matplotlib.pylab as plt
import functions as f

def flat_field_parameters():
  return np.array([1, 0, 0, -0.2, 0, -0.5])

def flat_field(x,y,par = flat_field_parameters()):
  return (par[0] + par[1]*x + par[2]*y + par[3]*x**2 + par[4]*x*y + par[5]*y**2) 

def generate_magnitudes(m_min,m_max,powerlaw,size):
  # Randomly generate source magnitudes according to probability distribution
  r = np.random.uniform(size=size)
  m = (1.0/powerlaw)*np.log10((10.0**(powerlaw*m_max)-10.0**(powerlaw*m_min))*r+10**(powerlaw*m_min))
  return m
  
class SourceCatalog:
    def __init__(self, size, m_min, m_max, powerlaw, limits, seed):
      np.random.seed(seed)
      self.star_ID = np.arange(size).astype(int)
      self.alpha = np.random.uniform(low=limits[0],high=limits[1],size=size) # α
      self.beta = np.random.uniform(low=limits[2],high=limits[3],size=size) # β
      self.mag = generate_magnitudes(m_min,m_max,powerlaw,size)
      self.size = size
      self.flux = f.mag2flux(self.mag)
      
def create_catalog(M, m_min, m_max, powerlaw, limits, seed, plots=None, verbose=None):
  if verbose != None: print "Generating God's Catalog..."
  # Calculate total number of stars in catalog
  number_stars = M * (limits[1]-limits[0])*(limits[3]-limits[2])
  # Create catalog
  catalog = SourceCatalog(number_stars, m_min, m_max,powerlaw, limits, seed)
  if verbose != None: print "...done!"

  #**********************************************************
  #********************* Plot Catalog ***********************
  #**********************************************************
    
  if plots != None:
    fontsize = 20
    plt.figure(figsize=(14,6))
    st = plt.suptitle("Sky Catalog", fontsize=fontsize)
    # Plot portion of sky
    plt.subplot(121)
    if verbose != None: print "Plotting portion of sky and histogram..."
    plt.plot(catalog.alpha,catalog.beta,'.k',markersize=1)
    plt.xlabel(ur'$\alpha$', fontsize=fontsize)
    plt.ylabel(ur'$\beta$', fontsize=fontsize)
    #plt.title("God's Sky", fontsize=fontsize)
    plt.xlim(limits[0] - (limits[1]-limits[0])*0.1, limits[1] + (limits[1]-limits[0])*0.1)
    plt.ylim(limits[2] - (limits[3]-limits[2])*0.1, limits[3] + (limits[3]-limits[2])*0.1)
    ax = plt.gca()
    for tick in ax.xaxis.get_major_ticks():
      tick.label1.set_fontsize(fontsize)
    for tick in ax.yaxis.get_major_ticks():
      tick.label1.set_fontsize(fontsize)
    # Histogram of source magnitude
    plt.subplot(122)
    bin=np.arange(m_min,m_max,0.05)
    if verbose != None: print "Plotting histogram of source magnitude..."
    plt.hist(catalog.mag,bins=bin, log=True)
    #plt.title("%i Sources in Sample" % catalog.size, fontsize=fontsize)
    plt.xlabel("Source Magnitude", fontsize=fontsize)
    plt.ylabel("log(N)", fontsize=fontsize)
    ax = plt.gca()
    for tick in ax.xaxis.get_major_ticks():
      tick.label1.set_fontsize(fontsize)
    for tick in ax.yaxis.get_major_ticks():
      tick.label1.set_fontsize(fontsize)
      
    filename = "Figures/Catalog.png"
    plt.savefig(filename,bbox_inches='tight',pad_inches=0.5)
    if verbose != None: print "...done!"
  
  return catalog
  # catalog.ID, catalog.mag, catalog.alpha, catalog.beta, catalog.size