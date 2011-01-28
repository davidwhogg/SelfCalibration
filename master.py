#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pylab as plt
import math
import os
# Set up LaTeX for plots
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)

# Custom Modules
import functions
import god
import parameters
import camera
import survey

# Remove all old plots
os.system('rm ./Figures/Camera_Images/*.png')
os.system('rm ./Figures/*.png')
os.system('rm ./Figures/Flat_Fields/*.png')

plots = None
verbose = None

#*************************************************************
#******************* Generate Sky Catalog ********************
#*************************************************************

M = parameters.density_of_stars() # density of stars (all magnitudes)
m_min = parameters.m_min() # minimum magnitude (i.e. saturation limit)
m_max = parameters.m_max() # maximum magnitude limit per dither (i.e. sensitivity limit 10σ)
powerlaw = parameters.powerlaw() #B in log10(dN/dm) = A + B*m
sky = 'sq' # create a catalog on a square (='sq') or sphere (='sp')
seed = 1 # seed for random number generation, providing same seed generate the same sky
limits = parameters.sky_limits() # limits of spacial/angular coordinates [α_min, α_max, β_min, β_max]

# Creates catalog of stars used for calibration
sky_catalog = god.create_catalog(M, m_min, m_max, powerlaw, sky, limits, seed, plots=plots, verbose = verbose)
# sky_catalog = [Star ID, magnitude, α, β]


#*************************************************************
#************************ Survey Sky *************************
#*************************************************************

observation_catalog = survey.survey(sky_catalog, "A.txt", plots=plots, verbose=verbose) 
# observed_catalog = [pointing ID, star ID, observed_flux, observed_invvar, focal plane x, focal plane y]

#*************************************************************
#********************* Ubercalibration ***********************
#*************************************************************

q = np.array([1,0,0,0,0,0])
order = 2
for iteration_number in range(0,50):
  s, s_invvar = functions.s_step(observation_catalog,q)
  q, q_invvar = functions.q_step(observation_catalog, s, order,iteration_number,plots=plots)

#*************************************************************
#*********************** Health Checks ***********************
#*************************************************************

if plots != None:
  # Flux Uncertainty variance check
  if verbose != None: print "Plotting flux uncertainty variance..."
  eta = parameters.eta()
  plt.figure(1002)
  # Calculate uncertainty for different magnitudes
  temp_mag = np.arange(m_min-0.5,m_max+0.5,0.1)
  temp_uncert = np.log10( (functions.flux_uncertainty_variance(temp_mag,eta))/(functions.mag2flux(temp_mag)**2))
  plt.plot(temp_mag,temp_uncert)
  plt.ylim(np.min(temp_uncert)-0.5,np.max(temp_uncert)+0.5)
  # Draw vertical (dashed lines) at the magnitude limits take 
  lower_limit_mag = [m_min, m_min]
  lower_limit_uncer = [np.min(temp_uncert),np.max(temp_uncert)]
  plt.plot(lower_limit_mag,lower_limit_uncer,'k--')
  upper_limit_mag = [m_max, m_max]
  upper_limit_uncer = [np.min(temp_uncert),np.max(temp_uncert)]
  plt.plot(upper_limit_mag,upper_limit_uncer,'k--')
  plt.xlabel(ur'Source Magnitude', fontsize=16)
  plt.ylabel(ur'$\log_{10}(\frac{{\sigma_f}^2}{f^2})$',fontsize=25)#, rotation='horizontal')
  # Label vertical lines
  plt.annotate(r"$m_{min}$", (m_min-0.1,np.max(temp_uncert)),fontsize=16) 
  plt.annotate(r"$m_{max}$", (m_max-0.1,np.max(temp_uncert)),fontsize=16) 
  plt.savefig('Figures/flux_uncertainty_variance.png',bbox_inches='tight',pad_inches=0.1)
  if verbose != None: print "...done!"

# Print sum of random numbers to check random seed is working (i.e. )
if verbose != None: print "Sum of all numbers in sky_catalog = %0.10lf" % np.sum(sky_catalog)