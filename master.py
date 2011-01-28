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
#import camera
#import survey

# Remove all old plots
os.system('rm ./Figures/Camera_Images/*.png')
os.system('rm ./Figures/*.png')
os.system('rm ./Figures/Flat_Fields/*.png')
os.system('rm ./Figures/Flat_Fields/*.gif')

if __name__ == "__main__":
  for strategy in ['D', 'A']: #['A', 'D']:

    plots = None
    ff_plots = None
    verbose = None

    #********************************************************
    #*************** Generate Sky Catalog *******************
    #********************************************************

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


    #********************************************************
    #******************** Survey Sky ************************
    #********************************************************
    survey_file = strategy + ".txt"
    observation_catalog = functions.survey(sky_catalog, survey_file, plots=plots, verbose=verbose) 
    # observed_catalog = [pointing ID, star ID, observed_flux, observed_invvar, focal plane x, focal plane y]

    #********************************************************
    #********************* Ubercalibration ******************
    #********************************************************

    order = 2
    if strategy == 'A': 
      niter = 1000
    else:
      niter = 15
      q = np.array([1,0,0,0,0,0])
    
    for iteration_number in range(0,niter):
      # Print flat-field parameters
      
      s, s_invvar = functions.s_step(observation_catalog,q)
      q, q_invvar, chi2 = functions.q_step(observation_catalog, s, order,iteration_number,plots=plots)
      
      # Calculate rms error in stars
      indx = [s != 0]
      obs_star_actual = sky_catalog[indx]
      rms = functions.rms_error(s[indx], functions.mag2flux(obs_star_actual[:,1]))
      
      print "%i: f(x,y) =%.2f%s%.2fx%s%.2fy%s%.2fxx%s%.2fxy%s%.2fyy   RMS = %.6f %%   chi2 = %0.2f (%i)" % (iteration_number+1, abs(q[0]),  functions.sign(q[1]), abs(q[1]), functions.sign(q[2]), abs(q[2]), functions.sign(q[3]), abs(q[3]), functions.sign(q[4]), abs(q[4]), functions.sign(q[5]), abs(q[5]), rms, chi2, len(observation_catalog[:,0]))
      
      #if plots != None: 
      if (ff_plots == 'y') or (iteration_number+1==niter): functions.plot_flat_fields(q, (iteration_number+1), plots=strategy)
    
    
    
      # Animate flat-field fitting
    if ff_plots == 'y': 
      os.system(("convert -delay 20 -loop 0 ./Figures/Flat_Fields/%s*.png ./Figures/Flat_Fields/%s_00_animation.gif" % (strategy,strategy)))
      
    #********************************************************
    #*********************** Health Checks ******************
    #********************************************************

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