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
import functions as f
import god
import parameters as p
#import camera
#import survey

# Remove all old plots
os.system('rm ./Figures/Camera_Images/*.png')
os.system('rm ./Figures/*.png')
os.system('rm ./Figures/*.pdf')
os.system('rm ./Figures/Flat_Fields/*.png')
os.system('rm ./Figures/Flat_Fields/*.gif')

if __name__ == "__main__":
  for strategy in ['D']: #['A', 'D']:

    catalog_plots = None
    survey_plots = None # strategy
    coverage_plots = None # strategy
    ff_plots = None # None, 'all'
    verbose = None
    health_plots = None

    #********************************************************
    #*************** Generate Sky Catalog *******************
    #********************************************************
    # Creates catalog of stars used for calibration
    sky_catalog = god.create_catalog(p.density_of_stars(), p.m_min(), p.m_max(), p.powerlaw(), p.sky_limits(), seed = 1, plots = catalog_plots, verbose = verbose)
        # powerlaw = B in log10(dN/dm) = A + B*m
        # sky_limits = [α_min, α_max, β_min, β_max]
        # returns sky_catalog = [Star ID, magnitude, α, β]
        # sky_catalog: *.star_ID, *.mag, *.alpha, *.beta, *.size

    
    #********************************************************
    #******************** Survey Sky ************************
    #********************************************************
    survey_file = strategy + ".txt"
    observation_catalog = f.survey(sky_catalog, survey_file, plots=survey_plots, verbose=verbose) 
        # observation_catalog: *.size, *.pointing_ID, *.star_ID, *.flux, *.invvar, *.x, *.y

    if coverage_plots != None: f.coverage(observation_catalog, strategy)

    #********************************************************
    #********************* Ubercalibration ******************
    #********************************************************
    f.ubercalibration(observation_catalog,sky_catalog,strategy,ff_plots=ff_plots)
    
    #********************************************************
    #*********************** Health Checks ******************
    #********************************************************

    if health_plots != None:
      # Flux Uncertainty variance check
      if verbose != None: print "Plotting flux uncertainty variance..."
      eta = p.eta()
      plt.figure(1002)
      # Calculate uncertainty for different magnitudes
      temp_mag = np.arange(m_min-0.5,m_max+0.5,0.1)
      temp_uncert = np.log10( (f.flux_uncertainty_variance(temp_mag,eta))/(f.mag2flux(temp_mag)**2))
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
