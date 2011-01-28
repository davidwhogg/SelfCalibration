#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import math
import matplotlib.pylab as plt

import parameters
import functions
import god

def camera(catalog, pointing, orientation, plots=None, verbose=None):
  # Measure the sources within the camera's FoV
  # catalog = [Star ID, magnitude, α, β]
  # pointing = [α, β] # telescope's pointing "sky" coordinates
  # orientation # telescope's orientation
  # if plots == None then do not make any plots, if plots == 'sdadfafdg' then make this the prefix on the plots
  
  #**********************************************************
  #***** Transform sky into camera coordinate systems *******
  #**********************************************************  

  # INEFFICIENT - transforming entire sky catalog into focal plane coordinates
  # Create a camera catalog [Star ID, magnitude, x, y]
  if verbose != None: print "Converting sky catalog to focal plane coordinates..."
  camera_catalog = 0*catalog
  
  # Covert orientation from degrees to radians
  theta = orientation*math.pi/180

  # Transform sky coordinates into focal plane coordinates
  camera_catalog[:,2] = (catalog[:,2] - pointing[0]) * math.cos(theta) - (catalog[:,3]-pointing[1]) * math.sin(theta)
  camera_catalog[:,3] = (catalog[:,2] - pointing[0]) * math.sin(theta) + (catalog[:,3] - pointing[1]) * math.cos(theta)

  # Leave star ID and magnitude unchanged
  camera_catalog[:,0] = catalog[:,0]
  camera_catalog[:,1] = catalog[:,1]
  
  if verbose != None: print "...done!"  
  
  #**********************************************************
  #************* Find stars within camera's FoV *************
  #**********************************************************
  
  if verbose != None: print "Finding stars within camera FoV..."
  # Define camera field-of-view
  FoV = parameters.FoV()
  x_min = -FoV[0]/2
  y_min = -FoV[1]/2
  x_max = FoV[0]/2
  y_max = FoV[1]/2

  # Find all stars within field-of-view
  inside_FoV = np.where((x_min<camera_catalog[:,2]) & (camera_catalog[:,2]<x_max) & (y_min<camera_catalog[:,3]) & (camera_catalog[:,3]<y_max))

  # Extract only stars within field-of-view
  stars_in_FoV = camera_catalog[inside_FoV[0]]
  if verbose != None: print "...done!"

  #**********************************************************
  #****************** "Measure" Sources *********************
  #**********************************************************
 
  if verbose != None: print "Measuring Stars..."
 
  # measured_sources = [star ID, observed_flux, observed_invvar, focal_position]
  measured_sources = np.zeros((len(inside_FoV[0]), 5))

  # Copy star ID
  measured_sources[:,0] = stars_in_FoV[:,0]

  # Copy focal plane position
  measured_sources[:,3] = stars_in_FoV[:,2]
  measured_sources[:,4] = stars_in_FoV[:,3]
  
  # Calculate inverse variance 1/(σ*σ)
  measured_sources[:,2] = 1.0/functions.flux_uncertainty_variance (functions.flux2mag(functions.mag2flux(stars_in_FoV[:,1])*god.flat_field(measured_sources[:,3],measured_sources[:,4])), parameters.eta())
  
  # Calculate observed_flux
  measured_sources[:,1] = functions.mag2flux(stars_in_FoV[:,1]) * god.flat_field(measured_sources[:,3],measured_sources[:,4]) + np.random.normal(size=len(measured_sources[:,0]))/np.sqrt(measured_sources[:,2])
  
  if verbose != None: print "...done!"

  #**********************************************************
  #********** Translate FoV into sky coordinates ************
  #**********************************************************

  # Plotting FoV on sky
  x = np.zeros(5)
  y = np.zeros(5)

  x[0] = x_min*math.cos(theta)+y_min*math.sin(theta)+pointing[0] 
  y[0] = -x_min*math.sin(theta)+y_min*math.cos(theta)+pointing[1] 

  x[1] = x_min*math.cos(theta)+y_max*math.sin(theta)+pointing[0] 
  y[1] = -x_min*math.sin(theta)+y_max*math.cos(theta)+pointing[1] 

  x[3] = x_max*math.cos(theta)+y_min*math.sin(theta)+pointing[0] 
  y[3] = -x_max*math.sin(theta)+y_min*math.cos(theta)+pointing[1] 

  x[2] = x_max*math.cos(theta)+y_max*math.sin(theta)+pointing[0] 
  y[2] = -x_max*math.sin(theta)+y_max*math.cos(theta)+pointing[1] 

  x[4] = x_min*math.cos(theta)+y_min*math.sin(theta)+pointing[0] 
  y[4] = -x_min*math.sin(theta)+y_min*math.cos(theta)+pointing[1]


  #**********************************************************
  #********************* Plot Output ************************
  #**********************************************************
  
  if plots != None:
    
    # Plot extracted sources over catalog
    plt.figure(figsize=(6,6))
    plt.plot(catalog[:,2],catalog[:,3],'o', markersize=2)
    temp = catalog[inside_FoV]
    plt.plot(temp[:,2],temp[:,3],'ro',markersize=2)
    plt.plot(x,y,'k', linewidth=2)
    title = ur'$\alpha$ = %.1lf, $\beta$ = %.1lf, $\theta$ = %.1lf$^\circ$' % (pointing[0],pointing[1], orientation)
    plt.title(title, fontsize=20)
    plt.xlabel(ur'$\alpha$', fontsize=20)
    plt.ylabel(ur'$\beta$', fontsize=20)
    filename = "Figures/Camera_Images/%s_full_sky_alpha_%.1lf_beta_%.1lf_rot_%.1lf.png" % (plots, pointing[0],pointing[1], orientation)
    print filename
    plt.savefig(filename)

    # Plot sources on focal plane
    plt.figure(figsize=(6,6))
    plt.plot(measured_sources[:,3],measured_sources[:,4],'o', markersize=2)
    title = ur'$\alpha$ = %.1lf, $\beta$ = %.1lf, $\theta$ = %.1lf$^\circ$' % (pointing[0],pointing[1], orientation)
    plt.title(title, fontsize=20)
    plt.xlabel(ur'Focal Plane $x$', fontsize=20)
    plt.ylabel(ur'Focal Plane $y$', fontsize=20)
    filename = "Figures/Camera_Images/%s_fp_alpha_%.1lf_beta_%.1lf_rot_%.1lf.png" % (plots, pointing[0],pointing[1], orientation)
    print filename
    plt.savefig(filename)
    
    # Plot FoV on Catalog
    plt.figure(2001, figsize=(6,6))
    #plt.plot(catalog[:,2],catalog[:,3],'o', markersize=2)
    sky_limits = parameters.sky_limits()
    plt.plot(x,y,'k', linewidth=2)
    plt.xlim(sky_limits[0]-FoV[0],sky_limits[1]+FoV[0])
    plt.ylim(sky_limits[2]-FoV[1],sky_limits[3]+FoV[1])
    plt.title('Fields-of-View on Sky', fontsize=20)
    plt.xlabel(ur'$\alpha$', fontsize=20)
    plt.ylabel(ur'$\beta$', fontsize=20)
    filename = './Figures/%s_FoV_on_sky.png' % plots
    plt.savefig((filename)) 

  return measured_sources
  # measured_sources = [star ID, observed_flux, observed_invvar, focal_position]
