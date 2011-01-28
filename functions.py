#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pylab as plt

import parameters
import god

def mag2flux(mag):
  return np.power(10,(-0.4*(mag-22.5)))
  
def flux2mag(flux):
  return 22.5-2.5*np.log10(flux)
  
def flux_uncertainty_variance(mag, eta):
  sky_error= 0.1*mag2flux(parameters.mag_at_ten_sigma())
  return (np.power(sky_error,2) + (eta**2)*np.power(mag2flux(mag),2))

def evaluate_flat_field_functions(x, y, order):
  L = (order + 1)*(order + 2)/2
  g = np.zeros((len(x),L))
  l = 0
  for n in range(order+1):
    for m in range(n+1):
      g[:,l] = (x**(n-m)) * (y**m)
      l += 1 
  return g

def evaluate_flat_field(x, y, q):
  # Function returns flat field at focal plane position specified. The function also calculates the order of the flat field equation from the number of flat field parameters
  
  # Calculate required order
  order = int(np.around(np.sqrt(0.25+2*len(q))-1.5))
  assert(len(q) == ((order + 1) * (order + 2) / 2)) # assert to break code if condition not met
  g = evaluate_flat_field_functions(x, y, order)
  return np.dot(g, q)
 
def normalize_flat_field(q):
  return (q / q[0])

def s_step(obs_cat, q):
  # First step in ubercalibration. This function calculates the new flux estimates.
  # obs_cat = observation_catalog = [pointing ID, star ID, observed_flux, observed_invvar, focal plane x, focal plane y]
  
  # Calculate flat-field based on source positions and current flat field parameters 
  
  ff = evaluate_flat_field(obs_cat[:,4], obs_cat[:,5], q)
  fcss = ff * obs_cat[:,2] * obs_cat[:,3]  
  ffss = ff * ff * obs_cat[:,3]
  star_ID = np.around(obs_cat[:,1]).astype(int)
  max_star = np.max(star_ID)
  s = np.zeros(max_star+1)
  s_invvar = np.zeros(max_star)
  for ID in range(max_star):
    indx = (star_ID == ID)
    denominator = np.sum(ffss[indx])
    s_invvar[ID] = denominator
    if denominator > 0.:
      s[ID] = np.sum(fcss[indx])/denominator
  return (s, s_invvar)
  
def q_step(obs_cat, s, order, iteration_number, plots=None):
  g = evaluate_flat_field_functions(obs_cat[:,4],obs_cat[:,5], order)
  ss = s[np.around(obs_cat[:,1]).astype(int)]
  q_invvar = np.dot(np.transpose(g)*ss*ss*obs_cat[:,3],g)
  numerator = np.sum(np.transpose(g)*ss*obs_cat[:,2]*obs_cat[:,3],axis=1)
  q = np.dot(np.linalg.inv(q_invvar),numerator)
  q = normalize_flat_field(q)
  np.set_printoptions(precision=2) 
  return (q, q_invvar)
  
def plot_flat_fields(our_q, iteration_number,plots=None):
  FoV = parameters.FoV() 
  
  plt.figure(3000,figsize=(13, 6)) # want square figure
  plt.subplot(121)

  plt.title(r"Flat-Fields (God's = Black; Fitted = Red) Iteration: %i" % iteration_number)
  plt.xlabel(r"\alpha")
  plt.ylabel(r"$\beta$")
  x = np.arange(-FoV[0]/2,FoV[0]/2,0.01)
  y = np.arange(-FoV[1]/2,FoV[1]/2,0.01)
  X, Y = np.meshgrid(x, y)
  
  # Have to reshape so that evaluate_flat_field() works 
  temp_x = np.reshape(X,-1)
  temp_y = np.reshape(Y,-1)
  temp_z = evaluate_flat_field(temp_x,temp_y,our_q)
  our_ff = np.reshape(temp_z, (len(X[:,0]),len(X[0,:])))
  god_q = god.flat_field_parameters()
  temp_z = evaluate_flat_field(temp_x,temp_y,god_q)
  god_ff = np.reshape(temp_z, (len(X[:,0]),len(X[0,:])))

  # Find parameters for contour plot
  god_ff_max = np.max(god_ff)
  god_ff_min = np.min(god_ff)
  step = (god_ff_max-god_ff_min)/10
  
  
  CS = plt.contour(X,Y,god_ff,np.arange(god_ff_min,god_ff_max,step),colors='k')
  plt.clabel(CS, fontsize=9, inline=1)
  CS2 = plt.contour(X,Y,our_ff,np.arange(god_ff_min,god_ff_max,step),colors='r',alpha=0.5)

  # Write formulas on plot
  our_formula = "$f(x,y) =%.2f%s%.2fx%s%.2fy%s%.2fx^2%s%.2fxy%s%.2fy^2$" % (abs(our_q[0]),  sign(our_q[1]), abs(our_q[1]), sign(our_q[2]), abs(our_q[2]), sign(our_q[3]), abs(our_q[3]), sign(our_q[4]), abs(our_q[4]), sign(our_q[5]), abs(our_q[5]) )
  plt.text(-.38,-.35,our_formula, color='r',bbox = dict(boxstyle="square",ec='w',fc='w', alpha=0.9), fontsize=9.5)
  god_formula = "$f(x,y) =%.2f%s%.2fx%s%.2fy%s%.2fx^2%s%.2fxy%s%.2fy^2$" % (abs(god_q[0]),  sign(god_q[1]), abs(god_q[1]), sign(god_q[2]), abs(god_q[2]), sign(god_q[3]), abs(god_q[3]), sign(god_q[4]), abs(god_q[4]), sign(god_q[5]), abs(god_q[5]) )
  plt.text(-.38,-.39,god_formula, color='k',bbox = dict(boxstyle="square",ec='w',fc='w', alpha=0.9), fontsize=9.5)

  '''
  from mpl_toolkits.mplot3d import Axes3D
  fig = plt.figure()
  ax = Axes3D(fig)
  ax.plot_wireframe(X, Y, god_ff,colors='k')
  ax.plot_wireframe(X, Y, our_ff,colors='r',alpha=0.0)
  ax.set_zlim3d(god_ff_min, god_ff_max)
  '''
  
  
  # Plot residual in flat-field
  plt.subplot(122)
  plt.title(r"Residual Error in Fitted Flat-Field (\%)")
  a = plt.imshow((100*(our_ff-god_ff)/god_ff),extent=(-FoV[0]/2,FoV[0]/2,-FoV[1]/2,FoV[1]/2), vmin = 0)#,vmax = 5)
  plt.colorbar(a,shrink=0.7)
  plt.xlabel(r"\alpha")
  plt.ylabel(r"$\beta$")
  
  # Save figure
  if iteration_number < 10:
    filename = 'Figures/Flat_Fields/%s_0%d_ff.png' % (plots,iteration_number)
  else:
    filename = 'Figures/Flat_Fields/%s_%d_ff.png' % (plots, iteration_number)
  
  plt.savefig(filename)
  plt.clf()  

def sign(x):
  if x >= 0:
    x = '+'
  else:
    x = '-'
  return x