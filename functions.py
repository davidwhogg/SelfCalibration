#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pylab as plt

# Custom Modules
import parameters
import god

#*****************************************************
#************ Transformation Functions ***************
#*****************************************************
def mag2flux(mag):
  return np.power(10,(-0.4*(mag-22.5)))
  
def flux2mag(flux):
  return 22.5-2.5*np.log10(flux)

def sky2fp(alpha, beta, pointing, orientation): 
  theta = - orientation*np.pi/180 # telescope rotates NOT sky
  x = (alpha - pointing[0]) * np.cos(theta) - (beta - pointing[1]) * np.sin(theta)
  y = (alpha - pointing[0]) * np.sin(theta) + (beta - pointing[1]) * np.cos(theta)
  return x,y 

def fp2sky(x, y, pointing, orientation):
  theta = - orientation*np.pi/180 # telescope rotates NOT sky
  alpha = x * np.cos(theta) + y * np.sin(theta) + pointing[0]
  beta = -x * np.sin(theta) + y * np.cos(theta) + pointing[1]
  return alpha, beta


#*****************************************************
#************* Single Image Functions ****************
#*****************************************************

def single_image(catalog, pointing, orientation, plots=None, verbose=None):
  # catalog = [Star ID, magnitude, α, β], pointing = [α, β] 
  theta = orientation*np.pi/180
  
  if verbose != None: print "Converting sky catalog to focal plane coordinates..."
  camera_catalog = 0*catalog
  camera_catalog[:,0] = catalog[:,0]
  camera_catalog[:,1] = catalog[:,1]  
  camera_catalog[:,2], camera_catalog[:,3] = sky2fp(catalog[:,2], catalog[:,3], pointing, orientation)
  if verbose != None: print "...done!"  
  
  if verbose != None: print "Finding stars within camera FoV..."
  FoV = parameters.FoV()
  x_min = -FoV[0]/2; y_min = -FoV[1]/2
  x_max = FoV[0]/2; y_max = FoV[1]/2
  inside_FoV = np.where((x_min<camera_catalog[:,2]) & (camera_catalog[:,2]<x_max) & (y_min<camera_catalog[:,3]) & (camera_catalog[:,3]<y_max))
  stars_in_FoV = camera_catalog[inside_FoV[0]]
  if verbose != None: print "...done!"

  if verbose != None: print "Measuring stars within FoV..."
  # measured_sources = [star ID, observed_flux, observed_invvar, focal_position]
  measured_sources = np.zeros((len(inside_FoV[0]), 5))
  measured_sources[:,0] = stars_in_FoV[:,0]
  measured_sources[:,3] = stars_in_FoV[:,2]
  measured_sources[:,4] = stars_in_FoV[:,3]  
  # Calculate inverse variance 1/(σ*σ)
  measured_sources[:,2] = 1.0/flux_uncertainty_variance (flux2mag(mag2flux(stars_in_FoV[:,1])*god.flat_field(measured_sources[:,3],measured_sources[:,4])), parameters.eta())
  # Calculate observed flux
  measured_sources[:,1] = mag2flux(stars_in_FoV[:,1]) * god.flat_field(measured_sources[:,3],measured_sources[:,4]) + np.random.normal(size=len(measured_sources[:,0]))/np.sqrt(measured_sources[:,2])
  if verbose != None: print "...done!"
  
  #***************************************************
  
  # Plotting FoV on sky
  x = np.zeros(5)
  y = np.zeros(5)
  x[0] = x_min; y[0] = y_min
  x[1] = x_min; y[1] = y_max
  x[2] = x_max; y[2] = y_max
  x[3] = x_max; y[3] = y_min    
  x[4] = x_min; y[4] = y_min

  alpha, beta = fp2sky(x,y,pointing, orientation)


  #**********************************************************
  #********************* Plot Output ************************
  #**********************************************************
  
  if plots !=None:
    
    # Plot extracted sources over catalog
    plt.figure(figsize=(6,6))
    plt.plot(catalog[:,2],catalog[:,3],'o', markersize=2)
    temp = catalog[inside_FoV]
    plt.plot(temp[:,2],temp[:,3],'ro',markersize=2)
    plt.plot(alpha,beta,'k', linewidth=2)
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
    
    '''
    # Plot FoV on Catalog
    plt.figure(2001, figsize=(6,6))
    #plt.plot(catalog[:,2],catalog[:,3],'o', markersize=2)
    sky_limits = parameters.sky_limits()
    plt.plot(alpha,beta,'k', linewidth=2)
    plt.xlim(sky_limits[0]-FoV[0],sky_limits[1]+FoV[0])
    plt.ylim(sky_limits[2]-FoV[1],sky_limits[3]+FoV[1])
    plt.title('Fields-of-View on Sky', fontsize=20)
    plt.xlabel(ur'$\alpha$', fontsize=20)
    plt.ylabel(ur'$\beta$', fontsize=20)
    filename = './Figures/%s_FoV_on_sky.png' % plots
    #plt.show()
    plt.savefig((filename)) 
    '''

  return measured_sources
  # measured_sources = [star ID, observed_flux, observed_invvar, focal_position]


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


#*****************************************************
#**************** Survey Functions *******************
#*****************************************************

def survey(catalog, survey_file, plots=None, verbose=None):
  # Routing loads the survey from file and then calls the camera function accordingly.
  
  if verbose != None: print "Loading survey..."
  # pointing = [pointing ID, α, β, orientation]
  pointing = np.loadtxt(survey_file)
  # Calculate number of pointings
  number_pointings = len(pointing[:,0])  
  if verbose != None: print "...done!"
  
  # Declare dictionary
  data_dic = {}
  
  # Perform survey
  if verbose != None: print "Surveying sky..."
  for i in range(0, number_pointings):
    data_dic[i] = single_image(catalog, [pointing[i,1],pointing[i,2]], pointing[i,3], plots=plots, verbose = verbose)
    # data_dic[pointing ID] = [star ID, observed_flux, observed_invvar, focal_position]
  if verbose != None: print "...done!"
  
  # Rearrange observations into an observation_catalog
  if verbose != None: print "Rearranging observation catalog"
  array_size = 0
  for i in range(0,len(data_dic)):
    array_size = array_size + len(data_dic[i][:,0])
  
  if verbose != None: print "Observation catalog requires %d rows)" % array_size
  
  # Copying dictionary into array
  observation_catalog = np.zeros((array_size,6))
  count = 0
  for i in range(0,len(data_dic)):
    single_exposure = data_dic[i]
    for j in range(0,len(single_exposure[:,0])):
      observation_catalog[count,0] = i
      observation_catalog[count,1:6] = single_exposure[j,:]
      count = count +1
      
  return observation_catalog
  # observed_catalog = [pointing ID, star ID, observed_flux, observed_invvar, focal plane x, focal plane y]


#*****************************************************
#**************** Ubercal Functions ******************
#*****************************************************

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
  chi2 = np.sum((np.dot(g,q) * ss - obs_cat[:,2])**2 * obs_cat[:,3])
  return (q, q_invvar, chi2)
  

#*****************************************************
#*************** Analysis Functions ******************
#*****************************************************

def rms_error(flux_estimate, true_flux):
  return 100*np.sqrt(np.mean((flux_estimate-true_flux)/true_flux)**2)
  
#*****************************************************
#*************** Plotting Functions ******************
#*****************************************************

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
  a = plt.imshow((100*(our_ff-god_ff)/god_ff),extent=(-FoV[0]/2,FoV[0]/2,-FoV[1]/2,FoV[1]/2), vmin = -1,vmax = 1, cmap='gray')
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