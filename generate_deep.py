#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Rory Holmes 2011
# script to generate deeep field survey strategy for ./cross_cal.py simulation 
# Requires: default_parameters.py

import numpy as np
import default_parameters as p
import matplotlib.pylab as plt

scale = 2
fig_width_pt = scale*415.55  # Get this from LaTeX using \showthe\columnwidth
inches_per_pt = 1.0/72.27               # Convert pt to inches
golden_mean = (np.sqrt(5)-1.0)/2.0         # Aesthetic ratio
fig_width = fig_width_pt*inches_per_pt  # width in inches
fig_height =fig_width*golden_mean       # height in inches
params = {'backend': 'pdf',
          'axes.labelsize': scale*10,
          'text.fontsize': scale*10,
          'legend.fontsize': scale*10,
          'xtick.labelsize': scale*9,
          'ytick.labelsize': scale*9,
          'linewidth' : scale*1,
          'text.usetex': True}
          #'font': {'family':'sans-serif','sans-serif':['Helvetica']}}
plt.rcParams.update(params)
plt.rc('font', family = 'serif', serif = 'cmr10')


# XX Dangerous hack!!
pdic = p.dic
FoV = pdic['FoV']
sky_limits = pdic['sky_limits']
field_overlap = 0.02
dither_overlap = 0.98
no_dithers = 4
square_num = 5.
passes = 12 # per year
angle_between_passes = 30

def single_pass(offset, theta):
  # Create Center Square
  x_low = -FoV[0]*(square_num+1)*0.5*(1-field_overlap)
  x_high = -x_low
  y_low = -FoV[1]*(square_num+1)*0.5*(1-field_overlap)
  y_high = -y_low
  x_centers = np.linspace(x_low,x_high, num = (square_num+2), endpoint=True)
  y_centers = np.linspace(y_low,y_high, num = (square_num+2), endpoint=True)
  X_centers,Y_centers = np.meshgrid(x_centers[1:-1],y_centers[1:-1])
  # Add Edges
  for ii in range(int(square_num-1)):
    X_centers = np.append(X_centers, x_centers[0])
    Y_centers = np.append(Y_centers, y_centers[0]+(ii+1.5) * FoV[1])
    X_centers = np.append(X_centers, x_centers[-1])
    Y_centers = np.append(Y_centers, y_centers[0]+(ii+1.5) * FoV[1])
    X_centers = np.append(X_centers, x_centers[0]+(ii+1.5) * FoV[0])
    Y_centers = np.append(Y_centers, y_centers[0])
    X_centers = np.append(X_centers, x_centers[0]+(ii+1.5) * FoV[0])
    Y_centers = np.append(Y_centers, y_centers[-1])
  # Create Dithers
  X = np.zeros((no_dithers*len(X_centers), 2))
  for ii in range(len(X_centers)):
    dx = 0.5*(1-dither_overlap) * FoV[0]
    dy = 0.5*(1-dither_overlap) * FoV[1]
    add_sub = np.array([[-1,-1,1,1], [-1,1,-1,1]])
    for jj in range(len(add_sub[0,:])):
      X[4*ii+jj,0] = X_centers[ii] + add_sub[0,jj] * dx
      X[4*ii+jj,1] = Y_centers[ii] + add_sub[1,jj] * dy
  # Rotate
  X_rot = np.zeros((no_dithers*len(X_centers), 3))
  X_rot[:,0] = X[:,0]*np.cos(theta*np.pi/180) - X[:,1]*np.sin(theta*np.pi/180) + offset[0]
  X_rot[:,1] = X[:,0]*np.sin(theta*np.pi/180) + X[:,1]*np.cos(theta*np.pi/180) + offset[1]
  X_rot[:,2] = theta
  return X_rot #[x, y, theta]

X = np.zeros((0,3))
for ii in range(passes):
  offset = np.random.uniform(FoV[0]*0.25, FoV[0]*0.75, size = 2)
  temp_X = single_pass(offset, ii*angle_between_passes)
  X = np.append(X, temp_X, axis = 0)

final_X = np.zeros((len(X[:,0]),4))

final_X[:,0] = np.arange(len(X[:,0]))
final_X[:,1:] = X
print len(final_X[:,0])/12
first_pass = final_X[0:len(final_X[:,0])/12,:]

def fp2sky(x, y, pointing, orientation):
  theta = - orientation*np.pi/180 # telescope rotates NOT sky
  alpha = x * np.cos(theta) + y * np.sin(theta) + pointing[0]
  beta = -x * np.sin(theta) + y * np.cos(theta) + pointing[1]
  return alpha, beta

def plot_survey(survey, survey_name):
  x_min = -FoV[0]/2; y_min = -FoV[1]/2
  x_max = FoV[0]/2; y_max = FoV[1]/2
  x = np.array([x_min, x_min, x_max, x_max, x_min])
  y = np.array([y_min, y_max, y_max, y_min, y_min])
  
  if survey_name == "Deep": Alpha = 0.1
  if survey_name == "Single Pass": Alpha = 0.3
  for image in survey:
    alpha, beta = fp2sky(x,y,image[1:3], image[3])
    plt.plot(alpha,beta,'k-',alpha=Alpha)

  plt.axis('scaled')
  #plt.xlim(sky_limits[0]-FoV[0], sky_limits[1]+FoV[0])
  #plt.ylim(sky_limits[2]-FoV[1], sky_limits[3]+FoV[1])
  plt.xlim(sky_limits[0]-1*FoV[0], sky_limits[1]+1*FoV[0])
  plt.ylim(sky_limits[2]-1*FoV[1], sky_limits[3]+1*FoV[1]) 
  print (r"Strategy %s: %i Pointings" % (survey_name, len(survey[:,0])))
  return

alpha = r'Sky Position $\alpha$ (deg$^2$)'
beta = r'Sky Position $\beta$ (deg$^2$)'

plt.figure(figsize = (fig_width, fig_height))
plt.subplot(121)
plot_survey(first_pass, "Single Pass")
plt.xlabel(alpha)
plt.ylabel(beta)

plt.subplot(122)
plot_survey(final_X, "Deep")
plt.gca().set_yticklabels([])
plt.xlabel(alpha)

np.savetxt("deep.txt", final_X)
plt.subplots_adjust(wspace=0,hspace=0.0)

plt.savefig("deep.pdf",bbox_inches='tight')
