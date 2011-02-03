#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import random
import tparameters as p
import functions
import matplotlib.pylab as plt

pdic = p.dic
FoV = pdic['FoV']
sky_limits = pdic['sky_limits']

def plot_survey(survey, survey_name):
  x_min = -FoV[0]/2; y_min = -FoV[1]/2
  x_max = FoV[0]/2; y_max = FoV[1]/2
  x = np.array([x_min, x_min, x_max, x_max, x_min])
  y = np.array([y_min, y_max, y_max, y_min, y_min])
  for image in survey:
    alpha, beta = functions.fp2sky(x,y,image[1:3], image[3])
    plt.plot(alpha,beta,'k-',alpha=0.25)
  plt.xlabel(r"$\alpha$", fontsize=25)
  plt.ylabel(r"$\beta$", fontsize=25)
  plt.title(r"Strategy %s: %i Pointings" % (survey_name, len(survey[:,0])))
  plt.xlim(sky_limits[0]-FoV[0], sky_limits[1]+FoV[0])
  plt.axis('equal')
  plt.ylim(sky_limits[2]-FoV[1], sky_limits[3]+FoV[1])
  return

def generate_uniform_survey(Ncovering):
  nx = np.ceil((sky_limits[1]-sky_limits[0])/(0.975*FoV[0])).astype(int)
  ny = np.ceil((sky_limits[3]-sky_limits[2])/(0.975*FoV[1])).astype(int)
  x_center_list = sky_limits[0] + 0.5*FoV[0] + (sky_limits[1]-sky_limits[0] - FoV[0])*np.arange(nx)/(nx-1.)
  y_center_list = sky_limits[2] + 0.5*FoV[1] + (sky_limits[3]-sky_limits[2] - FoV[1])*np.arange(ny)/(ny-1.)
  x = np.zeros((Ncovering*nx*ny,4))
  ii = 0
  for covering in range(Ncovering):
    for yy in range(ny):
      x[ii:ii+nx,0] = ii+np.arange(nx)
      x[ii:ii+nx,1] = x_center_list
      x[ii:ii+nx,2] = y_center_list[yy]
      x[ii:ii+nx,3] = 0.
      ii += nx

  return x
'''
def generate_dithered_survey():
  index = 0
  countx = 0.5
  county = 0.5
  for i in range(0,no_pointings):
    if index == 0:  
      x[i,1] = countx 
      x[i,2] = county 
      index = 1
    elif index == 1:
      x[i,1] = countx + 0.05
      x[i,2] = county 
      index = 2
    elif index == 2:
      x[i,1] = countx - 0.05
      x[i,2] = county + 0.05
      index = 3  
    elif index == 3:
      x[i,1] = countx + 0.05
      x[i,2] = county
      countx = countx + 0.7 - 0.05 
      if countx + 0.7 >= 9.5:
        countx = 0.7
        county=county + 0.7 - 2*0.05
      index = 0
    x[i,0] = i
    x[i,3] = 0
  return x
'''
def generate_random_survey(N):
  x = np.zeros((N,4))
  x[:,0] = range(N)
  x[:,1] = np.random.uniform(sky_limits[0]+FoV[0]/2,sky_limits[1]-FoV[0]/2,size=(N))
  x[:,2] = np.random.uniform(sky_limits[2]+FoV[1]/2,sky_limits[3]-FoV[1]/2,size=(N))
  x[:,3] = np.random.uniform(0,360, size =(N))
  return x

if __name__ == "__main__":
  plt.figure(figsize=(15,6))
  plt.subplot(121)
  xA = generate_uniform_survey(9)
  plot_survey(xA,'A')
  np.savetxt('A.txt',xA)
  xD = generate_random_survey(len(xA))
  plt.subplot(122)
  plot_survey(xD,'D')
  np.savetxt('D.txt',xD)
  plt.savefig("./Plotting_Data/surveys.png",bbox_inches='tight',pad_inches=0.)
  