#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import random
import parameters

FoV = parameters.FoV()
sky_limits = parameters.sky_limits()



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
  xA = generate_uniform_survey(9)
  np.savetxt('A.txt',xA)
  xD = generate_random_survey(len(xA))
  np.savetxt('D.txt',xD)
