#!/usr/bin/env python
# -*- coding: utf-8 -*-


import numpy as np
import random

survey = 2

no_pointings = 100
x = np.zeros((no_pointings,4))

countx = 0
county = 0

if survey == 1:
  # Generate Uniform Survey
  for i in range(0,36):
    x[i,0] = i
    x[i,1] = countx + 0.50
    x[i,2] = county + 0.50
    x[i,3] = random.uniform(0,360)
    countx = countx +0.35
    if countx >= 9.5:
      countx = 0
      county=county+0.35

elif survey == 2:
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

if survey == 3:
  # Generate Uniform Survey
  for i in range(0,36):
    x[i,0] = i
    x[i,1] = random.uniform(0.5,9.5)
    x[i,2] = random.uniform(0.5,9.5)
    x[i,3] = random.uniform(0,360)

#print x

np.savetxt('A.txt',x)
