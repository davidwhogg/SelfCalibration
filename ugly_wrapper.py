#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import time
import numpy as np
import matplotlib.pylab as plt

os.system('rm ./test.txt')

parameter = 'epsilon_max'
par_range = np.arange(0,1,0.5)

for i in par_range:
  command = './master.py %s %lf' % (parameter, par_range[i])
  print command
  os.system(command)
  time.sleep(0)

data = np.loadtxt('test.txt')

plt.plot(data[:,0], data[:,1])
plt.show()