#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import time
import numpy as np
import matplotlib.pylab as plt

os.system('rm ./test.txt')

parameter = 'epsilon_max'
par_range = np.arange(0,1,0.1)

for i in range(len(par_range)):
  command = './master.py %s %s' % (parameter, par_range[i])
  print command
  os.system(command)
  time.sleep(0)
  i+=1

data = np.loadtxt('test.txt')

plt.plot(data[:,0], data[:,1])
plt.xlabel('Epsilon')
plt.ylabel('Badness')
plt.show()