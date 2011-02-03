#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import time
import numpy as np
import matplotlib.pylab as plt
import glob
import pickle

os.system('mkdir -p ./Plotting_Data/Wrapped/')

def clean_up_old_files(parameter):
  existing_file = "./Plotting_Data/Wrapped/*%s.txt_old" % parameter
  os.system(('rm -f '+existing_file))
  os.chdir("./Plotting_Data/Wrapped")
  # Change all pickles to _old
  for files in glob.glob(("*%s.txt" % parameter)):
    newfilename = files+'_old'
    os.rename(files,newfilename)
  os.chdir("..")
  os.chdir("..")

parameter = 'density_of_stars'
par_range = np.array([3,10,20,30,40,50])#np.arange(10,30,)

clean_up_old_files(parameter)

for i in range(len(par_range)):
  command = './master.py %s %s' % (parameter, par_range[i])
  print command
  os.system(command)
  time.sleep(1)
  i+=1
