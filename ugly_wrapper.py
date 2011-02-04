#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import time
import string
import datetime
import pickle

from tparameters import dic as simulation_parameters

# Directory to save output data to
directory_path = 'temp'
if directory_path == None:
  directory_path = './' + string.replace(str(datetime.datetime.now()),' ', '.' )
print "Data saved --> ", directory_path
if os.path.exists(directory_path):
  os.system(('rm -r %s' % directory_path))
os.system(('mkdir -p %s' % directory_path))
ff_directory = directory_path+'/Flat_Fields/'
#os.system(('mkdir -p %s' % directory_path))

# Single run, or many runs?
modify_parameter = True
if modify_parameter:
  parameter = 'density_of_stars'
  parameter_values = [5,10,15]
  simulation_parameters[parameter] = parameter_values
else:
  parameter = None
  parameter_values = None

#pickle.dump(simulation_parameters, open(("%s/simulation_parameters.p" % directory_path), "wb" ))

# Call master program
if modify_parameter:
  for i in range(len(parameter_values)):
    new_directory_path = '%s/%s/%s' % (directory_path, parameter, str(parameter_values[i]))
    os.system(('./master.py %s %s %s' % (new_directory_path, parameter, parameter_values[i])))
    time.sleep(1)
    i+=1
  pickle.dump(simulation_parameters, open(("%s/%s/simulation_parameters.p" % (directory_path, str(parameter))), "wb" ))
else:
  new_directory_path = '%s/%s/%s' % (directory_path, 'None', '0')
  os.system(('./master.py %s' % new_directory_path))
  pickle.dump(simulation_parameters, open(("%s/%s/simulation_parameters.p" % (directory_path, 'None')), "wb" ))


    
