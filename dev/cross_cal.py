#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Rory Holmes 2011
# Master script to control cross-calibration simulation

# XXXXXX more

import numpy as np
import os
import pickle
import sys

# Custom Modules
import default_parameters
import functions as f
import god

plotdata = True
verbosemode = True

def run_sim(sky, params, strategy, out_dir):
  survey_file = strategy + ".txt"
  observation_catalog = f.survey(params,sky, survey_file, out_dir, plots=plotdata, verbose=verbosemode) 
  # observation_catalog: *.size, *.pointing_ID, *.star_ID, *.flux, *.invvar, *.x, *.y
  
  if plotdata:
    if verbosemode: print "Writing out coverage pickle..."
    f.coverage(params, observation_catalog, strategy, out_dir)
    if verbosemode: print "...done!"  
    if verbosemode: print "Writing out invvar pickle..."
    f.invvar_saveout(observation_catalog, out_dir)
    if verbosemode: print "...done!"  

  # Do cross-calibration
  f.ubercalibration(params, observation_catalog, sky_catalog, strategy, out_dir, plots=plotdata)  
  
# Read command line inputs
if len(sys.argv) > 1:
  out_dir = sys.argv[1] # directory to save output data to - always needed
  os.system(('mkdir -p %s' % out_dir))
  print "Output ==> %s" % out_dir
  if len(sys.argv) > 3: # additional options
    mod_param = sys.argv[2] # parameter to modify
    print ("Modifying Parameter: %s" % mod_param)
    if len(sys.argv) > 4:
      mod_value_low = float(sys.argv[3]) # low value for modified parameter
      mod_value_high = float(sys.argv[4]) # high value for modified parameter
      mod_value_spacing = float(sys.argv[5]) # spacing for modified parameter
      print ("Low Value: %.3f" % mod_value_low)
      print ("High Value: %.3f" % mod_value_high)
      print ("Spacing: %.3f" % mod_value_spacing)
      operating_mode = 3 # modified range of parameters
    else:
      mod_value = sys.argv[3] # OR modified single value
      if mod_value[0] == "[": 
	mod_value = mod_value[1:-1]
	mod_value = np.fromstring(mod_value, sep=',')
      else:
	mod_value = float(mod_value)
      print "New Value: %s" % mod_value 
      operating_mode = 2 # one modified parameter
  else:
    operating_mode = 1 # default parameters
else:
  print "Error - no output directory specified!"
  sys.exit()
  
# Load Default Parameters Dictionary and modify if required
params = default_parameters.dic
if operating_mode == 2:
  params[mod_param] = mod_value
  
# Run Simulations
if ((operating_mode == 1) or (operating_mode == 2)):
  sky_catalog = god.create_catalog(params, out_dir, plots = plotdata, verbose = verbosemode)
  for strategy in params['survey_strategies']:
    os.system('mkdir -p %s/%s' % (out_dir, strategy)) # create directory for survey
    run_sim(sky_catalog, params, strategy, out_dir)
'''
from multiprocessing import Pool

global counter
counter = 0
def cb(r):
    global counter
    print counter, r
    counter +=1
    
def det(M):
    return linalg.det(M)

po = Pool()
for i in xrange(1,300):
    j = random.normal(1,1,(100,100))
    po.apply_async(det,(j,),callback=cb)
po.close()
po.join()
print counter
'''


