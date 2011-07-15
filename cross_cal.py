#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Rory Holmes 2011
# Master script to control cross-calibration simulation
# Requires default_parameters.py, functions.py and plot.py

# Can call from the command line in 3 ways:
#	1. ./cross_cal [output_dir]
#	2. ./cross_cal [output_dir] [default_parameter_to_modify] [new_value]
#	3. ./cross_cal [output_dir] [default_parameter_to_modify] [low_range_value] [high_range_value] [number_of_samples_within_range]

import numpy as np
import os
import pickle
import sys
# Custom Modules
import default_parameters
import functions as f
import god

plotdata = True
verbosemode = False

def run_sim(sky, params, strategy, out_dir, data_dir):
  os.system('mkdir -p %s' % (data_dir))
  os.system('mkdir -p %s' % (data_dir + '/FF'))
  pickle.dump(params, open((data_dir+'/parameters.p'), "wb"))
  # Check to see if God and bestfit flat-fields have been saved already
  if os.path.isfile(out_dir+ '/bestfit_ff.py') == False:
    god.saveout_bestfit_ff(params, out_dir)
  if os.path.isfile(out_dir+ '/god_ff.py') == False:
    god.saveout_god_ff(params, out_dir)

  survey_file = strategy + ".txt"
  observation_catalog = f.survey(params,sky, survey_file, data_dir, plots=plotdata, verbose=verbosemode) 
  # observation_catalog: *.size, *.pointing_ID, *.star_ID, *.flux, *.invvar, *.x, *.y
  
  if plotdata:
    if verbosemode: print "Writing out coverage pickle..."
    f.coverage(params, observation_catalog, strategy, data_dir)
    if verbosemode: print "...done!"  
    if verbosemode: print "Writing out invvar pickle..."
    f.invvar_saveout(observation_catalog, data_dir)
    if verbosemode: print "...done!"  

  # Do cross-calibration
  sln = f.ubercalibration(params, observation_catalog, sky_catalog, strategy, data_dir, plots=plotdata)
  return sln # [iteration number, rms, badness, chi2]
  
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
      no_mod_value = int(sys.argv[5]) # spacing for modified parameter
      print ("Low Value: %.3f" % mod_value_low)
      print ("High Value: %.3f" % mod_value_high)
      print ("Number of Values in Range: %.3f" % no_mod_value)
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

if __name__ == '__main__':
  # Make output directory - remove old results
  os.system('mkdir -p %s' % (out_dir))
  os.system('rm -r %s/*' % (out_dir))

  # Load Default Parameters Dictionary and modify if required
  params = default_parameters.dic
  # Run Simulations
  if ((operating_mode == 1) or (operating_mode == 2)):
    for strategy in params['survey_strategies']:
      os.system('mkdir -p %s/%s' % (out_dir, strategy)) # create directory for survey
      if (operating_mode == 1):
	data_dir = ('%s/%s/default' % (out_dir,strategy))
      elif (operating_mode == 2):
	if mod_param in params:
	  params[mod_param] = mod_value
	  data_dir = ('%s/%s/%s=%.3f' % (out_dir,strategy, mod_param, mod_value))
	else: 
	  print "Error! Unknown parameter specified..." 
	  sys.exit()
      
      else: print "Error - operating mode not defined!"
      sky_catalog = god.create_catalog(params, out_dir, plots = plotdata, verbose = verbosemode)
      run_sim(sky_catalog, params, strategy, out_dir, data_dir)
  elif (operating_mode == 3):
    if (mod_param == 'FoV') or (mod_param == 'ff_samples') or (mod_param == 'flat_field_order'):
      print "Error! Unsuitable parameter selected..."
      sys.exit()
    if mod_param in params:
      param_range = np.linspace(mod_value_low, mod_value_high, num=no_mod_value, endpoint=True, retstep=False)
      if (mod_param == 'flat_field_order'): param_range = param_range.astype('int') 
      for strategy in params['survey_strategies']:
	os.system('mkdir -p %s/%s' % (out_dir, strategy)) # create directory for survey
	sln_dic = {}
	sln_dic['modified_parameter'] = mod_param
	sln_dic['parameter_range'] = param_range
	sln_it = 0*param_range
	sln_rms = 0*param_range
	sln_bdnss = 0*param_range
	sln_chi2 = 0*param_range
	for indx in range(len(param_range)):
	  params[mod_param] = param_range[indx]
	  sky_catalog = god.create_catalog(params, out_dir, plots = plotdata, verbose = verbosemode)
	  data_dir = ('%s/%s/%s=%.3f' % (out_dir, strategy, mod_param, param_range[indx]))
	  sln_it[indx], sln_bdnss[indx], sln_rms[indx], sln_chi2[indx] = run_sim(sky_catalog, params, strategy, out_dir, data_dir)      

	sln_dic['it_num'] = sln_it
	sln_dic['rms'] = sln_rms
	sln_dic['bdnss'] = sln_bdnss
	sln_dic['chi2'] = sln_chi2
	pickle.dump(sln_dic, open((out_dir + '/' + strategy +'/solution.p'), "wb"))
    else: 
      print "Error! Unknown parameter specified..." 
      sys.exit()
  else:
    print "Error - no operating mode defined!"

os.system('./plot.py %s' % out_dir)


