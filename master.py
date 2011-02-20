#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import os
import pickle
import sys

# Custom Modules
import functions as f
import god
import tparameters

def init_func():
  if len(sys.argv) < 3:
    pdic = tparameters.dic # parameter database, shared with other modules
    modified_value = None
    modified_parameter = None
  else:
    pdic = tparameters.dic # parameter database, shared with other module
    pdic[sys.argv[2]] = float(sys.argv[3])
    modified_parameter = str(sys.argv[2]) 
    modified_value = float(sys.argv[3])
  directory_path = str(sys.argv[1])

  return (pdic, directory_path, modified_parameter, modified_value)

pdic, directory_path, modified_parameter, modified_value = init_func()

if __name__ == "__main__":
  
  verbose = False
  sky_catalog_plots = True
  plot_invvar = True
  coverage_plots = True
  survey_plots = True
  ff_plots = True


  survey_strategies = pdic['survey_strategies']

  for ii in range(len(survey_strategies)): 
    os.system(('mkdir -p %s/%s' % (directory_path, survey_strategies[ii])))
    
    # Creates catalog of stars used for calibration
    sky_catalog = god.create_catalog(seed = 103, plots = sky_catalog_plots, verbose = verbose)
        # sky_catalog: *.star_ID, *.mag, *.alpha, *.beta, *.size

    # Survey sky
    survey_file = survey_strategies[ii] + ".txt"
    observation_catalog = f.survey(sky_catalog, survey_file, plots=survey_plots, verbose=verbose) 
        # observation_catalog: *.size, *.pointing_ID, *.star_ID, *.flux, *.invvar, *.x, *.y

    if coverage_plots:
      if verbose: print "Writing out coverage pickle..."
      f.coverage(observation_catalog, survey_strategies[ii])
      if verbose: print "...done!"  

    if plot_invvar: 
      if verbose: print "Writing out invvar pickle..."
      f.invvar_saveout(observation_catalog)
      if verbose: print "...done!"  

    # Ubercal
    f.ubercalibration(observation_catalog,sky_catalog,survey_strategies[ii],modified_parameter, modified_value, ff_plots=ff_plots)
    
