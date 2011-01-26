# -*- coding: utf-8 -*-

import numpy as np
import camera

def survey(catalog, survey_file, plots, verbose):
  # Routing loads the survey from file and then calls the camera function accordingly.
  
  if verbose != None: print "Loading survey..."
  # pointing = [pointing ID, α, β, orientation]
  pointing = np.loadtxt(survey_file)
  # Calculate number of pointings
  number_pointings = len(pointing[:,0])  
  if verbose != None: print "...done!"
  
  # Declare dictionary
  data_dic = {}
  
  # Perform survey
  if verbose != None: print "Surveying sky..."
  for i in range(0, number_pointings):
    data_dic[i] = camera.camera(catalog, [pointing[i,1],pointing[i,2]], pointing[i,3], plots=plots, verbose = verbose)
    # data_dic[pointing ID] = [star ID, observed_flux, observed_invvar, focal_position]
  if verbose != None: print "...done!"
  
  return data_dic
  # data_dic[pointing ID] = [star ID, observed_flux, observed_invvar, focal_position]