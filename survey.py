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
  
  # Rearrange observations into an observation_catalog
  if verbose != None: print "Rearranging observation catalog"
  array_size = 0
  for i in range(0,len(data_dic)):
    array_size = array_size + len(data_dic[i][:,0])
  
  if verbose != None: print "Observation catalog requires %d rows)" % array_size
  
  # Copying dictionary into array
  observation_catalog = np.zeros((array_size,6))
  count = 0
  for i in range(0,len(data_dic)):
    single_exposure = data_dic[i]
    for j in range(0,len(single_exposure[:,0])):
      observation_catalog[count,0] = i
      observation_catalog[count,1:6] = single_exposure[j,:]
      count = count +1
      
  return observation_catalog
  # observed_catalog = [pointing ID, star ID, observed_flux, observed_invvar, focal plane x, focal plane y]