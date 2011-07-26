#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np

dic = {
  'survey_strategies'     :     ['deep'], #['D', 'C', 'B', 'A']
  
  'mag_at_ten_sigma'    :     22.,

  'eta'                 :     0.001,
  
  'm_min'               :     17,
  
  'm_max'               :     22,
  
  'density_of_stars'    :     25, # all magnitudes, number per unit area on sky
  
  'powerlaw_constants'  :     np.array([-13.34863146, 1.25429311, -0.02122949]), # log10(dN/dm) = A + B*mag + C*mag**2
  
  'useful_fraction'     :     1.,
  
  'powerlaw'            :     0.25, # B in log10(dN/dm) = A + B*m
  
  'FoV'                 :     [0.76, 0.72], # [Δα, Δβ]
  
  'ff_samples'          :     [100, 100],
  
  'flat_field_order'    :     8,
  
  'epsilon_max'         :     1,
  
  'sky_limits'          :     [-4.0, 4.0, -4.0, 4.0], #[α_min, α_max, β_min, β_max]
  
  'seed'                :     1, # Random number seed
  
  'stop_condition'	:     1e-8, # Stop cross-cal when difference in X2 less than this
  
  'max_iterations'	:     1049 # the maximum number of cross cal iterations
  
  } 
