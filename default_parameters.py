#!/usr/bin/env python
# -*- coding: utf-8 -*-

dic = {
  'survey_strategies'     :     ['D'], #['D', 'C', 'B', 'A']
  
  'mag_at_ten_sigma'    :     22.5,

  'eta'                 :     0.01,
  
  'm_min'               :     17, #19.5,
  
  'm_max'               :     19, #22.5,
  
  'density_of_stars'    :     5, # all magnitudes
  
  'powerlaw'            :     0.25, # B in log10(dN/dm) = A + B*m
  
  'FoV'                 :     [0.7, 0.7], # [Δα, Δβ]
  
  'ff_samples'          :     [100, 100],
  
  'flat_field_order'    :     8,
  
  'epsilon_max'         :     1,
  
  'sky_limits'          :     [-4.0, 4.0, -4.0, 4.0], #[α_min, α_max, β_min, β_max]
  
  'seed'                :     1 # Random number seed
  
  } 
