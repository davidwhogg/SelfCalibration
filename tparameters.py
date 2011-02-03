#!/usr/bin/env python
# -*- coding: utf-8 -*-

dic = {
  'mag_at_ten_sigma'    :     22.5,

  'eta'                 :     0.01,
  
  'm_min'               :     19.5,
  
  'm_max'               :     22.5,
  
  'density_of_stars'    :     10, # all magnitudes
  
  'powerlaw'            :     0.25, # B in log10(dN/dm) = A + B*m
  
  'FoV'                 :     [0.7, 0.7], # [Δα, Δβ]
  
  'ff_samples'          :     [100, 100],
  
  'flat_field_order'    :     2,
  
  'epsilon_max'         :     1,
  
  'sky_limits'          :     [0.0, 10.0, 0.0, 10.0], #[α_min, α_max, β_min, β_max]
  
  }