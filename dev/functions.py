#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Rory Holmes 2011
# Contains all the functions used in the cross-calibration simulation [./cross_cal.py]

import numpy as np
import os
import pickle
import string
import scipy.optimize as sci
# Custom Modules
import god
import save_out
import transformations as tran 
import single_image



#*****************************************************
#*************** Best Fit of God's *******************
#*****************************************************



def bestfit_ff(params, out_dir):
  order = params['flat_field_order']
  ff_samples = params['ff_samples']
  FoV = params['FoV']
  temp_x = np.linspace(-0.5*FoV[0], 0.5*FoV[0], ff_samples[0])
  temp_y = np.linspace(-0.5*FoV[1], 0.5*FoV[1], ff_samples[1])
  X, Y = np.meshgrid(temp_x,temp_y)
  x = np.reshape(X,-1)
  y = np.reshape(Y,-1)
  g = evaluate_flat_field_functions(x,y, order)
  god_ff = god.flat_field(params,x,y)
  a = np.zeros((order + 1)*(order + 2)/2)
  a[0] = 1
  a[3] = -0.2
  a[5] = 0.5
  print "Fitting god's flat-field with basis..."
  fitted_parameters = sci.fmin_bfgs(compare_flats, a, args = (god_ff, g, x, y), gtol = 1e-8, maxiter = 1e16)
  print "Fitted Parameters: ", fitted_parameters
  "...done!"
  fitted_parameters = normalize_flat_field(params, fitted_parameters)
  bestfit_ff = np.reshape(evaluate_flat_field(params, x, y, fitted_parameters),
(len(X[:,0]),len(X[0,:])))
  dic = {}
  dic['x'] = X
  dic['y'] = Y
  dic['bestfit_ff'] = bestfit_ff
  dic['fit_parameters'] = fitted_parameters
  filename = '%s/bestfit_ff.p' % (out_dir)
  pickle.dump(dic, open(filename, "wb"))
