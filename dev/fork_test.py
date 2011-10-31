#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Todo:
#	change density of stars ==> density_of_sources 
#	

import numpy as np
import os
import pickle
import sys
from multiprocessing import Pool
import copy
# Custom Modules

def get_commandline_inputs(): # Read inputs from command line
	if len(sys.argv) == 3:
  		out_dir = sys.argv[1] # directory to save output data to - always needed
  		if out_dir[-1] == '/': out_dir = out_dir[0:-1]
  		param_file = sys.argv[2] # parameter file to use in simulation run
	else:
		print "Error - wrong number of command line arguements!"
		sys.exit() 	
  	return [out_dir, param_file]

if __name__ == '__main__':
	# Make output directory and load parameter dictionary
	out_dir, param_file = get_commandline_inputs()
	os.system(('mkdir -p %s' % out_dir))
	f = open(param_file,'r') 
	params = eval(f.read())

	number_of_runs = 0
	print "Parsing Parameter File: %s" % "test"
	for survey_strategies in params['survey_strategies']:
		for powerlaw_constants in params['powerlaw_constants']:
			for density_of_sources in params['density_of_sources']:
				print powerlaw_constants
				number_of_runs += 1
		
	print params
	print ("%d runs will now be launched..." % number_of_runs)

# Do i need to create new god's? if so do, if not generate one with long filename

# .god
# .results
# .params 
