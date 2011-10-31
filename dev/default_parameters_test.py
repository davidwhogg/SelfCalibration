# TODO
# 	Change alpha and beta
# 	Add mode = general, density compare, universal plot

# Parameters for self-calibration simulations
{
########################################################
############## Simulation Parameters ###################
########################################################

# Simulation mode: 
#		- 'general' = performs self-calibration for 
#		- 'density_compare' = creates density plots comparing TWO parameters -- XX which?, uses 1st value for each parameter 
#		- 'universal' = creates a universal plot?!
	'simulation_mode' 		:			'general'
# Survey Strategies to investigate  
  'survey_strategies'   :     ['D', 'C', 'A'], #['D', 'C', 'B', 'A']

########################################################
################## God Parameters ######################
########################################################

# Random number generation seed for sky simulation
  'seed'                :     1, # Random number seed

# Density of Sources to include in the simulations (only brightest sources are used)
  'density_of_sources'    :     np.linspace(1, 100, 10).astype('int'), # can run out of memory quite quickly! 100 is ok, 1000 can be run on broiler but is quite slow

# Distribution of source magnitudes: log10(dN/dm) = A + B*mag + C*mag**2
	'powerlaw_constants'  :     [np.array([-13.34863146, 1.25429311, -0.02122949]), np.array([1., 2., 3.])], # fit to Windhorst et al 2011 

# Region of sky to survey
	'sky_limits'          :     [np.array(-4.0, 4.0, -4.0, 4.0)], 
  
########################################################
############### Instrument Parameters ##################
########################################################

# Saturation Limit
	'm_min'               :     np.array([17]),

# 10 sigma detection limit
	'm_max'               :     np.array([22]),  

# Instrument Field-of-View
	'FoV'                 :     [np.array([0.76, 0.72])], 

# Noise Model
	'noise'								:			np.array([0.00173214, 0.1584655417]) # [eta, alpha]

# Additional noise term not taken into account in the analysis
   'epsilon_max'         :     np.array([1.]),
  
########################################################
################ Analysis Parameters ###################
########################################################

# Number of samples of instrument response used during the badness calculations  
  'ff_samples'          :     [300, 300],

# What order polynomial to fit instrument response?
  'flat_field_order'    :     8,
  
# Fraction of sources that can be used in the self-calibration procedure
	'useful_fraction'			:			1.,
  
# Stop self-cal when X^2 changes by less than:   
  'stop_condition'	:     1e-8, 

# Maximum number of self-calibration iterations before giving up 
  'max_iterations'	:     1049 
} 
