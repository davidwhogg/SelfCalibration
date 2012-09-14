# Rory Holmes (MPIA) / David Hogg (NYU)
# 2011 - 2012

# This file contains a sample parameter set for the self-calibration
# simulations. All the parameters are stored in a single dictionary.
# This file is parsed with the eval() function. All these parameters must be
# present for the simulation to run correctly.

{
    # Programmatic Parameters
    # =======================
    # Boolean / String: If set to False, the self-calibration simulations
    # do not save out any data (only returning the results to the calling
    # program). If a string is give, all data is saved out to this directory.
    # NB: old directories with this name are removed!
    'data_dir'            :     True,
    # Boolean: Set to True to run the simulation in verbose mode,
    'verbose'             :     True,

    # Sky Parameters
    # ==============
    # Float: the maximum number of sources (all magnitude) per unit area
    # on the sky. If more than this number are produced by the magnitude
    # distribution (see below), only the brightest are selected.
    'density_of_stars'    :     200.,
    # Float array: The parameters describing the magnitude distribution
    # of the sources in the sky, according to
    # log10(dN/dm) = A + B * mag + C * mag ** 2
    'powerlaw_constants'      :     [-13.34863146, 1.25429311, -0.02122949],
    # Float Array: The area of sky to generate sources in
    # [alpha_min, alpha_max, beta_min, beta_max]
    'sky_limits'          :     [-4.8, 4.8, -4.8, 4.8],
    # Float Array: The area of sky to use for analyses like rms
    # [alpha_min, alpha_max, beta_min, beta_max]
    'analysis_limits'          :     [-3.8, 3.8, -3.8, 3.8],

    # Instrument Parameters
    # =====================
    # Float Array: The simulate imager's field-of-view in degrees [alpha, beta]
    'FoV'                 :     [0.75, 0.75],
    # Float: The saturation limit of the simulated imager
    'm_min'               :     17.,
    # Float: The 10-sigma detection limit of the simulated imager
    'm_max'               :     22.,
    # Floats: The parameters used in the measurement noise model
    # sigma ** 2 = (1 + epsilon) * delta ** 2 + eta ** 2 * count_rate ** 2
    'eta'                 :     0.00173214,
    'delta'               :     0.1584655417,
    # where epsilon is a random number drawn uniformly in the range
    # [0.0, epsilon_max) for each measurement
    'epsilon_max'         :     1.,

    # Simulation Parameters
    # =====================
    # String: The file path to the survey strategy to be performed
    # during the simulation run. The file has the format:
    #       observation_number, RA (deg), Dec (deg), Orientation (deg)
    'survey_file'         :     'quasi_random.srvy',
    # Float array: The number of points to sample the focal plane on for the
    # badness calculations
    'ff_samples'          :     [300, 300],
    # Integer: The order of the flat-field used to fit the instrument response
    # in the self-calibration procedure
    'flat_field_order'    :     8,
    # Float Array: The seed for the random number generator (to ensure repeat
    # runs get the same answer!)
    'seed'                :     [1.],
    # Float: The stop condition for the self-calibration procedure and the
    # best-in-basis fitting (stop when difference is less than 2 times this)
    'stop_condition'	    :   1e-8,
    # Integer: The maximum number of iterations in the self-calibration
    # procedure and the best-in-basis fitting
    'max_iterations'	    :    1001
    }
