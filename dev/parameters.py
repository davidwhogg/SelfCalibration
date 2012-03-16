# This file contains the default parameters for the self-calibration simulations

{
    # Programmatic Parameters

    'data_dir'            :     'default_output' 
    
    'plotdata'            :     True
    
    'verbose'             :     False 

    # Simulation Parameters

    'survey_strategies'   :     ['deep'], #['D', 'C', 'B', 'A']

    'eta'                 :     0.00173214,

    'alpha'               :     0.1584655417,

    'm_min'               :     17,

    'm_max'               :     22,

    'density_of_stars'    :     50, # all magnitudes, number per unit area on sky

    'powerlaw_constants'  :     np.array([-13.34863146, 1.25429311, -0.02122949]), # log10(dN/dm) = A + B*mag + C*mag**2

    'useful_fraction'     :     1.,

    'powerlaw'            :     0.25, # B in log10(dN/dm) = A + B*m

    'FoV'                 :     [0.76, 0.72], # [Δα, Δβ]

    'ff_samples'          :     [100, 100],

    'flat_field_order'    :     8,

    'epsilon_max'         :     1.,

    'sky_limits'          :     [-4.0, 4.0, -4.0, 4.0], #[α_min, α_max, β_min, β_max]

    'seed'                :     1, # Random number seed

    'stop_condition'	    :     1e-5, # Stop cross-cal and best-in-basis fitting when difference in X2 less than this

    'max_iterations'	    :     1049 # the maximum number of cross cal or best-in-basis fitting iterations

    } 
