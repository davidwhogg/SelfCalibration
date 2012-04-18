# Rory Holmes
# Feb 2012

# This file contains the functions to save out data from the
# self-calibration simulations. All data is saved out as pickles.

# Make Python 3 compatible
from __future__ import division, print_function

# Standard Modules
import pickle as pickle

def parameters(data_dir, params, verbose=False):
    ''' Saves out the source catalog object
    
    Input
    -----
    data_dir            :   string
        The output directory path for the simulation run
    params              :   dictionary
        The parameter dictionary
    '''

    filename = '{0}/parameters.p'.format(data_dir)
    if verbose:
        print("Saving out parameters to {0}".format(filename))
    pickle.dump(params, open(filename, "wb"))
    if verbose:
        print("...done!")


def source_catalog(data_dir, source_catalog, verbose=False):
    ''' Saves out the source catalog object
    
    Input
    -----
    data_dir            :   string
        The output directory path for the simulation run
    source_catalog      :   obj
        The catalog object # *.k, *.alpha, *.beta, *.mag, *.size, *.flux,
    '''

    filename = '{0}/source_catalog.p'.format(data_dir)
    if verbose:
        print("Saving out source catalog to {0}".format(filename))
    pickle.dump(source_catalog, open(filename, "wb"))
    if verbose:
        print("...done!")
