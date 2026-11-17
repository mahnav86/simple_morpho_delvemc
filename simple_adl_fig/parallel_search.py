#!/usr/bin/env python
"""
Perform simple binning search
Modified for Astrolab by
Mahdieh Navabi
"""
__author__ = "Joanna Sakowska"

import os
import glob
import subprocess
import yaml
import healpy as hp

import simple_adl_fig.search

#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# ~~~~~~~~~~~~~~~~~~~~~ DELVE-MC DWARF SEARCH 2021  ~~~~~~~~~~~~~~~~~~~~~~~~~ #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#

"""
Python parallelisation for general clusters
The csub package is retired in favour of internal python packages
namely multiprocessing

Joanna Sakowska
"""
import pandas as pd
import multiprocessing
from multiprocessing import Pool

# Enter number of threads to use 

n_threads = 3

#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Sub-routines ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#
# JDS: To do: if local cat search, need to add 'cfg' argument I think
#

# Function to tell us job status
# Use if desired

#def info(title):
#    print(title)
#    print('Parent process:', os.getppid())
#    print('Process id:', os.getpid())


def submit_job(healpix_id):
    
     # Output
    outfile = 'results_{}.csv'.format(healpix_id)

    # Commands
    command = 'python {}/search.py --nside 32 --ipix {} --outfile {}'.format(os.path.dirname(simple_adl_fig.search.__file__), healpix_id, outfile)

    # Log
    # Change log directory here
    log_dir = '/user/HS501/mn00962/Documents/simple_adl/simple_adl/log_dir'
    log_name = '{}/results_{}.log'.format(log_dir, healpix_id)

    # Main function

    print(command)

    log = open(log_name, "w")

    search_output = subprocess.run(command.split(' '), shell=False)
 
    print(search_output)
    
    log.write(str(search_output))
    
    log.close()

    return

#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Main search ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--config',type=str,required=False,default='config.yaml',
                        help='config file [.yaml]')
    args = parser.parse_args()

    # Reading config
    with open(args.config, 'r') as ymlfile:
        cfg = yaml.load(ymlfile, Loader=yaml.SafeLoader)

    # # Reading in coordinates
    search_coordinates = pd.read_csv('{}.csv'.format('search_coordinates'))
    print(search_coordinates.columns)
    print('Ready to search the Magellanic Clouds!')
    # MN: To search directly from Astrolab
    for i in search_coordinates['healpix']:
         data_dir = submit_job(i) 

    # # Results directory
    results_dir = '/user/HS501/mn00962/Documents/simple_adl/simple_adl/results_dir_parallel/'

 