#!/usr/bin/python

#Author: Duncan Campbell
#Written: July 9, 2013
#Yale University
#Description: Read in fits group catalogues and save as HDF5 files.

###packages###
import numpy as np
from astropy.io import ascii
import h5py
import sys
import custom_utilities as cu

def main():

    catalogue = 'Mr19_age_distribution_matching_mock_radec_mock'
    
    filepath = cu.get_output_path()+'processed_data/tinker_groupcat/mock_runs/4th_run/'
    
    f = h5py.File(filepath+catalogue+'.hdf5', 'r')
    GC = f.get(catalogue)
    
    print GC.dtype.names
    print len(GC)
    print GC['gal_ID']
    print GC['group_ID']

if __name__ == '__main__':
  main()