#!/usr/bin/python

#Author: Duncan Campbell
#Written: August 14, 2013
#Yale University
#Description: Read in hdf5 tinker group catalogues and print out names

###packages###
import numpy as np
from astropy.io import ascii
import h5py
import sys
import glob
import custom_utilities as cu

def main():
  ###make sure to change these when running in a new enviorment!###
  #location of data directory
  filepath = cu.get_output_path() + 'processed_data/tinker_groupcat/'
  #################################################################

  catalogues=['group_dr7_m.star9.7','group_dr7_m.star10.1','group_dr7_m.star10.6','group_dr7_mag.r19']
  filepath = cu.get_output_path() + 'processed_data/tinker_groupcat/'

  for catalogue in catalogues:
      print catalogue
      f =  h5py.File(filepath+catalogue+'.hdf5', 'r')
      dset = f.get(catalogue)
      dset = np.array(dset)
      print 'length:', len(dset)
      for name in dset.dtype.names: print '\t', name

  

if __name__ == '__main__':
  main()
