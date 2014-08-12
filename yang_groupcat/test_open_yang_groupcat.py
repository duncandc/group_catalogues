#!/usr/bin/python

#Author: Duncan Campbell
#Written: August 14, 2013
#Yale University
#Description: Read in hdf5 wetzel groupcat catalogues and print out names

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
  filepath = cu.get_output_path() + 'processed_data/yang_groupcat/catalogues/'
  #################################################################

  catalogues=['sample1_L_petro','sample2_L_petro','sample3_L_petro',\
              'sample1_M_petro','sample2_M_petro','sample3_M_petro',\
              'sample1_L_model','sample2_L_model','sample3_L_model',\
              'sample1_M_model','sample2_M_model','sample3_M_model']

  for catalogue in catalogues:
      print 'reading in:', catalogue
      f =  h5py.File(filepath+catalogue+'.hdf5', 'r')
      dset = f.get(catalogue)
      dset = np.array(dset)
      print 'length:', len(dset)
      for name in dset.dtype.names: print '\t', name

  

if __name__ == '__main__':
  main()
