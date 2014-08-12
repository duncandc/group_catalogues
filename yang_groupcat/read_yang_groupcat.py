#!/usr/bin/python

#Author: Duncan Campbell
#Written: July 9, 2013
#Yale University
#Description: Read in ascii group catalogues and save as HDF5 files.

###packages###
import numpy as np
from astropy.io import ascii
import h5py
import gc


def main():
  ###make sure to change these when running in a new enviorment!###
  #location of data directory
  filepath = '/scratch/dac29/data/groupcat_DR7_v4.0/'
  #save data to directory...
  savepath = '/scratch/dac29/output/processed_data/yang_groupcat/'
  #################################################################
  
  names=['IGAL','k','ID','RAgal','DECgal','ZGAL','ZTYPE','LUMDIST','COMPL',\
         'VELDISP','VELDISP_ERR','R50_r','CONC_r','appmag_u','appmag_g',\
         'appmag_r','appmag_i','appmag_z','M_u,0.1','M_g,0.1','M_r,0.1',\
         'M_i,0.1','M_z,0.1','Mstar','CTYPE','GROUP_ID','GROUP_RA','GROUP_DEC',\
         'GROUP_Z','L19.5','STELMAS','Mgroup','Ngroup','Cgroup','RANK',\
         'Rproj_LW','Rproj_L','weight_L','zmax_L','weight_M','zmax_M',\
         'N_sersic']

  catalogues=['sample1_L_petro','sample2_L_petro','sample3_L_petro',\
              'sample1_M_petro','sample2_M_petro','sample3_M_petro',\
              'sample1_L_model','sample2_L_model','sample3_L_model',\
              'sample1_M_model','sample2_M_model','sample3_M_model']

  for catalogue in catalogues:
    print catalogue
    filename = catalogue+'.dat'
    data = ascii.read(filepath+filename, delimiter='\s', names=names, \
                      guess=False, Reader=ascii.Basic, data_start=0)
    f = h5py.File(savepath+catalogue+'.hdf5', 'w')
    dset = f.create_dataset(catalogue, data=data)
    f.close()
    gc.collect()


if __name__ == '__main__':
  main()
