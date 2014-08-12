#!/usr/bin/python

#Author: Duncan Campbell
#Written: July 9, 2013
#Yale University
#Description: Read in fits group catalogues and save as HDF5 files.

###packages###
import numpy as np
from astropy.io import fits
from astropy.io import ascii
import h5py
import gc
import custom_utilities as cu
import sys


def main():

    filepath = cu.get_data_path() + 'Berlind_groupcat/mock_runs/3rd_run/'
    savepath = cu.get_output_path() + 'processed_data/berlind_groupcat/mock_runs/3rd_run/'
    #################################################################

    catalogues=['Mr19_age_distribution_matching_mock_cen_shuffle_radec_mock.dat',\
                'Mr19_age_distribution_matching_mock_radec_mock.dat',\
                'Mr19_age_distribution_matching_mock_satsys_shuffle_radec_mock.dat',\
                'Mr19_age_distribution_matching_mock_sys_empty_shuffle_cen_shuffle_radec_mock.dat',\
                'Mr19_age_distribution_matching_mock_sys_empty_shuffle_radec_mock.dat',\
                'Mr19_age_distribution_matching_mock_sys_empty_shuffle_satrel_shuffle_radec_mock.dat',\
                'Mr19_age_distribution_matching_mock_sys_empty_shuffle_satsys_shuffle_radec_mock.dat']

    #names = ['IDgroup','richness','ra_center','dec_center','z_center','sigma_v','rms_size','foo1','foo1b','IDgal',\
    #         'ra','dec','z','foo2','foo3','foo4','foo5','R_proj']
    names = ['IDgroup','IDgal','ra','dec','z']

    N = len(catalogues)
    for i in range(0,N):
        filename = catalogues[i]
        data = ascii.read(filepath+filename,names=names)
        print data
        print savepath+filename[:-4]+'.hdf5'
        f = h5py.File(savepath+filename[:-4]+'.hdf5', 'w')
        dset = f.create_dataset(filename[:-4], data=data)
        f.close()


    

    
    


if __name__ == '__main__':
  main()
