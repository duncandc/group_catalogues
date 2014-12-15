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

    filepath = cu.get_data_path() + 'Berlind_groupcat/mock_runs/4th_run/'
    savepath = cu.get_output_path() + 'processed_data/berlind_groupcat/mock_runs/4th_run/'
    #################################################################

    catalogues=['Mr19_age_distribution_matching_mock_cen_shuffle_radec_mock.dat',\
                'Mr19_age_distribution_matching_mock_radec_mock.dat',\
                'Mr19_age_distribution_matching_mock_satsys_shuffle_radec_mock.dat',\
                'Mr19_age_distribution_matching_mock_sys_empty_shuffle_cen_shuffle_radec_mock.dat',\
                'Mr19_age_distribution_matching_mock_sys_empty_shuffle_radec_mock.dat',\
                'Mr19_age_distribution_matching_mock_sys_empty_shuffle_satrel_shuffle_radec_mock.dat',\
                'Mr19_age_distribution_matching_mock_sys_empty_shuffle_satsys_shuffle_radec_mock.dat']

    names = ['IDgroup','IDgal','ra','dec','cz']

    for i in range(0,len(catalogues)):
        filename = catalogues[i]
        catalogue = filename[:-4]
        data_1 = ascii.read(filepath+filename,names=names,format='no_header',comment='#')

        #create output file with galaxy ID and group ID
        dtype=[('gal_ID','>i8'),('group_ID','>i8')]
        dtype=np.dtype(dtype)
        data = np.recarray((len(data_1),),dtype=dtype)
        
        #fill in data structure.
        data['gal_ID']  = data_1['IDgal']
        data['group_ID']  = data_1['IDgroup']
    
        #save
        print 'saving:', savepath+catalogue+'.hdf5'
        f = h5py.File(savepath+catalogue+'.hdf5', 'w')
        dset = f.create_dataset(catalogue, data=data)
        f.close()
        print 'done.'
    

    

    
    


if __name__ == '__main__':
  main()
