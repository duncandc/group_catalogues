#!/usr/bin/env python

#Author: Duncan Campbell
#Written: July 9, 2013
#Yale University
#Description: Read in ascii group catalogues and save as HDF5 files.

###packages###
import numpy as np
from astropy.io import ascii
from astropy import table
from astropy import cosmology
import h5py
import sys
import math
from scipy.interpolate import interp1d
from scipy import interpolate
import custom_utilities as cu


def main():

    
    filepath = cu.get_data_path()+'Tinker_groupcat/mock_runs/4th_run/'
    savepath = cu.get_output_path()+'processed_data/tinker_groupcat/mock_runs/4th_run/'
    
    #############################################################################################
    i = int(sys.argv[1])
    catalogues = ['clf_groups_M19_1','clf_groups_M19_2','clf_groups_M19_3','clf_groups_M19_4',\
                  'clf_groups_M19_5','clf_groups_M19_6','clf_groups_M19_7']
    mocks      = ['Mr19_age_distribution_matching_mock_sys_empty_shuffle_satsys_shuffle',\
                  'Mr19_age_distribution_matching_mock_sys_empty_shuffle_satrel_shuffle',\
                  'Mr19_age_distribution_matching_mock_sys_empty_shuffle_cen_shuffle',\
                  'Mr19_age_distribution_matching_mock_sys_empty_shuffle',\
                  'Mr19_age_distribution_matching_mock',\
                  'Mr19_age_distribution_matching_mock_satsys_shuffle',\
                  'Mr19_age_distribution_matching_mock_cen_shuffle']
    catalogue = catalogues[i]
    mock = mocks[i]
    #############################################################################################

    #open group files
    names=['foo1','group_id','cen_id','group_mass','group_mass_prev','n_sat','l_tot',\
          'l_central','foo2','cen_ra','cen_dec','cen_cz','foo3','foo4']
    filename = catalogue+'.groups'
    groups = ascii.read(filepath+filename, delimiter='\s', names=names, \
                        data_start=0,format='no_header')
    groups = np.array(groups)
    print 'number of groups:', len(groups)

    #open the satellite probability files
    names = ['foo1','gal_id','group_id','cen_id','M_r','p_sat','lum','foo2','foo3','foo4',\
             'R_proj','d_gal','da_halo']
    filename = catalogue+'.prob'
    prob = ascii.read(filepath+filename, delimiter='\s', names=names, \
                      data_start=0,format='no_header')
    prob = np.array(prob)
    print 'nuber of galaxies:', len(prob)
    print 'number of satellites:', len(np.where(prob['p_sat']>0.5)[0])

    #open the index files
    names = ['foo1','ind','M_r']
    filename = catalogue+'.indx'
    indx = ascii.read(filepath+filename, delimiter='\s', names=names, \
                      data_start=0,format='no_header')
    indx = np.array(indx)
    print 'number of galaxies:', len(indx)
    #############################################################################################
    #open the radec mock
    filename   = mock+'_radec_mock.dat'
    filepath   = cu.get_output_path()+'processed_data/hearin_mocks/custom_catalogues/'
    radec_mock = ascii.read(filepath+filename, delimiter='\s', Reader=ascii.Basic)
    #open the full mock
    filename  = mock+'.hdf5'
    filepath  = cu.get_output_path()+'processed_data/hearin_mocks/custom_catalogues/'
    f         =  h5py.File(filepath+filename, 'r')
    full_mock = f.get(mock)
    print full_mock.dtype.names
    #############################################################################################

    ind_full_mock  = radec_mock['k'][indx['ind']-1] 

    #create output file with galaxy ID and group ID
    dtype=[('gal_ID','>i8'),('group_ID','>i8')]
    dtype=np.dtype(dtype)
    data = np.recarray((len(indx),), dtype=dtype)

    #fill in data structure.
    data['gal_ID']  = full_mock['ID_halo'][ind_full_mock]
    data['group_ID']  = prob['group_id']
    
    #save
    catalogue = mock
    print 'saving:', savepath+catalogue+'_radec_mock.hdf5'
    f = h5py.File(savepath+catalogue+'_radec_mock.hdf5', 'w')
    dset = f.create_dataset(catalogue+'_radec_mock', data=data)
    f.close()
    print 'done.'

if __name__ == '__main__':
    import argparse
    main() 
