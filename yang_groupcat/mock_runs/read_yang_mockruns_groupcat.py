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

def main():

    i = int(sys.argv[1])
    version = sys.argv[2]

    filepath = '/scratch/dac29/data/Yang_groups/mock_runs/4th_run/version_'+version+'/'
    savepath = '/scratch/dac29/output/processed_data/yang_groupcat/mock_runs/4th_run/version_'+version+'/'

    catalogues = ['Mr19_age_distribution_matching_mock_radec_mock',\
                  'Mr19_age_distribution_matching_mock_sys_empty_shuffle_satrel_shuffle_radec_mock']
    filenames_1 = ['imock1_1','imock2_1']
    filenames_2 = ['imock1_2','imock2_2']
    filenames_3 = ['mock1_L_m','mock2_L_m']

    names_1 = ['IDgal','ID','GROUP_ID','cen']
    names_2 = ['GROUP_ID','ID','IDgal']
    names_3 = ['GROUP_ID','ra','dec','z','L19','MGROUP']

    catalogue = catalogues[i]
    filename = filenames_1[i]
    print 'open file 1/3...'
    data_1 = ascii.read(filepath+filename,names=names_1)
    filename = filenames_2[i]
    print 'open file 2/3...'
    data_2 = ascii.read(filepath+filename,names=names_2)
    filename = filenames_3[i]
    print 'open file 3/3...'
    data_3 = ascii.read(filepath+filename,names=names_3)

    print len(data_1), len(data_2), len(data_3)

    #determine index from data_3 into data_1
    index = np.argsort(data_3['GROUP_ID'])
    sorted_3 = data_3['GROUP_ID'][index]
    sorted_index = np.searchsorted(sorted_3,data_1['GROUP_ID'])

    ind = np.take(index, sorted_index, mode="clip")
    mask = data_3['GROUP_ID'][ind] != data_1['GROUP_ID']
    result = np.ma.array(ind, mask=mask)

    #combine different files into one data structure
    dtype=[('gal_ID_1','>i8'),('gal_ID_2','>i8'),('group_ID','>i8'),('brightest','>i8'),\
           ('ra_cen','>f8'),('dec_cen','>f8'),('z_cen','>f8'),('group_L19','>f8'),('halo_mass','>f8')]
    dtype=np.dtype(dtype)
    data = np.recarray((len(data_1),),dtype=dtype)

    #fill in data structure.
    data['gal_ID_1']  = data_1['IDgal']
    data['gal_ID_2']  = data_1['ID']
    data['group_ID']  = data_1['GROUP_ID']
    data['brightest'] = data_1['cen']
    data['ra_cen']    = data_3['ra'][result]
    data['dec_cen']   = data_3['dec'][result]
    data['z_cen']     = data_3['z'][result]
    data['group_L19'] = data_3['L19'][result]
    data['halo_mass'] = data_3['MGROUP'][result]
    
    #save
    print 'saving:', savepath+catalogue+'.hdf5'
    f = h5py.File(savepath+catalogue+'.hdf5', 'w')
    dset = f.create_dataset(catalogue, data=data)
    f.close()
    print 'done.'


if __name__ == '__main__':
  main()
