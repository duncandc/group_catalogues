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
import custom_utilities as cu


def main():

    filepath = cu.get_data_path() + 'Wetzel_groupcat/'
    savepath = cu.get_output_path() + 'processed_data/Wetzel_groupcat/'
    #################################################################

    names=['ID','ra','dec','z','M_star','Mag_r','Mag_g','SSFR','p_sat','central_ind','M_halo','R_vir','d_proj']

    catalogues=['group_dr7_m.star9.7','group_dr7_m.star10.1','group_dr7_m.star10.6','group_dr7_mag.r19']

    for catalogue in catalogues:
        print catalogue
        filename = catalogue+'.txt'
        print filepath+filename
        data = ascii.read(filepath+filename, delimiter='\s', names=names, \
                      guess=False, Reader=ascii.Basic)
        f = h5py.File(savepath+catalogue+'.hdf5', 'w')
       dset = f.create_dataset(catalogue, data=data)
       f.close()
       gc.collect()


if __name__ == '__main__':
  main()
