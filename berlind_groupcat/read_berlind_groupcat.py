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


def main():

    filepath = cu.get_data_path() + 'Berlind_groupcat/'
    savepath = cu.get_output_path() + 'processed_data/berlind_groupcat/'
    #################################################################

    catalogues=['mr19_groups.fits', 'smthresh10.2.groups.dat', 'smthresh9.8.groups.dat']

    filename = catalogues[0]
    hdulist = fits.open(filepath+filename, memmap=True)
    data = hdulist[1].data
    print 'saving as:', savepath+filename[:-5]+'.hdf5'
    f = h5py.File(savepath+filename[:-5]+'.hdf5', 'w')
    dset = f.create_dataset(filename[:-5], data=data)
    f.close()
    gc.collect()

    names = ['ra','dec','z','groupID','rank','Mstar','Mr','SSFR','Mgroup']

    filename = catalogues[1]
    data = ascii.read(filepath+filename, guess=True, Reader=ascii.Basic, names=names,data_start=0)
    print 'saving as:', savepath+filename[:-4]+'.hdf5'
    f = h5py.File(savepath+filename[:-4]+'.hdf5', 'w')
    dset = f.create_dataset(filename[:-4], data=data)
    f.close()
    gc.collect()

    filename = catalogues[2]
    data = ascii.read(filepath+filename, guess=True, Reader=ascii.Basic, names=names,data_start=0)
    print 'saving as:', savepath+filename[:-4]+'.hdf5'
    f = h5py.File(savepath+filename[:-4]+'.hdf5', 'w')
    dset = f.create_dataset(filename[:-4], data=data)
    f.close()
    gc.collect()

    


if __name__ == '__main__':
  main()
