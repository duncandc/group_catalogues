#!/home/dac29/local/bin/python

#Author: Duncan Campbell
#Written: July 29, 2013
#Yale University
#Description: Match the CFHTLS catalogues to the group catalogue

import numpy as np
import h5py
import custom_utilities as cu
from collections import Counter

def main():
    ###make sure to change these when running in a new enviorment!###
    #location of data directory
    filepath_cat1 = cu.get_output_path() + 'processed_data/Berland_groupcat/'
    filepath_cat2 = cu.get_output_path() + 'processed_data/mpa_dr7/'
    #save data to directory...
    savepath1 = filepath_cat1+'mpa_dr7_match/'
    savepath2 = filepath_cat2+'berland_groupcat_match/'
    #################################################################

    catalogues_1=['mr19_groups', 'smthresh10.2.groups', 'smthresh9.8.groups']
    catalogues_2=['gal_info_gal_totspecsfr_dr7_v5_2']

    catalogue1=catalogues_1[0]
    catalogue2=catalogues_2[0]
    print catalogue1,'match into', catalogue2
    
    f1 =  h5py.File(filepath_cat1+catalogue1+'.hdf5', 'r')  #open catalogue file
    GC = f1.get(catalogue1)

    f2 =  h5py.File(filepath_cat2+catalogue2+'.hdf5', 'r')  #open catalogue file
    W = f2.get(catalogue2)

    da=2.0*1.0/3600.0 #matching length
    result = cu.spherematch(GC['RA'], GC['DEC'], W['RA'], W['DEC'], tol=da, nnearest=1)
    
    #check to see if anything was matched to the same object
    repeats = [item for item, count in Counter(result[1]).iteritems() if count > 1]
    if len(repeats)>0:
        print 'number of double matched objects:', len(repeats)
    else: print 'all matches are unique.'

    #check to see if every object has a match
    if len(result[0])==len(GC):
        print 'a match was found for every object.'
    else: print 'some objects do not have matches.', len(GC), len(result[0])

    filename1 = catalogue2+'_'+catalogue1+'_match'
    filename2 = catalogue1+'_'+catalogue2+'_match'

    np.save(savepath1+filename1, result[1])
    np.save(savepath2+filename2, result[0])

    catalogue1=catalogues_1[1]
    catalogue2=catalogues_2[0]
    print catalogue1,'match into', catalogue2
    
    f1 =  h5py.File(filepath_cat1+catalogue1+'.hdf5', 'r')  #open catalogue file
    GC = f1.get(catalogue1)

    f2 =  h5py.File(filepath_cat2+catalogue2+'.hdf5', 'r')  #open catalogue file
    W = f2.get(catalogue2)

    da=2.0*1.0/3600.0 #matching length
    result = cu.spherematch(GC['ra'], GC['dec'], W['RA'], W['DEC'], tol=da, nnearest=1)

    #check to see if anything was matched to the same object
    repeats = [item for item, count in Counter(result[1]).iteritems() if count > 1]
    if len(repeats)>0:
        print 'number of double matched objects:', len(repeats)

    #check to see if every object has a match
    if len(result[0])==len(GC):
        print 'a match was found for every object.'
    else: print 'some objects do not have matches.'

    filename1 = catalogue2+'_'+catalogue1+'_match'
    filename2 = catalogue1+'_'+catalogue2+'_match'

    np.save(savepath1+filename1, result[1])
    np.save(savepath2+filename2, result[0])

    catalogue1=catalogues_1[2]
    catalogue2=catalogues_2[0]
    print catalogue1,'match into', catalogue2
    
    f1 =  h5py.File(filepath_cat1+catalogue1+'.hdf5', 'r')  #open catalogue file
    GC = f1.get(catalogue1)

    f2 =  h5py.File(filepath_cat2+catalogue2+'.hdf5', 'r')  #open catalogue file
    W = f2.get(catalogue2)

    da=2.0*1.0/3600.0 #matching length
    result = cu.spherematch(GC['ra'], GC['dec'], W['RA'], W['DEC'], tol=da, nnearest=1)

    #check to see if anything was matched to the same object
    repeats = [item for item, count in Counter(result[1]).iteritems() if count > 1]
    if len(repeats)>0:
        print 'number of double matched objects:', len(repeats)

    #check to see if every object has a match
    if len(result[0])==len(GC):
        print 'a match was found for every object.'
    else: print 'some objects do not have matches.'

    filename1 = catalogue2+'_'+catalogue1+'_match'
    filename2 = catalogue1+'_'+catalogue2+'_match'

    np.save(savepath1+filename1, result[1])
    np.save(savepath2+filename2, result[0])


if __name__ == '__main__':
  main()
