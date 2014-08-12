#!/home/dac29/local/bin/python

#Author: Duncan Campbell
#Written: July 29, 2013
#Yale University
#Description: Match the CFHTLS catalogues to the group catalogue

import numpy as np
import h5py
import custom_utilities as cu
from matplotlib import pyplot as plt
import sys
from collections import Counter

def main():
    ###make sure to change these when running in a new enviorment!###
    #location of data directory
    filepath_cat1 = cu.get_output_path() + 'processed_data/yang_groupcat/'
    filepath_cat2 = cu.get_output_path() + 'processed_data/NYU_VAGC/'
    #save data to directory...
    savepath1 = filepath_cat1+'nyu_vagc_match/'
    savepath2 = filepath_cat2+'yang_groupcat_match/'
    #################################################################

    catalogues_1 = ['sample1_L_petro','sample2_L_petro','sample3_L_petro',\
                    'sample1_M_petro','sample2_M_petro','sample3_M_petro',\
                    'sample1_L_model','sample2_L_model','sample3_L_model',\
                    'sample1_M_model','sample2_M_model','sample3_M_model']
    catalogues_2 = ['nyu_vagc_dr7']


    for catalogue in catalogues_1:
        catalogue1 = catalogue
        catalogue2 = catalogues_2[0]
        print catalogue1,'match into', catalogue2
    
        f1 =  h5py.File(filepath_cat1+catalogue1+'.hdf5', 'r')  #open catalogue file
        GC = f1.get(catalogue1)

        f2 =  h5py.File(filepath_cat2+catalogue2+'.hdf5', 'r')  #open catalogue file
        W = f2.get(catalogue2)

        da=2.0*1.0/3600.0 #matching length in degrees
        result_1 = np.array(cu.spherematch(GC['RAgal'], GC['DECgal'], W['RA'], W['DEC'], tol=da, nnearest=1))
        #result_2 = np.array(cu.spherematch(GC['RAgal'], GC['DECgal'], W['RA'], W['DEC'], tol=da, nnearest=2)) #not used
        repeats = [item for item, count in Counter(result_1[1]).iteritems() if count > 1] #double matched objects
        if len(repeats) > 0:
            remove=np.zeros((0,),dtype=np.int) #which entries to remove
            for repeat in repeats:
                result_a = np.where(result_1[1]==repeat)[0] #find indices of the double matched object into catalogue2
                result_b = np.where(result_1[2][result_a]>np.min(result_1[2][result_a]))[0] #find less good matches
                result_c = result_a[result_b] #indices of less good matches into catalogue2
                remove = np.hstack((remove,result_c)) #indices which should be removed
            keep = np.arange(0,len(result_1[0]),1).astype(int)
            keep = np.in1d(keep,remove)==False
            result = result_1[:,keep]
            #unique = np.setdiff1d(result_2[0], result_1[0]) #not used
        else: result = result_1
        
        '''
        x = cu.spheredist(GC['RAgal'],GC['DECgal'],W['RA'][nyu_id],W['DEC'][nyu_id])
        bins = np.arange(0,0.1,da)
        h, bins = np.histogram(x,bins=bins)
        bins = bins*3600.0
        print sum(h[1:])
        print len(GC)-len(result_1[0])
        plt.plot(bins[:-1],h[:])
        plt.yscale('log')
        plt.xlabel('da (deg)')
        plt.ylabel('N')
        plt.title(catalogue1)
        plt.show(block=False)
        '''

        filename1 = catalogue2+'_'+catalogue1+'_match'
        filename2 = catalogue1+'_'+catalogue2+'_match'

        np.save(savepath1+filename1, result[1].astype(int))
        np.save(savepath2+filename2, result[0].astype(int))


if __name__ == '__main__':
  main()
