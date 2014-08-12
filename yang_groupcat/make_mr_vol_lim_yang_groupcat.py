#!/usr/bin/python

#Author: Duncan Campbell
#Written: July 30, 2013
#Yale University
#Description: make a volume lmited M_r groupcatalogue sample.

import numpy as np
import h5py
import matplotlib.pyplot as plt
import custom_utilities as cu
from astropy import cosmology
from scipy import interpolate
import sys


def main():
    filepath_cat = cu.get_output_path() + 'processed_data/yang_groupcat/custom_catalogues/'
    savepath = filepath_cat

    abs_mag_lim = float(sys.argv[1])
    print abs_mag_lim
    
    catalogues=['sample1_L_model', 'sample2_L_model', 'sample3_L_model'] #groupcat catalogues

    z = np.arange(0.01,0.21,0.01)
    lim=np.zeros(len(z))
    for i in range(0, len(z)):
        print z[i]
        lim[i] = abs_mag_lim_func(z[i])
    print z
    print lim
    func = interpolate.interp1d(lim[::-1], z[::-1], kind='cubic')

    for i in range(0,len(catalogues)):
        catalogue = catalogues[i]
        #open catalogue
        f =  h5py.File(filepath_cat+catalogue+'.hdf5', 'r')  #open catalogue file
        GC = f.get(catalogue)
        print 'read in',catalogue+'.'

        zmax = func(abs_mag_lim)

        selection_1 = GC['ZGROUP']<zmax
        selection_2 = GC['M_r,0.1']<abs_mag_lim
        selection = selection_1*selection_2
        sample = GC[selection]

        filename = catalogue+'.mr'+str(int(abs(abs_mag_lim)))
        f = h5py.File(savepath+filename+'.hdf5', 'w')
        dset = f.create_dataset(filename, data=sample)
        f.close()

        plt.plot(GC['ZGROUP'],GC['M_r,0.1'],'.', color='black', alpha=0.1)
        plt.plot(sample['ZGROUP'],sample['M_r,0.1'],'.')
        plt.plot(z,lim, color='red')
        plt.ylim([-15,-25])
        plt.xlabel('group z')
        plt.ylabel('M_r,0.1')
        plt.show(block=True)

def abs_mag_lim_func(z):
    #redshift dep absolute magnitude lim corrected to z=0.1
    dl = lum_dist(z)
    DM = dist_mod(dl)
    k = k_corr(z)
    lim = 17.77-DM-k+1.62*(z-0.1)-0.1
    return lim

def dist_mod(dl):
    DM = 5.0*np.log10(dl)+25.0
    return DM

def lum_dist(z):
    cosmo = cosmology.FlatLambdaCDM(H0=100, Om0=0.3) #h=1, Omega_m=0.3, Omega_Lambda=0.7
    dl = cosmology.funcs.luminosity_distance(z, cosmo=cosmo) #in Mpc
    return dl

def k_corr(z):
    k = 2.5*np.log10((z+0.9)/1.1)
    return k

def abs_mag(app_mag, z):    
    cosmo = cosmology.FlatLambdaCDM(H0=100, Om0=0.3) #h=1, Omega_m=0.3, Omega_Lambda=0.7
    ld = cosmology.funcs.luminosity_distance(z, cosmo=cosmo)
    Mag=app_mag-5.0*(np.log10(ld*1000.0*1000.0)-1.0)
    
    return Mag   
    

if __name__ == '__main__':
    main()  
