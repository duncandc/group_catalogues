#!/usr/bin/python

#Author: Duncan Campbell
#Written: August 14, 2013
#Yale University
#Description: Read in hdf5 berland groupcat catalogue and add/caclulate some quantities from nyu vagc, mpa
#this works on M_r limited samples!

###packages###
import numpy as np
from astropy.io import ascii
import h5py
import sys
import glob
import custom_utilities as cu
from astropy import cosmology
from astropy import table
from astropy.cosmology import FlatLambdaCDM

def main():
    ###make sure to change these when running in a new enviorment!###
    #location of data directory
    filepath1 = cu.get_output_path() + 'processed_data/NYU_VAGC/'
    savepath1 = filepath1 + 'custom_catalogues/'
    filepath2 = cu.get_output_path() + 'processed_data/mpa_dr7/'
    savepath2 = filepath2 + 'custom_catalogues/'
    filepath3 = cu.get_output_path() + 'processed_data/berlind_groupcat/'
    savepath3 = filepath3 + 'custom_catalogues/'
    #################################################################

    cosmo = FlatLambdaCDM(H0=100, Om0=0.3169) #h=1, Omega_Lambda=0.7

    catalogue1 = 'nyu_vagc_dr7'
    catalogue2 = 'gal_info_gal_totspecsfr_dr7_v5_2'
    try: sys.argv[1]
    except: mr = '19'
    else: mr = sys.argv[1]
    catalogue3 = 'mr'+mr+'_groups' 
    print 'reading in', catalogue3, 'catalogue...'
    catalogue3_new = 'sample3_L_model.mr'+mr
    print 'making', catalogue3_new, 'catalogue...'
    
    #open nyu vagc 
    print catalogue1
    f1 =  h5py.File(filepath1+catalogue1+'.hdf5', 'r')
    dset1 = f1.get(catalogue1)
    match13 = np.load(filepath1+'berlind_groupcat_match/'+catalogue3+'_'+catalogue1+'_match.npy')

    #open mpa 
    print catalogue2
    f2 =  h5py.File(filepath2+catalogue2+'.hdf5', 'r')
    dset2 = f2.get(catalogue2)
    match23 = np.load(filepath2+'berlind_groupcat_match/'+catalogue3+'_'+catalogue2+'_match.npy')

    #open groupcat
    print catalogue3
    f3 =  h5py.File(filepath3+catalogue3+'.hdf5', 'r')
    dset3 = f3.get(catalogue3)
    match31 = np.load(filepath3+'nyu_vagc_match/'+catalogue1+'_'+catalogue3+'_match.npy')
    match32 = np.load(filepath3+'mpa_dr7_match/'+catalogue2+'_'+catalogue3+'_match.npy')

    #print dset1.dtype.descr
    #print ' ' 
    #print dset2.dtype.descr
    #print ' ' 
    print dset3.dtype.descr
    print dset3.dtype
    
    #here is the data model for the new group catalogue
    dtype=[('ID','>i8'),('RA','>f8'),('DEC','>f8'),\
           ('Z','>f8'),('Z_ERR','>f8'),('Z_TYPE','>i8'),('VELDISP','>f8'),('VELDISP_ERR','>f8'),('FIBERCOL','>i8'),\
           ('M_u,0.1','>f8'),('M_g,0.1','>f8'),('M_r,0.1','>f8'),('M_i,0.1','>f8'),('M_z,0.1','>f8'),\
           ('N_SERSIC','>f8'),\
           ('MSTAR','>f8'),('SSFR','>f8'),\
           ('GROUP_ID','>i8'),('MGROUP','>f8'),('ZGROUP','>f8'),('R200','>f8'),('RPROJ','>f8'),('CEN_IND','>i8')]
    dtype = np.dtype(dtype)
    
    #define fiber collison galaxies
    result_GC = np.where(dset1['SDSS_SPECTRO_TAG'][match31]==-1)[0] #where in the group catalogue are the collisions
    collision = np.zeros(len(match31), dtype=int)
    collision[result_GC] = 1 #this flag==1 if this is a coollision galaxy

    #create array to store catalogue in
    data = np.recarray((len(dset3),), dtype=dtype)
    data.fill(-99.9) #if no value is available, set = -99.9

    #make galaxy ID's
    ID = np.arange(0,len(data),1).astype(int)

    #input basics
    data['ID']  = dset3['MEMID']
    data['RA']  = dset3['RA']
    data['DEC'] = dset3['DEC']
    data['Z'][match13]      = dset1['Z'][match31]     #redshift if available from anywhere
    data['Z_ERR'][match13]  = dset1['Z_ERR'][match31] #redshift err from sdss
    data['Z_TYPE'][match13] = dset1['ZTYPE'][match31] #source of reshift
    data['VELDISP'][match13]     = dset1['VDISP'][match31]
    data['VELDISP_ERR'][match13] = dset1['VDISP_ERR'][match31]
    data['FIBERCOL'][match13] = collision

    #do K+E corrected ABS magnitude
    AQ=[-4.22,-2.04,-1.62,-1.61,-0.76]
    EQ=AQ[0]*(data['Z']-0.1)
    x = dset1['ABSMAG_u.nearest.model.z0.10'][match31] - EQ
    data['M_u,0.1'] = x
    EQ=AQ[1]*(data['Z']-0.1)
    x = dset1['ABSMAG_g.nearest.model.z0.10'][match31] - EQ
    data['M_g,0.1'] = x
    EQ=AQ[2]*(data['Z']-0.1)
    x = dset1['ABSMAG_r.nearest.model.z0.10'][match31] - EQ
    data['M_r,0.1'] = x
    EQ=AQ[3]*(data['Z']-0.1)
    x = dset1['ABSMAG_i.nearest.model.z0.10'][match31] - EQ
    data['M_i,0.1'] = x
    EQ=AQ[4]*(data['Z']-0.1)
    x = dset1['ABSMAG_z.nearest.model.z0.10'][match31] - EQ
    data['M_z,0.1'] = x

    #apply some spectroscopiclly derived quantities
    data['N_SERSIC'] = dset1['SERSIC_N_r'][match31] 
    data['MSTAR'] = dset3['SM'] #take from Berlind
    data['SSFR'][match23] = dset2['MEDIAN'][match32]

    #add some group properties
    data['GROUP_ID'] = dset3['GROUPID']
    #x = np.column_stack(dset3['GROUPMASS'])
    #data['MGROUP'] = x[2] #0-richness, 1-vel disp, 2-tot r-band lum, 3-cen r-band lum, 4-tot stellar mass 
    data['MGROUP'] = np.log10(dset3['M200B'])

    #identify central galaxy in groups
    group_IDs = np.unique(data['GROUP_ID'])
    for group in group_IDs: #run through each group and identify central galaxy
        #print group
        members = np.where(data['GROUP_ID']==group)[0]
        smallest_mr = np.min(data['M_r,0.1'][members]) #central is the brightest in m_r
        central = np.where((data['GROUP_ID']==group) & (data['M_r,0.1']==smallest_mr))[0][0]
        data['CEN_IND'][members] = central
        data['ZGROUP'][members] = data['Z'][central]
        da = cu.spheredist(data['RA'][central],data['DEC'][central],data['RA'][members],data['DEC'][members])
        da = np.radians(da) #convert to radians
        chi =  cosmology.funcs.comoving_distance(data['ZGROUP'][central], cosmo=cosmo)*1000.0 #in kpc
        dl = cosmology.funcs.luminosity_distance(data['ZGROUP'][central], cosmo=cosmo)*1000.0 #in kpc
        data['RPROJ'][members]=chi/(1.0+data['ZGROUP'][members])*da #caclulate physical seperation
        data['RPROJ'][central]=0.0 #by defintion... this just helps with numerical erros....
        #print group, central, data['RPROJ'][central], len(members) 
    Omega_m = 0.3169   
    x = 258.1 * (10**data['MGROUP']/(10.0**12.0))**(1.0/3.0)*(Omega_m/0.25)**(1.0/3.0)*(1.0+data['ZGROUP'])**(-1.0) 
    data['R200'] = x

    print 'saving hdf5 version of the catalogue...'
    filename = catalogue3_new
    f = h5py.File(savepath3+filename+'.hdf5', 'w')
    dset = f.create_dataset(filename, data=data)
    f.close()

    print 'saving ascii version of the catalogue...'
    filename = catalogue3_new
    data_table = table.table.Table(data=data)
    ascii.write(data_table, savepath3+filename+'.dat')
    print data_table

    


if __name__ == '__main__':
  main()
