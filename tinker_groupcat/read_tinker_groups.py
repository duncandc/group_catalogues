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


def main():
    ###make sure to change these when running in a new enviorment!###############################
    #location of data directory
    filepath = '/scratch/dac29/data/Tinker_groups/mock_runs/4th_run/'
    #save data to directory...
    savepath = '/scratch/dac29/output/processed_data/tinker_groupcat/mock_runs/4th_run/custom_catalogues/'
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
                        data_start=0)
    groups = np.array(groups)
    print 'number of groups:', len(groups)

    #open the satellite probability files
    names = ['foo1','gal_id','group_id','cen_id','M_r','p_sat','lum','foo2','foo3','foo4',\
             'R_proj','d_gal','da_halo']
    filename = catalogue+'.prob'
    prob = ascii.read(filepath+filename, delimiter='\s', names=names, \
                      data_start=0)
    prob = np.array(prob)
    print 'nuber of galaxies:', len(prob)
    print 'number of satellites:', len(np.where(prob['p_sat']>0.5)[0])

    #open the index files
    names = ['foo1','ind','M_r']
    filename = catalogue+'.indx'
    indx = ascii.read(filepath+filename, delimiter='\s', names=names, \
                      data_start=0)
    indx = np.array(indx)
    print 'number of galaxies:', len(indx)
    #############################################################################################
    #open the radec mock
    filename   = mock+'_radec_mock.dat'
    filepath   = '/scratch/dac29/output/processed_data/hearin_mocks/custom_catalogues/'
    radec_mock = ascii.read(filepath+filename, delimiter='\s', Reader=ascii.Basic)
    #open the full mock
    filename  = mock+'.hdf5'
    filepath  = '/scratch/dac29/output/processed_data/hearin_mocks/custom_catalogues/'
    f         =  h5py.File(filepath+filename, 'r')
    full_mock = f.get(mock)
    print full_mock.dtype.names
    #############################################################################################


    #make group catalogue
    dtype=[('ID','>i8'),('RA','>f8'),('DEC','>f8'),\
           ('Z','>f8'),('Z_ERR','>f8'),('Z_TYPE','>i8'),('VELDISP','>f8'),('VELDISP_ERR','>f8'),('FIBERCOL','>i8'),\
           ('M_u,0.1','>f8'),('M_g,0.1','>f8'),('M_r,0.1','>f8'),('M_i,0.1','>f8'),('M_z,0.1','>f8'),\
           ('N_SERSIC','>f8'),\
           ('MSTAR','>f8'),('SSFR','>f8'),\
           ('GROUP_ID','>i8'),('MGROUP','>f8'),('ZGROUP','>f8'),('R200','>f8'),('RPROJ','>f8'),('CEN_IND','>i8'),\
           ('MGROUP_true','<f8'),('cen_true','>i8'),('N_sat','>i8'),('N_sat_red','>i8'),('N_sat_blue','>i8')]
    dtype = np.dtype(dtype)
    data = np.recarray((len(indx),), dtype=dtype)
    data.fill(-99.9) #if no value is available, set = -99.9
    
    ind_full_mock  = radec_mock['k'][indx['ind']-1] #location in full mock
    ind_radec_mock = indx['ind']-1

    data['ID']  = indx['ind']
    data['RA']  = radec_mock['ra'][ind_radec_mock]
    data['DEC'] = radec_mock['dec'][ind_radec_mock]
    data['Z']   = radec_mock['z'][ind_radec_mock]
    data['M_g,0.1'] = full_mock['g-r'][ind_full_mock]+indx['M_r']
    data['M_r,0.1'] = indx['M_r']
    data['MGROUP_true'] = full_mock['M_host'][ind_full_mock]
    ind = np.where(full_mock['ID_host'][ind_full_mock]==-1)[0]
    data['cen_true'][ind] = 1
    ind = np.where(full_mock['ID_host'][ind_full_mock]!=-1)[0]
    data['cen_true'][ind] = 0

    color = data['M_g,0.1']-data['M_r,0.1']
    LHS = 0.7 - 0.032*(data['M_r,0.1']+16.5) #Weinmann 2006
    blue = np.where(color<LHS)[0] #indices of blue galaxies
    red = np.where(color>LHS)[0] #indicies of red galaxies

    #calculate the project seperation in units of kpc/h
    cosmo = cosmology.core.FlatLambdaCDM(100, 0.27) #h=1.0
    data['RPROJ']    = prob['d_gal']*(180.0/math.pi) * 60.0 * (cosmology.funcs.kpc_proper_per_arcmin(data['Z'],cosmo=cosmo))
    data['GROUP_ID'] = prob['group_id']
    #determine group index and central galaxy index in 'groups' and 'prob' arrays
    group_ind = np.empty_like(data['GROUP_ID'])
    cen_ind   = np.empty_like(data['GROUP_ID'])
    Nsat   = np.empty_like(data['GROUP_ID'])
    for i in range(0,len(data)):
        group_ind[i] = np.where(groups['group_id']==data['GROUP_ID'][i])[0]
        cen_ind[i]   = np.where(data['ID']==prob['cen_id'][i])[0]
        sat_ind      = np.where((data['GROUP_ID']==data['GROUP_ID'][i]) & (i==cen_ind[i]))[0]
        Nsat[i]      = len(sat_ind)
        sat_ind_red  = np.where(np.in1d(sat_ind,red)==True)[0]                   
        data['N_sat_red'][i]  = len(sat_ind_red)
        sat_ind_blue = np.where(np.in1d(sat_ind,blue)==True)[0]
        data['N_sat_blue'][i] = len(sat_ind_blue)
        print i, Nsat[i], data['N_sat_red'][i], data['N_sat_blue'][i]
    #calculate group information
    data['MGROUP']   = np.log10(groups['group_mass'][group_ind])
    data['ZGROUP']   = groups['cen_cz'][group_ind]/(299792.458)
    omega_m = 0.27
    data['R200']     = 264.8*((10.0**data['MGROUP'])/(10.0**12.0))**(1.0/3.0)*(1.0+data['ZGROUP'])**(-1) 
    data['CEN_IND']  = cen_ind
    data['N_sat']  = Nsat

    data['ID']  = full_mock['ID_halo'][ind_full_mock]

    print 'saving hdf5 version of the catalogue...'
    filename = mock+'_clf_groups_M19'
    f = h5py.File(savepath+filename+'.hdf5', 'w')
    dset = f.create_dataset(filename, data=data)
    f.close()

    print 'saving ascii version of the catalogue...'
    filename = mock+'_clf_groups_M19'
    data_table = table.table.Table(data=data)
    ascii.write(data_table, savepath+filename+'.dat')
    print data_table




    

if __name__ == '__main__':
    import argparse
    main() 
