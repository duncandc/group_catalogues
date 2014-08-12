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

    #define some constants
    c = 299792.458 #km/s
    cosmo = cosmology.FlatLambdaCDM(H0=100, Om0=0.27) #h=1
    Omega_m = 0.27
    z_upper_lim = 0.068
    z_lower_lim = 0.020

    #make group catalogue
    dtype=[('ID','>i8'),('k_1','>i8'),('k_2','>i8'),('RA','>f8'),('DEC','>f8'),('Z','>f8'),('red','>i8'),\
           ('M_u,0.1','>f8'),('M_g,0.1','>f8'),('M_r,0.1','>f8'),('M_i,0.1','>f8'),('M_z,0.1','>f8'),('MSTAR','>f8'),\
           ('GROUP_ID','>i8'),('MGROUP','>f8'),('ZGROUP','>f8'),('R200','>f8'),\
           ('CEN_IND','>i8'),('RANK','>i8'),('RPROJ','>f8'),('N_sat','>i8'),('N_sat_red','>i8'),('N_sat_blue','>i8'),\
           ('HALO_M','>f8'),('HALO_RANK','>i8')]
    dtype = np.dtype(dtype)
    data = np.recarray((len(indx),), dtype=dtype)
    data.fill(-99.9) #if no value is available, set = -99.9
    
    #calculate index into xyz mock
    ind_full_mock  = radec_mock['k'][indx['ind']-1] 
    #caclulate index into ra-dec mock file
    ind_radec_mock = indx['ind']-1

    #grab values from ra-dec mock
    data['k_1'] = ind_radec_mock
    data['ID']  = indx['ind']
    data['RA']  = radec_mock['ra'][ind_radec_mock]
    data['DEC'] = radec_mock['dec'][ind_radec_mock]
    data['Z']   = radec_mock['z'][ind_radec_mock]

    #grab values from xyz mock
    data['ID']  = full_mock['ID_halo'][ind_full_mock]
    data['M_g,0.1'] = full_mock['g-r'][ind_full_mock]+indx['M_r']
    data['M_r,0.1'] = indx['M_r']
    data['HALO_M'] = full_mock['M_host'][ind_full_mock]
    ind = np.where(full_mock['ID_host'][ind_full_mock]==-1)[0]
    data['HALO_RANK'][ind] = 1
    ind = np.where(full_mock['ID_host'][ind_full_mock]!=-1)[0]
    data['HALO_RANK'][ind] = 0

    color = data['M_g,0.1']-data['M_r,0.1']
    LHS = 0.7 - 0.032*(data['M_r,0.1']+16.5) #Weinmann 2006
    blue = np.where(color<LHS)[0] #indices of blue galaxies
    red = np.where(color>LHS)[0]  #indicies of red galaxies

    #record color designation
    data['red'][red]  = 1
    data['red'][blue] = 0

    #calculate the project seperation in units of kpc/h
    data['RPROJ']    = prob['d_gal']*(180.0/math.pi) * 60.0 * (cosmology.funcs.kpc_proper_per_arcmin(data['Z'],cosmo=cosmo))
    data['GROUP_ID'] = prob['group_id']
    
    group_ind = np.empty_like(data['GROUP_ID'])
    for i in range(0,len(data)):
        group_ind[i] = np.where(groups['group_id']==data['GROUP_ID'][i])[0]
        group_id = data['GROUP_ID'][i]
        members  = np.where(data['GROUP_ID']==group_id)[0]
        central  = np.where(data['M_r,0.1'][members]==min(data['M_r,0.1'][members]))[0][0]
        central  = members[central]
        satellites = np.where(members!=central)[0]
        satellites = members[satellites]
        #record rank
        data['RANK'][central]    = 1
        data['RANK'][satellites] = 0
        #record number of satellites in the group
        data['N_sat'][members]   = len(satellites)
        sat_red  = np.where(np.in1d(satellites,red)==True)[0]                   
        data['N_sat_red'][members]  = len(sat_red)
        sat_blue = np.where(np.in1d(satellites,blue)==True)[0]
        data['N_sat_blue'][members] = len(sat_blue)
        #record other group information
        data['CEN_IND'][members] = central
    #calculate group information
    data['MGROUP']   = np.log10(groups['group_mass'][group_ind])
    data['ZGROUP']   = groups['cen_cz'][group_ind]/(299792.458)
    data['R200']     = 258.1 * (10.0**data['MGROUP']/(10.0**12.0))**(1.0/3.0)*(Omega_m/0.25)**(1.0/3.0)*(1.0+data['ZGROUP'])**(-1.0) 

    print 'check:', min(data['Z'])>=z_lower_lim
    print 'check:', max(data['Z'])<=z_upper_lim

    #read in mass halo function
    filepath = '/scratch/dac29/fortran_code/mass_functions/'
    filename = 'Bolshoi_Massfunc.dat'
    names = ['dM','dn','nsum']
    dndM = ascii.read(filepath+filename, delimiter='\s', names=names, data_start=0)
    dndM = np.array(dndM)

    #idenditify centrals and satellites
    centrals   = np.where(data['RPROJ']==0)[0]
    satellites = np.where(data['RPROJ']>0)[0]

    #calculate group total r-band luminosities
    S_r = 4.64
    group_L = np.zeros((len(data),),dtype=np.float)

    for i  in range(0,len(centrals)):
        gal = np.where(data['GROUP_ID']==data['GROUP_ID'][centrals[i]])[0]
        group_L[gal] = np.log10(np.sum(10.0**(solar_lum(data['M_r,0.1'][gal], S_r))))
    tot_lum = group_L[centrals]

    #calculate abundance matched masses for groups
    geo_f = 1.0/8.0 #gemoetric factor: spherical octant
    r_max = cosmology.funcs.comoving_distance(z_upper_lim, cosmo=cosmo).value #in Mpc
    r_min = cosmology.funcs.comoving_distance(z_lower_lim, cosmo=cosmo).value #in Mpc
    mock_volume = (4.0/3.0)*math.pi*(r_max**3.0-r_min**3)*geo_f

    #caclulate the group luminosity function
    N_gal = np.cumsum(np.zeros(len(centrals))+1) #cumulative number of groups
    n_gal = N_gal/mock_volume #number density
    L_gal = np.sort(tot_lum)[::-1] #group luminosity
    ind   = np.argsort(tot_lum)[::-1]

    #integrate halo mass function
    n_halo  = dndM['nsum'][::-1] #cumulative number desnity
    M_halo  = dndM['dM'][::-1] #halo mass

    #interpolate the halo mass function
    x = np.log10(n_halo)
    y = M_halo
    f = interpolate.interp1d(x, y, kind='cubic', bounds_error='False', fill_value=0.0)

    data['MGROUP'][centrals[ind]] = f(np.log10(n_gal))

    for i  in range(0,len(centrals)):
        gal = np.where(data['GROUP_ID']==data['GROUP_ID'][centrals[i]])[0]
        data['MGROUP'][gal] = data['MGROUP'][centrals[i]]

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


def solar_lum(M,Msol):
    L = ((Msol-M)/2.5)
    return L 

    

if __name__ == '__main__':
    import argparse
    main() 
