#!/usr/bin/python

#Author: Duncan Campbell
#Written: August 14, 2013
#Yale University
#Description: process group finder results on mock runs

###packages###
import numpy as np
from astropy.io import ascii
from astropy import cosmology
from astropy.io import ascii
from astropy import table
import h5py
import sys
import custom_utilities as cu
import math
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy import interpolate

def main():
    savepath = '/scratch/dac29/output/processed_data/berlind_groupcat/mock_runs/4th_run/custom_catalogues/'
    #############################################################################################
    catalogue = sys.argv[1]
    #open the mock group cat
    filepath    = cu.get_output_path()+'processed_data/berlind_groupcat/mock_runs/4th_run/'
    catalogue_1 = catalogue+'_radec_mock'
    f  = h5py.File(filepath+catalogue_1+'.hdf5', 'r')
    GC = f.get(catalogue_1)
    GC = np.array(GC)
    print 'length:', len(GC)
    for name in GC.dtype.names: print '\t', name

    #open the true groups
    filepath    = cu.get_output_path()+'processed_data/hearin_mocks/custom_catalogues/'
    catalogue_2 = catalogue
    f    =  h5py.File(filepath+catalogue_2+'.hdf5', 'r')
    mock = f.get(catalogue_2)
    mock = np.array(mock)
    print 'length:', len(mock)
    for name in mock.dtype.names: print '\t', name

    #open the ra,dec mock fed into the group finder
    filepath   = cu.get_output_path()+'processed_data/hearin_mocks/custom_catalogues/'
    filename   = catalogue+'_radec_mock.dat'
    mock_radec = ascii.read(filepath+filename, delimiter='\s', Reader=ascii.Basic, data_start=1)
    print 'length:', len(mock_radec)
    print mock_radec
    mock_radec = np.array(mock_radec)
    #############################################################################################

    #define some things....
    c = 299792.458 #km/s
    cosmo = cosmology.FlatLambdaCDM(H0=100, Om0=0.27) #h=1
    Omega_m = 0.27
    z_upper_lim = 0.068
    z_lower_lim = 0.020

    #create a new catalogue to store results
    dtype=[('ID','>i8'),('k_1','>i8'),('k_2','>i8'),('RA','>f8'),('DEC','>f8'),('Z','>f8'),('red','>i8'),\
           ('M_u,0.1','>f8'),('M_g,0.1','>f8'),('M_r,0.1','>f8'),('M_i,0.1','>f8'),('M_z,0.1','>f8'),('MSTAR','>f8'),\
           ('GROUP_ID','>i8'),('MGROUP','>f8'),('ZGROUP','>f8'),('R200','>f8'),\
           ('CEN_IND','>i8'),('RANK','>i8'),('RPROJ','>f8'),('N_sat','>i8'),('N_sat_red','>i8'),('N_sat_blue','>i8'),\
           ('HALO_M','>f8'),('HALO_RANK','>i8')]
    dtype = np.dtype(dtype)
    data = np.recarray((len(GC),), dtype=dtype)
    data.fill(-99.9) #empty value indicator

    #catalgues matched by row index
    data['GROUP_ID'] = GC['IDgroup']
    
    #caclulate index into ra-dec mock file
    index = np.argsort(mock_radec['ID'])
    sorted_x = mock_radec['ID'][index]
    ind = np.searchsorted(sorted_x,GC['IDgal'])
    ind = index[ind]

    #grab data from ra-dec mock
    data['k_1']  = ind
    data['ID']  = mock_radec['ID'][ind]
    data['RA']  = mock_radec['ra'][ind]
    data['DEC'] = mock_radec['dec'][ind]
    data['Z']   = mock_radec['z'][ind]

    #calculate index into xyz mock
    ind = mock_radec['k'][ind]
    #grab values from xyz mock
    data['k_2']     = ind
    data['M_g,0.1'] = mock['M_r,0.1'][ind]+mock['g-r'][ind]
    data['M_r,0.1'] = mock['M_r,0.1'][ind]
    data['HALO_M']  = mock['M200b_host'][ind]

    #determine cen/sat designation in xyz mock
    result = np.where(mock['ID_host'][ind]==-1)[0]
    data['HALO_RANK'][result] = 1 #central
    result = np.where(mock['ID_host'][ind]!=-1)[0]
    data['HALO_RANK'][result] = 0 #satellite

    #calculate galaxy colors
    color = data['M_g,0.1']-data['M_r,0.1']
    LHS = 0.7 - 0.032*(data['M_r,0.1']+16.5) #Weinmann 2006
    blue = np.where(color<LHS)[0] #indices of blue galaxies
    red = np.where(color>LHS)[0]  #indicies of red galaxies
    
    #record color designation
    data['red'][red]  = 1
    data['red'][blue] = 0

    for i in range(0,len(data)):
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
        data['ZGROUP'][members]  = data['Z'][central]
        #calculate projected distance from central
        da  = cu.spheredist(data['RA'][central],data['DEC'][central],data['RA'][members],data['DEC'][members])
        da  = np.radians(da) #convert to radians
        chi = cosmology.funcs.comoving_distance(data['ZGROUP'][central], cosmo=cosmo).value*1000.0 #in kpc
        data['RPROJ'][members]=chi/(1.0+data['ZGROUP'][members])*da #caclulate physical seperation
        data['RPROJ'][central]=0.0 #==0 if it is the central

    for i in range(0,len(data)):
        group_id = data['GROUP_ID'][i]
        members = np.where(data['GROUP_ID']==group_id)[0]
        central = np.where(data['M_r,0.1'][members]==min(data['M_r,0.1'][members]))[0][0]
        central = members[central]
        data['CEN_IND'][members] = central

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

    R_200 = 258.1 * (10.0**data['MGROUP']/(10.0**12.0))**(1.0/3.0)*(Omega_m/0.25)**(1.0/3.0)*(1.0+data['ZGROUP'])**(-1.0) 
    data['R200'] = R_200

    print 'saving hdf5 version of the catalogue...'
    filename = catalogue+'_groups'
    f = h5py.File(savepath+filename+'.hdf5', 'w')
    dset = f.create_dataset(filename, data=data)
    f.close()

    print 'saving ascii version of the catalogue...'
    data_table = table.table.Table(data=data)
    ascii.write(data_table, savepath+filename+'.dat')
    print data_table

    
def solar_lum(M,Msol):
    L = ((Msol-M)/2.5)
    return L    

def dist_mod(z, cosmo):
    ld = cosmology.funcs.luminosity_distance(z, cosmo=cosmo).value
    dist_mod = 5.0*(np.log10(ld*1000.0*1000.0)-1.0)
    return dist_mod


if __name__ == '__main__':
  main()
