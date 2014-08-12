#!/usr/bin/python

#Author: Duncan Campbell
#Written: August 14, 2013
#Yale University
#Description: Read in hdf5 yang groupcat catalogue and add quantities

###packages###
import numpy as np
import h5py
import sys
import custom_utilities as cu
from astropy import cosmology
from astropy import table
from astropy.io import ascii

def main():
    ###make sure to change these when running in a new enviorment!###
    #location of data directory
    filepath1 = cu.get_output_path() + 'processed_data/NYU_VAGC/'
    savepath1 = filepath1 + 'custom_catalogues/'
    filepath2 = cu.get_output_path() + 'processed_data/mpa_dr7/'
    savepath2 = filepath2 + 'custom_catalogues/'
    filepath3 = cu.get_output_path() + 'processed_data/yang_groupcat/'
    savepath3 = filepath3 + 'custom_catalogues/'
    #################################################################

    cosmo = cosmology.FlatLambdaCDM(H0=100, Om0=0.3169) #set cosmology for all calculations, h=1

    catalogue1 = 'nyu_vagc_dr7'
    catalogue2 = 'gal_info_gal_totspecsfr_dr7_v5_2'
    try: sys.argv[1]
    except: catalogue3 = 'sample3_M_model'
    else: catalogue3 = sys.argv[1]
    print 'making', catalogue3, 'catalogue...'

    #open nyu vagc 
    print catalogue1
    f1 =  h5py.File(filepath1+catalogue1+'.hdf5', 'r')
    dset1 = f1.get(catalogue1)
    match13 = np.load(filepath1+'yang_groupcat_match/'+catalogue3+'_'+catalogue1+'_match.npy') #match into dset3

    #open mpa  catalogue
    print catalogue2
    f2 =  h5py.File(filepath2+catalogue2+'.hdf5', 'r')
    dset2 = f2.get(catalogue2)
    match23 = np.load(filepath2+'yang_groupcat_match/'+catalogue3+'_'+catalogue2+'_match.npy') #match into dset3

    #open groupcat
    print catalogue3
    f3 =  h5py.File(filepath3+catalogue3+'.hdf5', 'r')
    dset3 = f3.get(catalogue3)
    dset3 = np.array(dset3)
    match31 = np.load(filepath3+'nyu_vagc_match/'+catalogue1+'_'+catalogue3+'_match.npy') #match into dset1
    match32 = np.load(filepath3+'mpa_dr7_match/'+catalogue2+'_'+catalogue3+'_match.npy') #match into dset2
    
    dtype=[('ID','>i8'),('RA','>f8'),('DEC','>f8'),\
           ('Z','>f8'),('Z_ERR','>f8'),('Z_TYPE','>i8'),('VELDISP','>f8'),('VELDISP_ERR','>f8'),('FIBERCOL','>i8'),\
           ('M_u,0.1','>f8'),('M_g,0.1','>f8'),('M_r,0.1','>f8'),('M_i,0.1','>f8'),('M_z,0.1','>f8'),\
           ('N_SERSIC','>f8'),\
           ('MSTAR','>f8'),('SSFR','>f8'),\
           ('GROUP_ID','>i8'),('MGROUP','>f8'),('ZGROUP','>f8'),('R200','>f8'),('RPROJ','>f8'),('CEN_IND','>i8')]
    dtype = np.dtype(dtype)

    #create array to store catalogue in
    data = np.recarray((len(dset3),), dtype=dtype)
    data.fill(-99.9) #empty value indicator
    #input basics
    data['ID']  = dset3['IGAL']
    data['RA']  = dset3['RAgal']
    data['DEC'] = dset3['DECgal']
    data['Z']              = dset3['ZGAL'] 
    data['Z_ERR'][match13] = dset1['Z_ERR'][match31] #redshift err from sdss
    data['Z_TYPE']         = dset3['ZTYPE']
    data['VELDISP']     = dset3['VELDISP']
    data['VELDISP_ERR'] = dset3['VELDISP_ERR']
    #collions are where the nearest neighbor has been used
    collision = np.zeros((len(data),)).astype(int)
    collision[np.where(data['Z_TYPE']==4)] = 1
    data['FIBERCOL'] = collision    
    data['M_u,0.1'] = dset3['M_u,0.1'] #take magnitudes(model/petro) from group catalogue
    data['M_g,0.1'] = dset3['M_g,0.1']
    data['M_r,0.1'] = dset3['M_r,0.1']
    data['M_i,0.1'] = dset3['M_i,0.1']
    data['M_z,0.1'] = dset3['M_z,0.1']
    #add some derived quantities
    data['N_SERSIC'] = dset3['N_sersic']
    data['MSTAR'] = dset3['Mstar']
    data['SSFR'][match23] = dset2['MEDIAN'][match32]
    #add some group properties
    data['GROUP_ID'] = dset3['GROUP_ID']
    data['MGROUP'] = dset3['Mgroup']
    #data['ZGROUP'] = dset3['GROUP_Z']
    #identify central galaxy in groups
    group_IDs = np.unique(data['GROUP_ID'])
    
    for group in group_IDs: #run through each group and identify central galaxy
        members = np.where(data['GROUP_ID']==group)[0]
        largest_mass = np.min(data['MSTAR'][members]) #central is the brightest in m_r
        central = np.where((data['GROUP_ID']==group) & (data['MSTAR']==largest_mass))[0][0]
        data['CEN_IND'][members] = central
        data['ZGROUP'][members] = data['Z'][central] #group redshift is the central rtedshift
        #find the angular distance between the central galaxy and member galaxies
        da = cu.spheredist(data['RA'][central],data['DEC'][central],data['RA'][members],data['DEC'][members])
        da = np.radians(da) #convert to radians
        chi =  cosmology.funcs.comoving_distance(data['ZGROUP'][central], cosmo=cosmo)*1000.0 #in kpc
        dl = cosmology.funcs.luminosity_distance(data['ZGROUP'][central], cosmo=cosmo)*1000.0 #in kpc
        data['RPROJ'][members]=chi/(1.0+data['ZGROUP'][members])*da #caclulate physical seperation
    Omega_m = 0.3169
    x = 258.1 * (10**data['MGROUP']/(10.0**12.0))**(1.0/3.0)*(Omega_m/0.25)**(1.0/3.0)*(1.0+data['ZGROUP'])**(-1.0) 
    data['R200'] = x
    

    #replace any remaining empty flags with my empty flag
    print 'checking empty value flags...'
    for name in data.dtype.names:
        print name, len(np.where(data[name]==-99.0)[0])
        data[name][np.where(data[name]==-99.0)[0]]=-99.9

    print 'saving hdf5 version of the catalogue...'
    filename = catalogue3
    f = h5py.File(savepath3+filename+'.hdf5', 'w')
    dset = f.create_dataset(filename, data=data)
    f.close()

    print 'saving ascii version of the catalogue...'
    filename = catalogue3
    data_table = table.table.Table(data=data)
    ascii.write(data_table, savepath3+filename+'.dat')
    print data_table
    


if __name__ == '__main__':
  main()
