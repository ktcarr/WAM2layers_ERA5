# -*- coding: utf-8 -*-
"""
Created on Tue Jun 21 12:38:45 2016

@author: Ent00002
"""
import numpy as np
import xarray as xr

### Physical constants ###



def get_longitude(lon_min, lon_max, dlon):
    '''Get longitude'''
    return np.arange(lon_min, lon_max+1, dlon)

def get_latitude(lat_min, lat_max, dlat):
    '''Get latitude (in descending order)'''
    latitude  = np.arange(lat_min, lat_max+1, dlat)
    return latitude[::-1]

def get_gridcell_dims(latitude, dlon, dlat): 
    '''Get longitude/latitude range, and gridcell dimensions,
    based on user inputs'''

    ## constants
    dg = 111089.56 # [m] length of 1 degree latitude
    Erad = 6.371e6 # [m] Earth radius 
     
    # new area size calculation:
    lat_n_bound = np.minimum(90.0,  latitude + 0.5*dlat)
    lat_s_bound = np.maximum(-90.0, latitude - 0.5*dlat)
    
    A_gridcell      = np.zeros([len(latitude),1])
    A_gridcell[:,0] = (np.pi/180.0) * Erad**2 * dlat
                       abs( np.sin(lat_s_bound*np.pi/180.0) -\
                            np.sin(lat_n_bound*np.pi/180.0) )
    
    # Get length of N/S gridcell boundaries (in meters)
    L_N_gridcell = dlon * np.cos((latitude + dlat / 2.0) * np.pi / 180.0) * dg
    L_S_gridcell = dlon * np.cos((latitude - dlat / 2.0) * np.pi / 180.0) * dg

    # Get length of east/west bondary
    L_EW_gridcell = dlat * dg

    return A_gridcell, L_N_gridcell, L_S_gridcell, L_EW_gridcell

def getconstants(lon_min, lon_max, dlon, lat_min, lat_max, dlat):
    '''Wrapper function to get constants based on specified inputs.
    Returns dictionary with constants'''
    latitude = get_latitude(lat_min, lat_max, dlat)
    longitude = get_longitude(lon_min, lon_max, dlon)
    A_gridcell,L_N_gridcell,L_S_gridcell,L_EW_gridcell = get_gridcell_dims(
                            latitude=latitude,
                            dlon=dlon,
                            dlat=dlat)

    ## Put results in dictionary
    constants = {'g': 9.80665, # [m/s2] from ERA-interim archive
                 'density_water' : 1000, # [kg/m3]
                 'longitude':longitude,
                 'latitude':latitude,
                 'A_gridcell':A_gridcell,
                 'L_N_gridcell':L_N_gridcell,
                 'L_S_gridcell':L_S_gridcell,
                 'L_EW_gridcell':L_EW_gridcell}
    
    return constants


# def get_lsm(lsm_path):
#     '''Tentative function to get LSM'''
#     # Create land-sea-mask (in this model lakes are considered part of the land)
#     # 0 = sea, 1 = land
#     lsm = xr.open_dataset(lsm_path)['lsm'].isel(time=0, drop=True)
#     lsm = lsm.sel(latitude=latitude, longitude=longitude) # determine lon/lat values to use
#     lsm       = lsm.values
#         
#     lsm[0,:] = 0 # the northern boundary is always oceanic = 0
#     lsm[-1,:] = 0 # the southern boundary is always oceanic = 0