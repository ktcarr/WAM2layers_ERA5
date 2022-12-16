# -*- coding: utf-8 -*-
"""
Created on Tue Jun 21 12:38:45 2016

@author: Ent00002
"""
import numpy as np
import xarray as xr

def getconstants(latitude,longitude,lake_mask,lsm_path): # def getconstants in Python is the same as function in MATLAB. 
   
    # Create land-sea-mask (in this model lakes are considered part of the land)
    # 0 = sea, 1 = land
    lsm = xr.open_dataset(lsm_path)['lsm'].isel(time=0, drop=True)
    lsm = lsm.sel(latitude=latitude, longitude=longitude) # determine lon/lat values to use
    lsm       = lsm.values
        
    lsm[0,:] = 0 # the northern boundary is always oceanic = 0
    lsm[-1,:] = 0 # the southern boundary is always oceanic = 0
    
    # Constants 
    g = 9.80665 # [m/s2] from ERA-interim archive
    density_water = 1000 # [kg/m3]
    dg = 111089.56 # [m] length of 1 degree latitude
    timestep = 3*3600 # [s] timestep of three hours (= 3*3600 seconds) for all data except e and p
    Erad = 6.371e6 # [m] Earth radius
    
    # Semiconstants
    gridcell = np.abs(longitude[1] - longitude[0]) # [degrees] grid cell size
    
    # new area size calculation:
    lat_n_bound = np.minimum(90.0 , latitude + 0.5*gridcell)
    lat_s_bound = np.maximum(-90.0 , latitude - 0.5*gridcell)
    
    A_gridcell = np.zeros([len(latitude),1])
    A_gridcell[:,0] = (np.pi/180.0)*Erad**2 * abs( np.sin(lat_s_bound*np.pi/180.0) - np.sin(lat_n_bound*np.pi/180.0) ) * gridcell
    
    L_N_gridcell = gridcell * np.cos((latitude + gridcell / 2.0) * np.pi / 180.0) * dg # [m] length northern boundary of a cell
    L_S_gridcell = gridcell * np.cos((latitude - gridcell / 2.0) * np.pi / 180.0) * dg # [m] length southern boundary of a cell
    L_EW_gridcell = gridcell * dg # [m] length eastern/western boundary of a cell 
    
    return lsm, g, density_water, timestep,  A_gridcell, L_N_gridcell, L_S_gridcell, L_EW_gridcell, gridcell
