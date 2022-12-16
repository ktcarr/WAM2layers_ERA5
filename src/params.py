####################################### IMPORTS #########################################
import numpy as np
from os.path import join

# ## Filepaths for moisture tracking model
# data_fp   = "/vortexfs1/scratch/kcarr/moisture_track"
# input_fp  = join(data_fp,  "input")
# output_fp = join(data_fp,  "output")
# inter_fp  = join(data_fp,  "interdata")
# lsm_fp    = join(input_fp, "lsm.nc")

# File paths
data_fp   = '/vortexfs1/scratch/kcarr/'
input_fp  = join(data_fp, 'moisture-track_input/')
output_fp = join(data_fp, 'moisture-track_output/divt45/output/')
inter_fp  = join(data_fp, 'moisture-track_output/divt45/interdata/')
lsm_fp    = join(input_folder, 'lsm.nc')

####################################### PARAMETERS  for tracking ######################################
years     = np.arange(1979, 2021)
longitude = np.arange(0,360)     
latitude  = np.arange(80,-80,-1) #latitude must be in descending order!! 
dlat      = 1. # horizontal resolution of data
dlon      = 1.

### Parameters for fluxes/tracking
boundary   = 29 # vertical separation is at ~816 hPa for surface pressure = 1013 hPa
divt       = 45 # division of the timestep, 30 means a calculation timestep of 180 min/30 min = 6 min (match Imme)
count_time = 8 # number of time-steps to process at one time

# Other parameters
isglobal = 1 # fill in 1 for global computations (i.e. Earth round), fill in 0 for a local domain with boundaries

