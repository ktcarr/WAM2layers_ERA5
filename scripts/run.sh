#!/bin/bash
#SBATCH --partition=compute
#SBATCH --mail-type=NONE
#SBATCH --ntasks=1
#SBATCH --mem=64gb
#SBATCH --time=23:59:00
#SBATCH --job-name global
#SBATCH --output=log_%j.log


###### Parameters #######

## Set filepaths ##
LOCAL_FP="/vortexfs1/home/kcarr/WAM2layers_ERA5"
DATA_FP="$LOCAL_FP/data"

SCRATCH_FP="/vortexfs1/scratch/kcarr/WAM2layer-data"
INPUT_FP=$SCRATCH_FP/input
FLUXES_FP=$SCRATCH_FP/interdata
TRACKED_MOISTURE_FP=$SCRATCH_FP/tracked_moisture
OUTPUT_FP=$SCRATCH_FP/output

## Outline of region from which to track moisture
REGION_FP=$DATA_FP/midwest_outline.csv

## Grid specific constants (lon/lat range and horiz. resolution) ##
LON_MIN="   0" 
LON_MAX=" 359"
LAT_MIN=" -80"
LAT_MAX="  80"
DLAT=1
DLON=1

## Numerical parameters ##
TIMESTEP=10800 # units: seconds
DIVT=45        # interpolation factor (for time)
KVF=2          # vertical diffusivity
COUNT_TIME=8   # number of timesteps to process at once
BOUNDARY=29    # index of pressure level which divides model layers
IS_GLOBAL=1    # does the eastern boundary touch the western boundary?

module load mambaforge
source activate $LOCAL_FP/envs

python -u $LOCAL_FP/src/get_fluxes.py --year $1 \
                                      --lon_min $LON_MIN \
                                      --lon_max $LON_MAX \
                                      --lat_min $LAT_MIN \
                                      --lat_max $LAT_MAX \
                                      --dlat $DLAT \
                                      --dlon $DLON \
                                      --timestep $TIMESTEP \
                                      --divt $DIVT \
                                      --boundary $BOUNDARY \
                                      --count_time $COUNT_TIME \
                                      --is_global $IS_GLOBAL \
                                      --input_fp $INPUT_FP \
                                      --fluxes_fp $FLUXES_FP
                                      # --doy_start $2 \
                                      # --doy_end $3 \

python -u $LOCAL_FP/src/backtrack.py --year $1 \
                                     --lon_min $LON_MIN \
                                     --lon_max $LON_MAX \
                                     --lat_min $LAT_MIN \
                                     --lat_max $LAT_MAX \
                                     --dlat $DLAT \
                                     --dlon $DLON \
                                     --divt $DIVT \
                                     --count_time $COUNT_TIME \
                                     --is_global $IS_GLOBAL \
                                     --input_fp $INPUT_FP \
                                     --fluxes_fp $FLUXES_FP
                                     --tracked_moisture_fp $TRACKED_MOISTURE_FP
                                     # --doy_start $2 \
                                     # --doy_end $3 \
