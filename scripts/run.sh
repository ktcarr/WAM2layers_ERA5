#!/bin/bash
#SBATCH --partition=compute
#SBATCH --mail-type=NONE
#SBATCH --ntasks=1
#SBATCH --mem=64gb
#SBATCH --time=23:59:00
#SBATCH --job-name track
#SBATCH --output=log_%j.log

## Parse arguments to script
YEAR=$1
DOY_START=$2
DOY_END=$3

## Set filepaths ##
LOCAL_FP="/vortexfs1/home/kcarr/WAM2layers_ERA5"
DATA_FP="$LOCAL_FP/data"

SCRATCH_FP="/vortexfs1/scratch/kcarr/WAM2layer-data"
INPUT_FP=$SCRATCH_FP/input
FLUXES_FP=$SCRATCH_FP/fluxes
TRACKED_MOISTURE_FP=$SCRATCH_FP/tracked_moisture
OUTPUT_FP=$SCRATCH_FP/output

## Outline of region from which to track moisture
REGION_FP=$DATA_FP/midwest_outline.csv

## Grid specific constants (lon/lat range and horiz. resolution) ##
LON_MIN="   0" 
LON_MAX=" 359"
LAT_MIN=" -30"
LAT_MAX="  80"
DLAT=1
DLON=1

## Numerical parameters ##
KVF=2          # vertical diffusivity
FREQ=8         # frequency of data (timesteps per day)
FREQ_EP=24     # frequency of E and P data (timesteps per day)
DIVT=30        # time interpolation factor (increase FREQ by this factor)
BOUNDARY=29    # index of pressure level which divides model layers
IS_GLOBAL=1    # does the eastern boundary touch the western boundary?

# load environment
# module load mambaforge
# source activate $LOCAL_FP/envs

# run the model
python -u $LOCAL_FP/src/get_fluxes.py --year $YEAR \
                                      --lon_min $LON_MIN \
                                      --lon_max $LON_MAX \
                                      --lat_min $LAT_MIN \
                                      --lat_max $LAT_MAX \
                                      --dlat $DLAT \
                                      --dlon $DLON \
                                      --divt $DIVT \
                                      --boundary $BOUNDARY \
                                      --freq $FREQ \
                                      --freq_ep $FREQ_EP \
                                      --is_global $IS_GLOBAL \
                                      --input_fp $INPUT_FP \
                                      --fluxes_fp $FLUXES_FP \
                                      --doy_start $DOY_START \
                                      --doy_end $DOY_END

# python -u $LOCAL_FP/src/backtrack.py --year $YEAR \
#                                      --lon_min $LON_MIN \
#                                      --lon_max $LON_MAX \
#                                      --lat_min $LAT_MIN \
#                                      --lat_max $LAT_MAX \
#                                      --dlat $DLAT \
#                                      --dlon $DLON \
#                                      --divt $DIVT \
#                                      --kvf $KVF \
#                                      --freq $FREQ \
#                                      --is_global $IS_GLOBAL \
#                                      --fluxes_fp $FLUXES_FP \
#                                      --input_fp $INPUT_FP \
#                                      --tracked_moisture_fp $TRACKED_MOISTURE_FP \
#                                      --region_fp $REGION_FP \
#                                      --doy_start $DOY_START \
#                                      --doy_end $DOY_END
# 
# python -u $LOCAL_FP/src/postprocess.py --year $YEAR \
#                                        --lon_min $LON_MIN \
#                                        --lon_max $LON_MAX \
#                                        --lat_min $LAT_MIN \
#                                        --lat_max $LAT_MAX \
#                                        --dlat $DLAT \
#                                        --dlon $DLON \
#                                        --fluxes_fp $FLUXES_FP \
#                                        --tracked_moisture_fp $TRACKED_MOISTURE_FP \
#                                        --output_fp $OUTPUT_FP \
#                                        --doy_start $DOY_START \
#                                        --doy_end $DOY_END
