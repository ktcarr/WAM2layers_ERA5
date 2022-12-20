#!/bin/bash
#SBATCH --partition=compute
#SBATCH --mail-type=NONE
#SBATCH --ntasks=1
#SBATCH --mem=64gb
#SBATCH --time=23:59:00
#SBATCH --job-name global
#SBATCH --output=log_%j.log

LOCAL_FP="/vortexfs1/home/kcarr/WAM2layers_ERA5"
DATA_FP="/vortexfs1/scratch/kcarr/WAM2layer-data"

INPUT_FP=$DATA_FP/input
FLUXES_FP=$DATA_FP/interdata
TRACKED_MOISTURE_FP=$DATA_FP/tracked_moisture
OUTPUT_FP=$DATA_FP/output

module load mambaforge
source activate $LOCAL_FP/envs

python -u $LOCAL_FP/src/get_fluxes.py --year $1 \
                                      --lon_min "0" \
                                      --lon_max "359" \
                                      --lat_min  "-80" \
                                      --lat_max " 80" \
                                      --dlat 1 \
                                      --dlon 1 \
                                      --timestep "10800" \
                                      --divt "45" \
                                      --boundary "29" \
                                      --count_time "8" \
                                      --is_global "1" \
                                      --input_folder $INPUT_FP \
                                      --interdata_folder $FLUXES_FP
                                      # --doy_start $2 \
                                      # --doy_end $3 \
