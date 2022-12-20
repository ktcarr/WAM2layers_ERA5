#!/bin/bash
#SBATCH --partition=compute
#SBATCH --ntasks=1
#SBATCH --mem=4gb
#SBATCH --time=04:00:00
#SBATCH --job-name=move-data
#SBATCH --output=log-%j.log

module load cdo

## Arguments to the script:
## $1: year (int)
## $2: directory in which to save the data (string)

INPUT_FP=$1
YEAR=$2
ERA_FP="/vortexfs1/share/cmip6/data/era5/reanalysis"

echo "Merging single-level fields"
ERA_3HR_FP="$ERA_FP/single-levels/3hr"
for VARNAME in vert-int sp tcw
    do
        cdo -b F64 -mergetime ${ERA_3HR_FP}/${VARNAME}/${VARNAME}_${YEAR}-*.nc ${INPUT_FP}/${YEAR}-${VARNAME}.nc
    done

echo "Getting E and P"
ERA_1HR_FP="$ERA_FP/single-levels/1hr"
for VARNAME in evaporation total_precipitation
    do
        cdo -b F64 -mergetime ${ERA_1HR_FP}/${VARNAME}/${YEAR}-*_${VARNAME}.nc ${INPUT_FP}/${YEAR}-${VARNAME}_temp.nc
        cdo remapbil,r360x181 ${INPUT_FP}/${YEAR}-${VARNAME}_temp.nc ${INPUT_FP}/${YEAR}-${VARNAME}.nc
        rm ${INPUT_FP}/${YEAR}-${VARNAME}_temp.nc
    done

echo "Merging pressure-level fields..."
UVQ_FP="$ERA_FP/pressure-levels/3hr/uvq"
cdo -b F64 -mergetime $UVQ_FP/uvq-$YEAR-*.nc $INPUT_FP/$YEAR-uvq.nc
echo "Done."
