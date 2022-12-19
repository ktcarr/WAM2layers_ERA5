#!/bin/bash
#SBATCH --partition=compute
#SBATCH --mail-type=NONE
#SBATCH --ntasks=1
#SBATCH --mem=64gb
#SBATCH --time=23:59:00
#SBATCH --job-name global
#SBATCH --output=log_%j.log

# Arguments are:
# $1 params_fp (specifies divt and filepaths)
# $2 Kvf       (vertical dispersion factor)
# $3 year

RUN_FP="/vortexfs1/home/kcarr/moisture_track3.0/run"
VERTICES_FP="$RUN_FP/params/midwest_global/midwest_vertices.npy"

module load anaconda
source activate torch_env

echo "Begin backtracking"
python -u $RUN_FP/backtrack.py --params_fp $1 \
                               --Kvf $2 \
                               --year1 $3 \
                               --year2 $3 \
                               --veryfirstrun 1 \
                               --region_fp $VERTICES_FP \
                               --list_of_days_fp "none" 

echo "\nBegin postprocessing"
python -u $RUN_FP/postprocess.py --params_fp $1 \
                                 --Kvf $2 \
                                 --year $3 \
                                 --doy_end 365 
