#!/bin/bash
#SBATCH --partition=compute
#SBATCH --mail-type=NONE
#SBATCH --ntasks=1
#SBATCH --mem=64gb
#SBATCH --time=23:59:00
#SBATCH --job-name global
#SBATCH --output=log_%j.log

RUN_FP="/vortexfs1/home/kcarr/moisture_track3.0/run"
# PARAMS_FP="$RUN_FP/params/midwest_new"

module load anaconda
source activate torch_env
python -u $RUN_FP/get_fluxes.py --params_fp $1 \
                                --year1 $2 \
                                --year2 $2
