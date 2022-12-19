#!/bin/bash

## Arguments to the script:
## $1: year (int)

INPUT_FP=$1

for YEAR in {1979..2021..1}
    do
        sbatch merge-monthly-data.sh $INPUT_FP $YEAR
    done

