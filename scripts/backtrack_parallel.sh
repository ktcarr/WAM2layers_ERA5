#!/bin/bash

# Arguments are:
# $1 params_fp (specifies divt and filepaths)
# $2 Kvf       (vertical dispersion factor)

for i in {1979..2020..1}
    do
        sbatch backtrack.sh $1 $2 $i
    done

