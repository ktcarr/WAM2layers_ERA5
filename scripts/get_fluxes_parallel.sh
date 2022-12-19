#!/bin/bash

# first argument to the script is the path of the parameters folder

for i in {1979..2020..1}
    do
        sbatch get-fluxes.sh $1 $i
    done

