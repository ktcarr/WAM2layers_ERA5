#!/bin/bash

for YEAR in {1979..2020..1}
    do
        sbatch run.sh $YEAR 1 -1
    done

# Handle 2021 separately (only three months of data)
sbatch run.sh 2021 1 90

