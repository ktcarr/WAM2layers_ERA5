#!/bin/bash

for YEAR in {1980..2020..1}
    do
        sbatch run.sh $YEAR 1 -1
    done

# Handle 1979 separately (first 7 hours of E & P data missing)
sbatch run.sh 1979 2 365
    
# Handle 2021 separately (only three months of data)
sbatch run.sh 2021 1 89

