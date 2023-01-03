# Moisture tracking for ERA5

## File structure  
```
|-- scripts  
     -- run.sh  
     -- run_parallel.sh  
    |-- move_to_scratch  
         -- move_to_scratch.sh  
         -- move_to_scratch_parallel.sh   
|-- src  
     -- __init__.py  
     -- get_fluxes.py  
     -- backtrack.py  
     -- postprocess.py  
     -- utils.py  
|-- data
     -- midwest_outline.csv  
|-- envs
 -- .gitignore  
 -- environment.yml  
 -- README.md  
 -- setup.py  
```

#### Description of files in ```./scripts/``` 
--```move_to_scratch/move_to_scratch.sh```: Move single year of ERA-5 data from file server to scratch directory.  
--```run.sh```: bash script to run model which contains parameters for given experiment.  
Scripts ending in ```_parallel.sh``` are used to parallelize across different years of data.


#### Description of files in ```./src/```
--```get_fluxes.py```: Compute 2-layer, vertically-integrated fluxes.  
--```backtrack.py```: Trace moisture back in time to source location.  
--```postprocess.py```: Compute tracked evaporation, and put data in nicer netcdf format.  
--```utils.py```: Compute constants (e.g., gridcell area)  


## Instructrions

#### Download data (1-degree lon/lat)
1. Download ERA-5 data:
    - $u,v,q$ (3-hourly resolution at all 37 pressure levels)    
    - $E,P$ (1-hourly resolution)
    - surface pressure, total column water (3-hourly)  
    - vertically-integrated water fluxes (3-hourly)    
2. Activate virtual environment with ```source activate ./envs```
3. Move data to scratch (machine-specific). Use ```./scripts/move_to_scratch_parallel.sh $INPUT_FP```.

## Running the model
1. Edit parameters as desired in the ```run.sh``` script.  
2. Navigate to the ```./scripts``` folder. Running options:  
-- Run the model for a single year with ```./run.sh ${YEAR} ${DOY_START} ${DOY_END}```.
-- Run the model for 1979-2021 with ```./run_parallel.sh```.
