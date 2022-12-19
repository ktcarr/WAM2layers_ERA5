# Moisture tracking for ERA5

## File structure  
```
|-- scripts  
     -- get_fluxes.sh  
     -- backtrack.sh  
     -- backtrack_parallel.sh  
     -- get_fluxes_parallel.sh   
    |-- move_to_scratch  
         -- move_to_scratch.sh  
         -- move_to_scratch_parallel.sh  
|-- src  
     -- __init__.py  
     -- get_fluxes.py  
     -- backtrack.py  
     -- postprocess.py  
     -- utils.py    
|-- tests  
 -- .gitignore  
 -- environment.yml  
 -- README.md  
 -- setup.py  
```

## Instructrions

#### Download data (1-degree lon/lat)
0. Activate virtual environment with ```source activate ./envs```
1. Download ERA-5 data:
    - $u,v,q$ (3-hourly resolution at all 37 pressure levels)    
    - $E,P$ (1-hourly resolution)
    - surface pressure, total column water (3-hourly)  
    - vertically-integrated water fluxes (3-hourly)    
2. Move data to scratch (machine-specific). Use ```./scripts/move_to_scratch_parallel.sh $INPUT_FP```.


## Description of files in ```./scripts/``` 

```move-data.sh```: Move single year of ERA-5 data from file server to local directory.  

```move_to_scratch.sh```: Move single year of ERA-5 data from file server to scratch directory.  

```move_to_scratch_parallel.sh```: Move all years in parallel.  

```get_fluxes.sh```: bash script to compute fluxes. Contains parameters for given experiment; this inclues lon/lat ranges, filepaths, and numerical parameters.  


## Description of files in ```./src/```

```get_fluxes.py```: Compute 2-layer, vertically-integrated fluxes.  

```backtrack.py```: Trace moisture back in time to source location.

```postprocess.py```: Compute tracked evaporation, and put data in nicer netcdf format. 

```utils.py```: Compute constants (e.g., gridcell area)  

## Running the model
1. Edit parameters as desired in the ```get_fluxes.sh``` script. After editing, run a single year with ```./get_fluxes --year <YEAR>```  
