# Moisture tracking for ERA5

## File structure  
```
|-- scripts  
     -- backtrack_parallel.sh  
     -- backtrack.sh  
     -- get_fluxes_parallel.sh   
     -- get_fluxes.sh  
    |-- move_to_scratch  
         -- move_to_scratch.sh  
         -- move_to_scratch_parallel.sh  
|-- src  
     -- __init__.py  
     -- backtrack.py  
     -- getconstants.py  
     -- get_fluxes.py  
     -- params.py  
     -- postprocess.py  
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


## Description of files in ```./src/```

```get_fluxes.py```: Compute 2-layer, vertically-integrated fluxes.  

```backtrack.py```: Trace moisture back in time to source location.

```postprocess.py```: Compute tracked evaporation, and put data in nicer netcdf format. 

```getconstants.py```: Contains constants (e.g., value of $g$).
   
```params.py```: parameters for given experiment, including lon/lat ranges, and directories to save to.

## Running the model
1. Run `get_fluxes_new.sh' to compute fluxes with specified timestep value (set by divt in parameter file). Parameter files are found in `run/params/midwest_global' folder.  

2. Run `backtrack_and_postprocess' to track moisture back in time and postprocess output. Need to specify flux parameter file and kvf. Note: this version of the script processes each year separately, meaning that rainfall in jan/feb may not be accurately tracked.  

3. Run `python /vortexfs1/home/kcarr/iap-2021/aggregate_output.py' to put .mat files output from moisture tracking into a nicer NetCDF4 format.  
