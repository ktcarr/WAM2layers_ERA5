# Moisture tracking for ERA5
This code is adapted from Ruud van der Ent's original implementation (https://github.com/ruudvdent/WAM2layersPython) and is based on the two-layer moisture tracking model described in the following paper:  

Van der Ent, R. J., L. Wang-Erlandsson, P. W. Keys, and H. H. Savenije, 2014: Contrasting roles of interception and transpiration in the hydrological cycle - Part 2: Moisture recycling. _Earth System Dynamics_, 5 (2), 471–489, https://doi.org/10.5194/esd-5-471-2014.

The version of the model in this repository was used for the paper "Impact of atmospheric circulation variability on U.S. Midwest moisture sources" by Carr and Ummenhofer (2023), accepted for publication in _Journal of Climate_ (https://doi.org/10.1175/JCLI-D-23-0178.1). This version of the code differs from the original implementation in that vertically-integrated moisture fluxes are estimated from pressure-level, rather than model-level, data. Those interested in applying the model should use the actively-maintained version, located here: https://github.com/WAM2layers/WAM2layers. 

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
