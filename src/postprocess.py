import numpy as np
import xarray as xr


def data_path(year, doy_idx, fluxes_fp, tracked_moisture_fp, output_fp):
    load_Sa_track = os.path.join(
        tracked_moisture_fp, str(year) + "-" + str(doy_idx) + "Sa_track.mat"
    )

    load_fluxes_and_storages = os.path.join(
        fluxes_fp, str(year) + "-" + str(doy_idx) + "fluxes_storages.mat"
    )

    save_path = os.path.join(
        output_fp,
        "E_track_continental" + str(year) + ".mat",
    )

    return (
        load_Sa_track,
        load_fluxes_and_storages,
        save_path,
    )

def get_xr(loading_FS, loading_ST, latitude, longitude):
    """function for converting scipy.mat output to xarray"""

    ## get data 
    Sa_track_top = loading_ST["Sa_track_top"]
    Sa_track_down = loading_ST["Sa_track_down"] 
    north_loss = loading_ST["north_loss"]
    south_loss = loading_ST["south_loss"]
    down_to_top = loading_ST["down_to_top"]
    top_to_down = loading_ST["top_to_down"]
    water_lost = loading_ST["water_lost"]
    
    # load the total moisture data 
    Fa_E_top = loading_FS["Fa_E_top"]
    Fa_N_top = loading_FS["Fa_N_top"]
    Fa_E_down = loading_FS["Fa_E_down"]
    Fa_N_down = loading_FS["Fa_N_down"]
    Fa_Vert = loading_FS["Fa_Vert"]
    E = loading_FS["E"]
    P = loading_FS["P"]
    W_top = loading_FS["W_top"]
    W_down = loading_FS["W_down"]

    ## Stack data with two layers
    Sa_track = np.stack([Sa_track_down, Sa_track_top], axis=1)
    W = np.stack([W_down, W_top], axis=1)
    Fa_E = np.stack([Fa_E_down, Fa_E_top], axis=1)
    Fa_N = np.stack([Fa_N_down, Fa_N_top], axis=1)

    ## define coordinates and dimensions for xarray
    coords=dict(
                longitude=(["longitude"], longitude),
                latitude=(["latitude"], latitude),
                time=(["time"], np.arange(E.shape[0])),
                level=(["level"], ["down","top"])
    )

    dims_1d = ["time","longitude"]
    dims_2d = ["time","latitude","longitude"]
    dims_3d = ["time","level","latitude","longitude"]

    ## Convert to dataset
    storage = xr.Dataset(
        data_vars=dict(
                Sa_track=(dims_3d, Sa_track[1:]),
                W=(dims_3d, W[1:]),
                ),
        coords=coords
    )

    fluxes = xr.Dataset(
        data_vars=dict(
                Fa_E=(dims_3d, Fa_E),
                Fa_N=(dims_3d, Fa_N),
                Fa_Vert=(dims_2d, Fa_Vert),
                E=(dims_2d, E),
                P=(dims_2d, P),
                north_loss=(dims_1d, north_loss.squeeze()),
                south_loss=(dims_1d, south_loss.squeeze()),
                water_lost=(dims_2d, water_lost)
            ),
        coords=coords
    )

    return storage, fluxes



#%% Runtime & Results
if __name__ == "__main__":

    from timeit import default_timer as timer
    import argparse
    import os.path
    import src.utils
    import pandas as pd
    import scipy.io as sio

    #### Read parameters #####
    parser = argparse.ArgumentParser()

    ## Specify year and days of year for tracking
    parser.add_argument("--year", dest="year", type=int)
    parser.add_argument("--doy_start", dest="doy_start", type=int)
    parser.add_argument("--doy_end", dest="doy_end", type=int)

    ## Specify lon/lat grid for tracking
    parser.add_argument("--lon_min", dest="lon_min", type=float, default=0.0)
    parser.add_argument("--lon_max", dest="lon_max", type=float, default=359.0)
    parser.add_argument("--lat_min", dest="lat_min", type=float, default=-90.0)
    parser.add_argument("--lat_max", dest="lat_max", type=float, default=90.0)
    parser.add_argument("--dlat", dest="dlat", type=float, default=1.0)
    parser.add_argument("--dlon", dest="dlon", type=float, default=1.0)

    ## Data folders
    parser.add_argument("--fluxes_fp", dest="fluxes_fp", type=str)
    parser.add_argument("--tracked_moisture_fp", dest="tracked_moisture_fp", type=str)
    parser.add_argument("--output_fp", dest="output_fp", type=str)

    args = parser.parse_args()

    # Get lon/lat, and gridcell dimensions
    constants = src.utils.get_constants_from_args(args)
    longitude = constants["longitude"]
    latitude = constants["latitude"]

    doy_indices = src.utils.get_doy_indices(args.doy_start, args.doy_end, args.year)

    ##### Begin the actual postprocessing ####
    start1 = timer()

    fluxes_daily = []
    storage_daily = []

    for i, doy_idx in enumerate(doy_indices):
        ## i is the loop iteration, doy_idx is the index of the DOY

        start = timer()
        datapath = data_path(year=args.year, doy_idx=doy_idx, fluxes_fp=args.fluxes_fp, tracked_moisture_fp=args.tracked_moisture_fp, output_fp=args.output_fp)

        # load tracked data
        loading_ST = sio.loadmat(datapath[0], verify_compressed_data_integrity=False)
        loading_FS = sio.loadmat(datapath[1], verify_compressed_data_integrity=False) 
        storage, fluxes = get_xr(loading_ST=loading_ST, loading_FS=loading_FS, latitude=latitude, longitude=longitude)

        # compute tracked evaporation
        tracked_to_total_frac = storage["Sa_track"] / storage["W"]
        fluxes["E_track"] = fluxes["E"] * tracked_to_total_frac.sel(level="down")

        # compute daily avg/totals
        fluxes_daily.append(fluxes.sum("time"))
        storage_daily.append(storage.mean("time"))

        end = timer()
        print(f"Postprocess runtime for day {doy_idx+1}: {end-start:.2f} seconds")
   
    ## concatenate data from individual days
    time_idx = src.utils.get_time_idx(doy_indices, args.year)
    fluxes_daily = xr.concat(fluxes_daily, dim=time_idx)
    storage_daily= xr.concat(storage_daily,dim=time_idx)

    ## save to file
    fluxes_daily.to_netcdf(os.path.join(args.output_fp, f"fluxes_daily_{args.year}.nc"))
    storage_daily.to_netcdf(os.path.join(args.output_fp,f"storage_daily_{args.year}.nc"))

    end1 = timer()
    print(f"The total runtime of Con_E_Recyc_Output is {end1 - start1:.2f} seconds.")
