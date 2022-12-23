import numpy as np
import scipy.io as sio
import datetime
import os


def data_path(year, doy_idx, fluxes_fp, tracked_moisture_fp):
    load_Sa_track = os.path.join(
        tracked_moisture_fp, str(year) + "-" + str(doy_idx) + "Sa_track.mat"
    )

    load_fluxes_and_storages = os.path.join(
        fluxes_fp, str(year) + "-" + str(doy_idx) + "fluxes_storages.mat"
    )

    save_path = os.path.join(
        output_folder,
        "E_track_continental" + str(year) + ".mat",
    )

    return (
        load_Sa_track,
        load_fluxes_and_storages,
        save_path,
    )


#%% Runtime & Results
if __name__ == "__main__":

    from timeit import default_timer as timer
    import argparse

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
    constants = utils.get_constants_from_args(args)

    # Parse constants
    g = constants["g"]
    density_water = constants["density_water"]
    longitude = constants["longitude"]
    latitude = constants["latitude"]
    A_gridcell = constants["A_gridcell"]
    L_N_gridcell = constants["L_N_gridcell"]
    L_S_gridcell = constants["L_S_gridcell"]
    L_EW_gridcell = constants["L_EW_gridcell"]

    doy_indices = utils.get_doy_indices(args.doy_start, args.doy_end, args.year)

    ##### Begin the actual postprocessing ####
    start1 = timer()

    N = len(doy_indices)  # number of days to process

    E_per_day = np.zeros((N, len(latitude), len(longitude)))
    E_track_per_day = np.zeros((N, len(latitude), len(longitude)))
    P_per_day = np.zeros((N, len(latitude), len(longitude)))
    Sa_track_down_per_day = np.zeros((N, len(latitude), len(longitude)))
    Sa_track_top_per_day = np.zeros((N, len(latitude), len(longitude)))
    W_down_per_day = np.zeros((N, len(latitude), len(longitude)))
    W_top_per_day = np.zeros((N, len(latitude), len(longitude)))
    north_loss_per_day = np.zeros((N, 1, len(longitude)))
    south_loss_per_day = np.zeros((N, 1, len(longitude)))
    down_to_top_per_day = np.zeros((N, len(latitude), len(longitude)))
    top_to_down_per_day = np.zeros((N, len(latitude), len(longitude)))
    water_lost_per_day = np.zeros((N, len(latitude), len(longitude)))

    for i, doy_idx in enumerate(doy_indices):
        ## i is the loop iteration, doy_idx is the index of the DOY

        start = timer()
        datapath = data_path(y, doy_idx, years, timetracking)

        # load tracked data
        loading_ST = sio.loadmat(datapath[0], verify_compressed_data_integrity=False)
        Sa_track_top = loading_ST["Sa_track_top"]
        Sa_track_down = loading_ST["Sa_track_down"]
        north_loss = loading_ST["north_loss"]
        south_loss = loading_ST["south_loss"]
        down_to_top = loading_ST["down_to_top"]
        top_to_down = loading_ST["top_to_down"]
        water_lost = loading_ST["water_lost"]
        Sa_track = Sa_track_top + Sa_track_down

        # load the total moisture data
        loading_FS = sio.loadmat(datapath[1], verify_compressed_data_integrity=False)
        Fa_E_top = loading_FS["Fa_E_top"]
        Fa_N_top = loading_FS["Fa_N_top"]
        Fa_E_down = loading_FS["Fa_E_down"]
        Fa_N_down = loading_FS["Fa_N_down"]
        Fa_Vert = loading_FS["Fa_Vert"]
        E = loading_FS["E"]
        P = loading_FS["P"]
        W_top = loading_FS["W_top"]
        W_down = loading_FS["W_down"]

        W = W_top + W_down

        # compute tracked evaporation
        E_track = E[:, :, :] * (Sa_track_down[1:, :, :] / W_down[1:, :, :])

        # save per day
        E_per_day[i, :, :] = np.sum(E, axis=0)
        E_track_per_day[i, :, :] = np.sum(E_track, axis=0)
        P_per_day[i, :, :] = np.sum(P, axis=0)
        Sa_track_down_per_day[i, :, :] = np.mean(Sa_track_down[1:, :, :], axis=0)
        Sa_track_top_per_day[i, :, :] = np.mean(Sa_track_top[1:, :, :], axis=0)
        W_down_per_day[i, :, :] = np.mean(W_down[1:, :, :], axis=0)
        W_top_per_day[i, :, :] = np.mean(W_top[1:, :, :], axis=0)

        north_loss_per_day[i, :, :] = np.sum(north_loss, axis=0)
        south_loss_per_day[i, :, :] = np.sum(south_loss, axis=0)
        down_to_top_per_day[i, :, :] = np.sum(down_to_top, axis=0)
        top_to_down_per_day[i, :, :] = np.sum(top_to_down, axis=0)
        water_lost_per_day[i, :, :] = np.sum(water_lost, axis=0)

        end = timer()
        print(
            "Runtime output for day " + str(doy_idx + 1) + " in year " + str(y) + " is",
            (end - start),
            " seconds.",
        )

    if daily == 1:
        if timetracking == 0:  # create dummy values
            Sa_time_down_per_day = 0
            Sa_time_top_per_day = 0
            E_time_per_day = 0

        sio.savemat(
            datapath[2],
            {
                "E_per_day": E_per_day,
                "E_track_per_day": E_track_per_day,
                "P_per_day": P_per_day,
                "Sa_track_down_per_day": Sa_track_down_per_day,
                "Sa_track_top_per_day": Sa_track_top_per_day,
                "W_down_per_day": W_down_per_day,
                "W_top_per_day": W_top_per_day,
                "E_time_per_day": E_time_per_day,
                "water_lost_per_day": water_lost_per_day,
            },
            do_compression=True,
        )

    end1 = timer()
    print("The total runtime of Con_E_Recyc_Output is", (end1 - start1), " seconds.")
