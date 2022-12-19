import numpy as np
import scipy.io as sio
import calendar
import datetime
from getconstants import getconstants
from timeit import default_timer as timer
import os
import importlib.util
import argparse


# END OF INPUT
#%% Datapaths (FILL THIS IN)


def data_path(y, a, years, timetracking):
    load_Sa_track = os.path.join(
        sub_interdata_folder, str(y) + "-" + str(a) + "Sa_track.mat"
    )

    load_Sa_time = os.path.join(
        sub_interdata_folder, str(y) + "-" + str(a) + "Sa_time.mat"
    )

    load_fluxes_and_storages = os.path.join(
        interdata_folder, str(y) + "-" + str(a) + "fluxes_storages.mat"
    )

    save_path = os.path.join(
        output_folder,
        "E_track_continental_full"
        + str(years[0])
        + "-"
        + str(years[-1])
        + "-timetracking"
        + str(timetracking)
        + ".mat",
    )

    save_path_daily = os.path.join(
        output_folder,
        "E_track_continental_daily_full"
        + str(y)
        + "-timetracking"
        + str(timetracking)
        + ".mat",
    )

    return (
        load_Sa_track,
        load_Sa_time,
        load_fluxes_and_storages,
        save_path,
        save_path_daily,
    )


#%% Runtime & Results
if __name__ == "__main__":
    #### Read parameters #####
    parser = argparse.ArgumentParser()
    parser.add_argument("--params_fp", dest="params_path")
    parser.add_argument("--year", dest="year", type=int)
    parser.add_argument("--doy_end", dest="doy_end", type=int)
    parser.add_argument("--Kvf", dest="Kvf")
    args = parser.parse_args()

    params_path = args.params_path
    years = [args.year]  # years plural for backwards compatibility with old version..
    doy_end = args.doy_end
    Kvf = args.Kvf

    ###  Use this code to import specified params.py file ##############
    ####################################################################
    spec = importlib.util.spec_from_file_location("module.name", params_path)
    params = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(params)

    ######### Parameters (shared for all experiments) ###################
    latitude = params.latitude
    longitude = params.longitude
    lsm_path = params.lsm_path
    interdata_folder = params.interdata_folder
    output_folder = os.path.join(params.output_folder, f"kvf_{int(Kvf)}")
    sub_interdata_folder = os.path.join(output_folder, "back")
    lake_mask = np.nan  # dummy argument for compatibility with other stuff...

    #%%BEGIN OF INPUT1 (FILL THIS IN)
    yearpart = np.arange(doy_end)
    daily = 1  # 1 for writing out daily data, 0 for only monthly data
    timetracking = 0  # 0 for not tracking time and 1 for tracking time

    # obtain the constants
    (
        lsm,
        g,
        density_water,
        timestep,
        A_gridcell,
        L_N_gridcell,
        L_S_gridcell,
        L_EW_gridcell,
        gridcell,
    ) = getconstants(latitude, longitude, lake_mask, lsm_path)

    start1 = timer()
    startyear = years[0]

    E_per_year_per_month = np.zeros((len(years), 12, len(latitude), len(longitude)))
    E_track_per_year_per_month = np.zeros(
        (len(years), 12, len(latitude), len(longitude))
    )
    P_per_year_per_month = np.zeros((len(years), 12, len(latitude), len(longitude)))
    Sa_track_down_per_year_per_month = np.zeros(
        (len(years), 12, len(latitude), len(longitude))
    )
    Sa_track_top_per_year_per_month = np.zeros(
        (len(years), 12, len(latitude), len(longitude))
    )
    W_down_per_year_per_month = np.zeros(
        (len(years), 12, len(latitude), len(longitude))
    )
    W_top_per_year_per_month = np.zeros((len(years), 12, len(latitude), len(longitude)))
    north_loss_per_year_per_month = np.zeros((len(years), 12, 1, len(longitude)))
    south_loss_per_year_per_month = np.zeros((len(years), 12, 1, len(longitude)))
    down_to_top_per_year_per_month = np.zeros(
        (len(years), 12, len(latitude), len(longitude))
    )
    top_to_down_per_year_per_month = np.zeros(
        (len(years), 12, len(latitude), len(longitude))
    )
    water_lost_per_year_per_month = np.zeros(
        (len(years), 12, len(latitude), len(longitude))
    )

    if timetracking == 1:
        Sa_time_down_per_year_per_month = np.zeros(
            (len(years), 12, len(latitude), len(longitude))
        )
        Sa_time_top_per_year_per_month = np.zeros(
            (len(years), 12, len(latitude), len(longitude))
        )
        E_time_per_year_per_month = np.zeros(
            (len(years), 12, len(latitude), len(longitude))
        )

    for y in years:
        ly = int(calendar.isleap(y))
        final_time = 364 + ly
        # this_yearpart = np.arange(90) if y==2021 else yearpart
        this_yearpart = np.arange(91)
        # this_yearpart = np.arange(183,214)

        E_per_day = np.zeros((365 + ly, len(latitude), len(longitude)))
        E_track_per_day = np.zeros((365 + ly, len(latitude), len(longitude)))
        P_per_day = np.zeros((365 + ly, len(latitude), len(longitude)))
        Sa_track_down_per_day = np.zeros((365 + ly, len(latitude), len(longitude)))
        Sa_track_top_per_day = np.zeros((365 + ly, len(latitude), len(longitude)))
        W_down_per_day = np.zeros((365 + ly, len(latitude), len(longitude)))
        W_top_per_day = np.zeros((365 + ly, len(latitude), len(longitude)))
        north_loss_per_day = np.zeros((365 + ly, 1, len(longitude)))
        south_loss_per_day = np.zeros((365 + ly, 1, len(longitude)))
        down_to_top_per_day = np.zeros((365 + ly, len(latitude), len(longitude)))
        top_to_down_per_day = np.zeros((365 + ly, len(latitude), len(longitude)))
        water_lost_per_day = np.zeros((365 + ly, len(latitude), len(longitude)))
        if timetracking == 1:
            Sa_time_down_per_day = np.zeros((365 + ly, len(latitude), len(longitude)))
            Sa_time_top_per_day = np.zeros((365 + ly, len(latitude), len(longitude)))
            E_time_per_day = np.zeros((365 + ly, len(latitude), len(longitude)))

        for a in this_yearpart:
            start = timer()
            datapath = data_path(y, a, years, timetracking)
            if a <= final_time:
                # load tracked data
                loading_ST = sio.loadmat(
                    datapath[0], verify_compressed_data_integrity=False
                )
                Sa_track_top = loading_ST["Sa_track_top"]
                Sa_track_down = loading_ST["Sa_track_down"]
                north_loss = loading_ST["north_loss"]
                south_loss = loading_ST["south_loss"]
                down_to_top = loading_ST["down_to_top"]
                top_to_down = loading_ST["top_to_down"]
                water_lost = loading_ST["water_lost"]
                Sa_track = Sa_track_top + Sa_track_down
                if timetracking == 1:
                    loading_STT = sio.loadmat(
                        datapath[1], verify_compressed_data_integrity=False
                    )
                    Sa_time_top = loading_STT["Sa_time_top"]
                    Sa_time_down = loading_STT["Sa_time_down"]

                # load the total moisture data
                loading_FS = sio.loadmat(
                    datapath[2], verify_compressed_data_integrity=False
                )
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
                E_per_day[a, :, :] = np.sum(E, axis=0)
                E_track_per_day[a, :, :] = np.sum(E_track, axis=0)
                P_per_day[a, :, :] = np.sum(P, axis=0)
                Sa_track_down_per_day[a, :, :] = np.mean(
                    Sa_track_down[1:, :, :], axis=0
                )
                Sa_track_top_per_day[a, :, :] = np.mean(Sa_track_top[1:, :, :], axis=0)
                W_down_per_day[a, :, :] = np.mean(W_down[1:, :, :], axis=0)
                W_top_per_day[a, :, :] = np.mean(W_top[1:, :, :], axis=0)

                north_loss_per_day[a, :, :] = np.sum(north_loss, axis=0)
                south_loss_per_day[a, :, :] = np.sum(south_loss, axis=0)
                down_to_top_per_day[a, :, :] = np.sum(down_to_top, axis=0)
                top_to_down_per_day[a, :, :] = np.sum(top_to_down, axis=0)
                water_lost_per_day[a, :, :] = np.sum(water_lost, axis=0)

                if timetracking == 1:
                    # compute tracked evaporation time
                    E_time = 0.5 * (
                        Sa_time_down[:-1, :, :] + Sa_time_down[1:, :, :]
                    )  # seconds

                    # save per day
                    Sa_time_down_per_day[a, :, :] = np.mean(
                        Sa_time_down[:-1, :, :], axis=0
                    )  # seconds
                    Sa_time_top_per_day[a, :, :] = np.mean(
                        Sa_time_top[:-1, :, :], axis=0
                    )  # seconds
                    E_time_per_day[a, :, :] = (
                        np.sum((E_time * E_track), axis=0) / E_track_per_day[a, :, :]
                    )  # seconds

                    # remove nans
                    where_are_NaNs = np.isnan(E_time_per_day)
                    E_time_per_day[where_are_NaNs] = 0

            end = timer()
            print(
                "Runtime output for day " + str(a + 1) + " in year " + str(y) + " is",
                (end - start),
                " seconds.",
            )

        if daily == 1:
            if timetracking == 0:  # create dummy values
                Sa_time_down_per_day = 0
                Sa_time_top_per_day = 0
                E_time_per_day = 0

            sio.savemat(
                datapath[4],
                {
                    "E_per_day": E_per_day,
                    "E_track_per_day": E_track_per_day,
                    "P_per_day": P_per_day,
                    "Sa_track_down_per_day": Sa_track_down_per_day,
                    "Sa_track_top_per_day": Sa_track_top_per_day,
                    "Sa_time_down_per_day": Sa_time_down_per_day,
                    "Sa_time_top_per_day": Sa_time_top_per_day,
                    "W_down_per_day": W_down_per_day,
                    "W_top_per_day": W_top_per_day,
                    "E_time_per_day": E_time_per_day,
                    "water_lost_per_day": water_lost_per_day,
                },
                do_compression=True,
            )

        # values per month
        for m in range(12):
            first_day = int(datetime.date(y, m + 1, 1).strftime("%j"))
            last_day = int(
                datetime.date(y, m + 1, calendar.monthrange(y, m + 1)[1]).strftime("%j")
            )
            days = (
                np.arange(first_day, last_day + 1) - 1
            )  # -1 because Python is zero-based

            E_per_year_per_month[y - startyear, m, :, :] = np.squeeze(
                np.sum(E_per_day[days, :, :], axis=0)
            )
            E_track_per_year_per_month[y - startyear, m, :, :] = np.squeeze(
                np.sum(E_track_per_day[days, :, :], axis=0)
            )
            P_per_year_per_month[y - startyear, m, :, :] = np.squeeze(
                np.sum(P_per_day[days, :, :], axis=0)
            )
            Sa_track_down_per_year_per_month[y - startyear, m, :, :] = np.squeeze(
                np.mean(Sa_track_down_per_day[days, :, :], axis=0)
            )
            Sa_track_top_per_year_per_month[y - startyear, m, :, :] = np.squeeze(
                np.mean(Sa_track_top_per_day[days, :, :], axis=0)
            )
            W_down_per_year_per_month[y - startyear, m, :, :] = np.squeeze(
                np.mean(W_down_per_day[days, :, :], axis=0)
            )
            W_top_per_year_per_month[y - startyear, m, :, :] = np.squeeze(
                np.mean(W_top_per_day[days, :, :], axis=0)
            )
            north_loss_per_year_per_month[y - startyear, m, :, :] = np.squeeze(
                np.sum(north_loss_per_day[days, :, :], axis=0)
            )
            south_loss_per_year_per_month[y - startyear, m, :, :] = np.squeeze(
                np.sum(south_loss_per_day[days, :, :], axis=0)
            )
            down_to_top_per_year_per_month[y - startyear, m, :, :] = np.squeeze(
                np.sum(down_to_top_per_day[days, :, :], axis=0)
            )
            top_to_down_per_year_per_month[y - startyear, m, :, :] = np.squeeze(
                np.sum(top_to_down_per_day[days, :, :], axis=0)
            )
            water_lost_per_year_per_month[y - startyear, m, :, :] = np.squeeze(
                np.sum(water_lost_per_day[days, :, :], axis=0)
            )

            if timetracking == 1:
                Sa_time_down_per_year_per_month[y - startyear, m, :, :] = np.squeeze(
                    np.mean(Sa_time_down_per_day[days, :, :], axis=0)
                )
                Sa_time_top_per_year_per_month[y - startyear, m, :, :] = np.squeeze(
                    np.mean(Sa_time_top_per_day[days, :, :], axis=0)
                )
                E_time_per_year_per_month[y - startyear, m, :, :] = np.squeeze(
                    np.sum(
                        E_time_per_day[days, :, :] * E_track_per_day[days, :, :], axis=0
                    )
                ) / np.squeeze(E_track_per_year_per_month[y - startyear, m, :, :])

                # remove nans
                where_are_NaNs = np.isnan(E_time_per_year_per_month)
                E_time_per_year_per_month[where_are_NaNs] = 0

            elif timetracking == 0:
                Sa_time_down_per_year_per_month = 0
                Sa_time_top_per_year_per_month = 0
                E_time_per_year_per_month = 0

    # save monthly data
    sio.savemat(
        datapath[3],
        {
            "E_per_year_per_month": E_per_year_per_month,
            "E_track_per_year_per_month": E_track_per_year_per_month,
            "P_per_year_per_month": P_per_year_per_month,
            "Sa_track_down_per_year_per_month": Sa_track_down_per_year_per_month,
            "Sa_track_top_per_year_per_month": Sa_track_top_per_year_per_month,
            "Sa_time_down_per_year_per_month": Sa_time_down_per_year_per_month,
            "Sa_time_top_per_year_per_month": Sa_time_top_per_year_per_month,
            "E_time_per_year_per_month": E_time_per_year_per_month,
            "W_down_per_year_per_month": W_down_per_year_per_month,
            "W_top_per_year_per_month": W_top_per_year_per_month,
            "north_loss_per_year_per_month": north_loss_per_year_per_month,
            "south_loss_per_year_per_month": south_loss_per_year_per_month,
            "down_to_top_per_year_per_month": down_to_top_per_year_per_month,
            "top_to_down_per_year_per_month": top_to_down_per_year_per_month,
            "water_lost_per_year_per_month": water_lost_per_year_per_month,
        },
    )

    end1 = timer()
    print("The total runtime of Con_E_Recyc_Output is", (end1 - start1), " seconds.")
