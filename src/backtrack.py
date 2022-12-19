"""backtrack.py: Track moisture back in time, using pre-computed fluxes
   This script takes two arguments: the path to the flux and tracking parameter files"""

import numpy as np
import argparse
import scipy.io as sio
import calendar
from getconstants import getconstants
from timeit import default_timer as timer
import os
from matplotlib.path import Path
import pickle
import importlib.util


def makeMask(path, lat, lon):
    #     Function returns a mask with with 1s representing the area inside of path
    lon_lat_grid = np.meshgrid(lon, lat)
    t = zip(
        lon_lat_grid[0].flatten(), lon_lat_grid[1].flatten()
    )  # get pairs of lat/lon
    t = np.array(list(t))  # convert to array
    mask = path.contains_points(t).reshape(len(lat), len(lon))  # create mask
    return mask


# def data_path_ea(years,yearpart):
#     save_empty_arrays_ld_track = os.path.join(sub_interdata_folder, str(years[0]+1) + '-' + str(0) + 'Sa_track.mat')
#     save_empty_arrays_ld_time = os.path.join(sub_interdata_folder, str(years[0]+1) + '-' + str(0) + 'Sa_time.mat')
#
#     save_empty_arrays_track = os.path.join(sub_interdata_folder, str(years[0]) + '-' + str(yearpart[0]+1) + 'Sa_track.mat')
#     save_empty_arrays_time = os.path.join(sub_interdata_folder, str(years[0]) + '-' + str(yearpart[0]+1) + 'Sa_time.mat')
#
#     return save_empty_arrays_ld_track,save_empty_arrays_ld_time,save_empty_arrays_track,save_empty_arrays_time


def data_path(previous_data_to_load, yearnumber, a):
    load_Sa_track = os.path.join(
        sub_interdata_folder, previous_data_to_load + "Sa_track.mat"
    )
    load_fluxes_and_storages = os.path.join(
        interdata_folder, str(yearnumber) + "-" + str(a) + "fluxes_storages.mat"
    )
    load_Sa_time = os.path.join(
        sub_interdata_folder, previous_data_to_load + "Sa_time.mat"
    )

    save_path_track = os.path.join(
        sub_interdata_folder, str(yearnumber) + "-" + str(a) + "Sa_track.mat"
    )
    save_path_time = os.path.join(
        sub_interdata_folder, str(yearnumber) + "-" + str(a) + "Sa_time.mat"
    )
    return (
        load_Sa_track,
        load_fluxes_and_storages,
        load_Sa_time,
        save_path_track,
        save_path_time,
    )


#%% Code (no need to look at this for running)


def get_Sa_track_backward(
    latitude,
    longitude,
    count_time,
    divt,
    Kvf,
    Region,
    Fa_E_top,
    Fa_N_top,
    Fa_E_down,
    Fa_N_down,
    Fa_Vert,
    E,
    P,
    W_top,
    W_down,
    Sa_track_top_last,
    Sa_track_down_last,
):

    # make P_region matrix
    Region3D = np.tile(
        np.reshape(Region, [1, len(latitude), len(longitude)]), [len(P[:, 0, 0]), 1, 1]
    )
    P_region = Region3D * P  # Get precipitation occuring in the sink region

    # Total moisture in the column
    W = W_top + W_down

    # separate the direction of the vertical flux and make it absolute
    Fa_upward = np.zeros(np.shape(Fa_Vert))
    Fa_upward[Fa_Vert <= 0] = Fa_Vert[Fa_Vert <= 0]
    Fa_downward = np.zeros(np.shape(Fa_Vert))
    Fa_downward[Fa_Vert >= 0] = Fa_Vert[Fa_Vert >= 0]
    Fa_upward = np.abs(Fa_upward)

    # include the vertical dispersion
    if Kvf == 0:
        pass
        # do nothing
    else:
        Fa_upward = (1.0 + Kvf) * Fa_upward
        Fa_upward[Fa_Vert >= 0] = Fa_Vert[Fa_Vert >= 0] * Kvf
        Fa_downward = (1.0 + Kvf) * Fa_downward
        Fa_downward[Fa_Vert <= 0] = np.abs(Fa_Vert[Fa_Vert <= 0]) * Kvf

    # define the horizontal fluxes over the boundaries
    # fluxes over the eastern boundary
    Fa_E_top_boundary = np.zeros(np.shape(Fa_E_top))
    Fa_E_top_boundary[:, :, :-1] = 0.5 * (Fa_E_top[:, :, :-1] + Fa_E_top[:, :, 1:])
    if isglobal == 1:
        Fa_E_top_boundary[:, :, -1] = 0.5 * (Fa_E_top[:, :, -1] + Fa_E_top[:, :, 0])
    Fa_E_down_boundary = np.zeros(np.shape(Fa_E_down))
    Fa_E_down_boundary[:, :, :-1] = 0.5 * (Fa_E_down[:, :, :-1] + Fa_E_down[:, :, 1:])
    if isglobal == 1:
        Fa_E_down_boundary[:, :, -1] = 0.5 * (Fa_E_down[:, :, -1] + Fa_E_down[:, :, 0])

    # find out where the positive and negative fluxes are
    Fa_E_top_pos = np.ones(np.shape(Fa_E_top))
    Fa_E_down_pos = np.ones(np.shape(Fa_E_down))
    Fa_E_top_pos[Fa_E_top_boundary < 0] = 0
    Fa_E_down_pos[Fa_E_down_boundary < 0] = 0
    Fa_E_top_neg = Fa_E_top_pos - 1
    Fa_E_down_neg = Fa_E_down_pos - 1

    # separate directions west-east (all positive numbers)
    Fa_E_top_WE = Fa_E_top_boundary * Fa_E_top_pos
    Fa_E_top_EW = Fa_E_top_boundary * Fa_E_top_neg
    Fa_E_down_WE = Fa_E_down_boundary * Fa_E_down_pos
    Fa_E_down_EW = Fa_E_down_boundary * Fa_E_down_neg

    # fluxes over the western boundary
    Fa_W_top_WE = np.nan * np.zeros(np.shape(P))
    Fa_W_top_WE[:, :, 1:] = Fa_E_top_WE[:, :, :-1]
    Fa_W_top_WE[:, :, 0] = Fa_E_top_WE[:, :, -1]
    Fa_W_top_EW = np.nan * np.zeros(np.shape(P))
    Fa_W_top_EW[:, :, 1:] = Fa_E_top_EW[:, :, :-1]
    Fa_W_top_EW[:, :, 0] = Fa_E_top_EW[:, :, -1]
    Fa_W_down_WE = np.nan * np.zeros(np.shape(P))
    Fa_W_down_WE[:, :, 1:] = Fa_E_down_WE[:, :, :-1]
    Fa_W_down_WE[:, :, 0] = Fa_E_down_WE[:, :, -1]
    Fa_W_down_EW = np.nan * np.zeros(np.shape(P))
    Fa_W_down_EW[:, :, 1:] = Fa_E_down_EW[:, :, :-1]
    Fa_W_down_EW[:, :, 0] = Fa_E_down_EW[:, :, -1]

    # fluxes over the northern boundary
    Fa_N_top_boundary = np.nan * np.zeros(np.shape(Fa_N_top))
    Fa_N_top_boundary[:, 1:, :] = 0.5 * (Fa_N_top[:, :-1, :] + Fa_N_top[:, 1:, :])
    Fa_N_down_boundary = np.nan * np.zeros(np.shape(Fa_N_down))
    Fa_N_down_boundary[:, 1:, :] = 0.5 * (Fa_N_down[:, :-1, :] + Fa_N_down[:, 1:, :])

    # find out where the positive and negative fluxes are
    Fa_N_top_pos = np.ones(np.shape(Fa_N_top))
    Fa_N_down_pos = np.ones(np.shape(Fa_N_down))
    Fa_N_top_pos[Fa_N_top_boundary < 0] = 0
    Fa_N_down_pos[Fa_N_down_boundary < 0] = 0
    Fa_N_top_neg = Fa_N_top_pos - 1
    Fa_N_down_neg = Fa_N_down_pos - 1

    # separate directions south-north (all positive numbers)
    Fa_N_top_SN = Fa_N_top_boundary * Fa_N_top_pos
    Fa_N_top_NS = Fa_N_top_boundary * Fa_N_top_neg
    Fa_N_down_SN = Fa_N_down_boundary * Fa_N_down_pos
    Fa_N_down_NS = Fa_N_down_boundary * Fa_N_down_neg

    # fluxes over the southern boundary
    Fa_S_top_SN = np.nan * np.zeros(np.shape(P))
    Fa_S_top_SN[:, :-1, :] = Fa_N_top_SN[:, 1:, :]
    Fa_S_top_NS = np.nan * np.zeros(np.shape(P))
    Fa_S_top_NS[:, :-1, :] = Fa_N_top_NS[:, 1:, :]
    Fa_S_down_SN = np.nan * np.zeros(np.shape(P))
    Fa_S_down_SN[:, :-1, :] = Fa_N_down_SN[:, 1:, :]
    Fa_S_down_NS = np.nan * np.zeros(np.shape(P))
    Fa_S_down_NS[:, :-1, :] = Fa_N_down_NS[:, 1:, :]

    # defining size of output
    Sa_track_down = np.zeros(np.shape(W_down))
    Sa_track_top = np.zeros(np.shape(W_top))

    # assign begin values of output == last (but first index) values of the previous time slot
    Sa_track_down[-1, :, :] = Sa_track_down_last
    Sa_track_top[-1, :, :] = Sa_track_top_last

    # defining sizes of tracked moisture
    Sa_track_after_Fa_down = np.zeros(np.shape(Sa_track_down_last))
    Sa_track_after_Fa_P_E_down = np.zeros(np.shape(Sa_track_down_last))
    Sa_track_E_down = np.zeros(np.shape(Sa_track_down_last))
    Sa_track_W_down = np.zeros(np.shape(Sa_track_down_last))
    Sa_track_N_down = np.zeros(np.shape(Sa_track_down_last))
    Sa_track_S_down = np.zeros(np.shape(Sa_track_down_last))
    Sa_track_after_Fa_top = np.zeros(np.shape(Sa_track_top_last))
    Sa_track_after_Fa_P_E_top = np.zeros(np.shape(Sa_track_top_last))
    Sa_track_E_top = np.zeros(np.shape(Sa_track_top_last))
    Sa_track_W_top = np.zeros(np.shape(Sa_track_top_last))
    Sa_track_N_top = np.zeros(np.shape(Sa_track_top_last))
    Sa_track_S_top = np.zeros(np.shape(Sa_track_top_last))

    # define sizes of total moisture
    Sa_E_down = np.zeros(np.shape(Sa_track_down_last))
    Sa_W_down = np.zeros(np.shape(Sa_track_down_last))
    Sa_N_down = np.zeros(np.shape(Sa_track_down_last))
    Sa_S_down = np.zeros(np.shape(Sa_track_down_last))
    Sa_E_top = np.zeros(np.shape(Sa_track_top_last))
    Sa_W_top = np.zeros(np.shape(Sa_track_top_last))
    Sa_N_top = np.zeros(np.shape(Sa_track_top_last))
    Sa_S_top = np.zeros(np.shape(Sa_track_top_last))

    # define variables that find out what happens to the water
    north_loss = np.zeros((int(count_time * divt), 1, len(longitude)))
    south_loss = np.zeros((int(count_time * divt), 1, len(longitude)))
    down_to_top = np.zeros(np.shape(P))
    top_to_down = np.zeros(np.shape(P))
    water_lost = np.zeros(np.shape(P))
    water_lost_down = np.zeros(np.shape(P))
    water_lost_top = np.zeros(np.shape(P))

    # Sa calculation backward in time
    for t in np.arange(int(count_time * divt), 0, -1):
        # down: define values of total moisture
        Sa_E_down[0, :, :-1] = W_down[
            t, :, 1:
        ]  # Atmospheric storage of the cell to the east [m3]
        # to make dependent on isglobal but for now kept to avoid division by zero errors
        Sa_E_down[0, :, -1] = W_down[
            t, :, 0
        ]  # Atmospheric storage of the cell to the east [m3]
        Sa_W_down[0, :, 1:] = W_down[
            t, :, :-1
        ]  # Atmospheric storage of the cell to the west [m3]
        # to make dependent on isglobal but for now kept to avoid division by zero errors
        Sa_W_down[0, :, 0] = W_down[
            t, :, -1
        ]  # Atmospheric storage of the cell to the west [m3]
        Sa_N_down[0, 1:, :] = W_down[
            t, 0:-1, :
        ]  # Atmospheric storage of the cell to the north [m3]
        Sa_S_down[0, :-1, :] = W_down[
            t, 1:, :
        ]  # Atmospheric storage of the cell to the south [m3]

        # top: define values of total moisture
        Sa_E_top[0, :, :-1] = W_top[
            t, :, 1:
        ]  # Atmospheric storage of the cell to the east [m3]
        # to make dependent on isglobal but for now kept to avoid division by zero errors
        Sa_E_top[0, :, -1] = W_top[
            t, :, 0
        ]  # Atmospheric storage of the cell to the east [m3]
        Sa_W_top[0, :, 1:] = W_top[
            t, :, :-1
        ]  # Atmospheric storage of the cell to the west [m3]
        # to make dependent on isglobal but for now kept to avoid division by zero errors
        Sa_W_top[0, :, 0] = W_top[
            t, :, -1
        ]  # Atmospheric storage of the cell to the west [m3]
        Sa_N_top[0, 1:, :] = W_top[
            t, :-1, :
        ]  # Atmospheric storage of the cell to the north [m3]
        Sa_S_top[0, :-1, :] = W_top[
            t, 1:, :
        ]  # Atmospheric storage of the cell to the south [m3]

        # down: define values of tracked moisture of neighbouring grid cells
        Sa_track_E_down[0, :, :-1] = Sa_track_down[
            t, :, 1:
        ]  # Atmospheric tracked storage of the cell to the east [m3]
        if isglobal == 1:
            Sa_track_E_down[0, :, -1] = Sa_track_down[
                t, :, 0
            ]  # Atmospheric tracked storage of the cell to the east [m3]
        Sa_track_W_down[0, :, 1:] = Sa_track_down[
            t, :, :-1
        ]  # Atmospheric storage of the cell to the west [m3]
        if isglobal == 1:
            Sa_track_W_down[0, :, 0] = Sa_track_down[
                t, :, -1
            ]  # Atmospheric storage of the cell to the west [m3]
        Sa_track_N_down[0, 1:, :] = Sa_track_down[
            t, :-1, :
        ]  # Atmospheric storage of the cell to the north [m3]
        Sa_track_S_down[0, :-1, :] = Sa_track_down[
            t, 1:, :
        ]  # Atmospheric storage of the cell to the south [m3]

        # down: calculate with moisture fluxes
        Sa_track_after_Fa_down[0, 1:-1, :] = (
            Sa_track_down[t, 1:-1, :]
            + Fa_E_down_WE[t - 1, 1:-1, :]
            * (
                Sa_track_E_down[0, 1:-1, :] / Sa_E_down[0, 1:-1, :]
            )  #         (above) Look at cell to the East: compute ratio of tracked to total moisture, and multiply by the flux going west-to-east across the eastern boundary (we are going backwards, so adding this quantity gives us the tracked moisture in the center cell from the previous time step. Repeat for the other boundaries
            - Fa_E_down_EW[t - 1, 1:-1, :]
            * (
                Sa_track_down[t, 1:-1, :] / W_down[t, 1:-1, :]
            )  #         (above) Similary, look at center cell: compute ratio of tracked moisture to total moisture, and look at flux going east to west across the eastern boundary. In reverse, this flux goes west to east, so multiply this ratio by the ratio of moisture in the center cell to get the moisture in the current cell for the the previous time step.
            - Fa_W_down_WE[t - 1, 1:-1, :]
            * (Sa_track_down[t, 1:-1, :] / W_down[t, 1:-1, :])
            + Fa_W_down_EW[t - 1, 1:-1, :]
            * (Sa_track_W_down[0, 1:-1, :] / Sa_W_down[0, 1:-1, :])
            + Fa_N_down_SN[t - 1, 1:-1, :]
            * (Sa_track_N_down[0, 1:-1, :] / Sa_N_down[0, 1:-1, :])
            - Fa_N_down_NS[t - 1, 1:-1, :]
            * (Sa_track_down[t, 1:-1, :] / W_down[t, 1:-1, :])
            - Fa_S_down_SN[t - 1, 1:-1, :]
            * (Sa_track_down[t, 1:-1, :] / W_down[t, 1:-1, :])
            + Fa_S_down_NS[t - 1, 1:-1, :]
            * (Sa_track_S_down[0, 1:-1, :] / Sa_S_down[0, 1:-1, :])
            - Fa_downward[t - 1, 1:-1, :]
            * (Sa_track_down[t, 1:-1, :] / W_down[t, 1:-1, :])
            + Fa_upward[t - 1, 1:-1, :] * (Sa_track_top[t, 1:-1, :] / W_top[t, 1:-1, :])
        )

        # top: define values of tracked moisture of neighbouring grid cells
        Sa_track_E_top[0, :, :-1] = Sa_track_top[
            t, :, 1:
        ]  # Atmospheric tracked storage of the cell to the east [m3]
        if isglobal == 1:
            Sa_track_E_top[0, :, -1] = Sa_track_top[
                t, :, 0
            ]  # Atmospheric tracked storage of the cell to the east [m3]
        Sa_track_W_top[0, :, 1:] = Sa_track_top[
            t, :, :-1
        ]  # Atmospheric tracked storage of the cell to the west [m3]
        if isglobal == 1:
            Sa_track_W_top[0, :, 0] = Sa_track_top[
                t, :, -1
            ]  # Atmospheric tracked storage of the cell to the west [m3]
        Sa_track_N_top[0, 1:, :] = Sa_track_top[
            t, :-1, :
        ]  # Atmospheric tracked storage of the cell to the north [m3]
        Sa_track_S_top[0, :-1, :] = Sa_track_top[
            t, 1:, :
        ]  # Atmospheric tracked storage of the cell to the south [m3]

        # top: calculate with moisture fluxes
        Sa_track_after_Fa_top[0, 1:-1, :] = (
            Sa_track_top[t, 1:-1, :]
            + Fa_E_top_WE[t - 1, 1:-1, :]
            * (Sa_track_E_top[0, 1:-1, :] / Sa_E_top[0, 1:-1, :])
            - Fa_E_top_EW[t - 1, 1:-1, :]
            * (Sa_track_top[t, 1:-1, :] / W_top[t, 1:-1, :])
            - Fa_W_top_WE[t - 1, 1:-1, :]
            * (Sa_track_top[t, 1:-1, :] / W_top[t, 1:-1, :])
            + Fa_W_top_EW[t - 1, 1:-1, :]
            * (Sa_track_W_top[0, 1:-1, :] / Sa_W_top[0, 1:-1, :])
            + Fa_N_top_SN[t - 1, 1:-1, :]
            * (Sa_track_N_top[0, 1:-1, :] / Sa_N_top[0, 1:-1, :])
            - Fa_N_top_NS[t - 1, 1:-1, :]
            * (Sa_track_top[t, 1:-1, :] / W_top[t, 1:-1, :])
            - Fa_S_top_SN[t - 1, 1:-1, :]
            * (Sa_track_top[t, 1:-1, :] / W_top[t, 1:-1, :])
            + Fa_S_top_NS[t - 1, 1:-1, :]
            * (Sa_track_S_top[0, 1:-1, :] / Sa_S_top[0, 1:-1, :])
            + Fa_downward[t - 1, 1:-1, :]
            * (Sa_track_down[t, 1:-1, :] / W_down[t, 1:-1, :])
            - Fa_upward[t - 1, 1:-1, :] * (Sa_track_top[t, 1:-1, :] / W_top[t, 1:-1, :])
        )

        # losses to the north and south
        north_loss[t - 1, 0, :] = Fa_N_top_NS[t - 1, 1, :] * (
            Sa_track_top[t, 1, :] / W_top[t, 1, :]
        ) + Fa_N_down_NS[t - 1, 1, :] * (Sa_track_down[t, 1, :] / W_down[t, 1, :])
        south_loss[t - 1, 0, :] = Fa_S_top_SN[t - 1, -2, :] * (
            Sa_track_top[t, -2, :] / W_top[t, -2, :]
        ) + Fa_S_down_SN[t - 1, -2, :] * (Sa_track_down[t, -2, :] / W_down[t, -2, :])

        # down: add precipitation and subtract evaporation
        Sa_track_after_Fa_P_E_down[0, 1:-1, :] = (
            Sa_track_after_Fa_down[0, 1:-1, :]
            + P_region[t - 1, 1:-1, :] * (W_down[t, 1:-1, :] / W[t, 1:-1, :])
            - E[t - 1, 1:-1, :] * (Sa_track_down[t, 1:-1, :] / W_down[t, 1:-1, :])
        )

        # top: add precipitation
        Sa_track_after_Fa_P_E_top[0, 1:-1, :] = Sa_track_after_Fa_top[
            0, 1:-1, :
        ] + P_region[t - 1, 1:-1, :] * (W_top[t, 1:-1, :] / W[t, 1:-1, :])

        # down and top: redistribute unaccounted water that is otherwise lost from the sytem
        down_to_top[t - 1, :, :] = np.reshape(
            np.maximum(
                0,
                np.reshape(
                    Sa_track_after_Fa_P_E_down, (np.size(Sa_track_after_Fa_P_E_down))
                )
                - np.reshape(W_down[t - 1, :, :], (np.size(W_down[t - 1, :, :]))),
            ),
            (len(latitude), len(longitude)),
        )
        top_to_down[t - 1, :, :] = np.reshape(
            np.maximum(
                0,
                np.reshape(
                    Sa_track_after_Fa_P_E_top, (np.size(Sa_track_after_Fa_P_E_top))
                )
                - np.reshape(W_top[t - 1, :, :], (np.size(W_top[t - 1, :, :]))),
            ),
            (len(latitude), len(longitude)),
        )
        Sa_track_after_all_down = (
            Sa_track_after_Fa_P_E_down
            - down_to_top[t - 1, :, :]
            + top_to_down[t - 1, :, :]
        )
        Sa_track_after_all_top = (
            Sa_track_after_Fa_P_E_top
            - top_to_down[t - 1, :, :]
            + down_to_top[t - 1, :, :]
        )

        # down and top: water lost to the system:
        water_lost_down[t - 1, :, :] = np.reshape(
            np.maximum(
                0,
                np.reshape(Sa_track_after_all_down, (np.size(Sa_track_after_all_down)))
                - np.reshape(W_down[t - 1, :, :], (np.size(W_down[t - 1, :, :]))),
            ),
            (len(latitude), len(longitude)),
        )
        water_lost_top[t - 1, :, :] = np.reshape(
            np.maximum(
                0,
                np.reshape(Sa_track_after_all_top, (np.size(Sa_track_after_all_top)))
                - np.reshape(W_top[t - 1, :, :], (np.size(W_top[t - 1, :, :]))),
            ),
            (len(latitude), len(longitude)),
        )
        water_lost = water_lost_down + water_lost_top

        # down: determine Sa_region of this next timestep 100% stable
        Sa_track_down[t - 1, 1:-1, :] = np.reshape(
            np.maximum(
                0,
                np.minimum(
                    np.reshape(W_down[t - 1, 1:-1, :], np.size(W_down[t - 1, 1:-1, :])),
                    np.reshape(
                        Sa_track_after_all_down[0, 1:-1, :],
                        np.size(Sa_track_after_all_down[0, 1:-1, :]),
                    ),
                ),
            ),
            (len(latitude[1:-1]), len(longitude)),
        )
        # top: determine Sa_region of this next timestep 100% stable
        Sa_track_top[t - 1, 1:-1, :] = np.reshape(
            np.maximum(
                0,
                np.minimum(
                    np.reshape(W_top[t - 1, 1:-1, :], np.size(W_top[t - 1, 1:-1, :])),
                    np.reshape(
                        Sa_track_after_all_top[0, 1:-1, :],
                        np.size(Sa_track_after_all_top[0, 1:-1, :]),
                    ),
                ),
            ),
            (len(latitude[1:-1]), len(longitude)),
        )

    return (
        Sa_track_top,
        Sa_track_down,
        north_loss,
        south_loss,
        down_to_top,
        top_to_down,
        water_lost,
    )

    ##### EXPLANATION OF WATER LOST ####
    ## Note: water lost could probably be reduced by decreasing dt and/or dx
    ## # down and top: water lost to the system:
    ## # Can't have more tracked moisture than total moisture
    ## water_lost_down[t-1,:,:] = np.maximum(0, Sa_track - W_down)
    ##
    ## # Stabilize: if more tracked moisture than total moisture, set tracked moisture to total moisture
    ## # otherwise, take tracked moisture
    ## Sa_track_down[t-1,1:-1,:] = np.maximum(0, np.minimum(W_down, Sa_track))


#%% Code


def get_Sa_track_backward_TIME(
    latitude,
    longitude,
    count_time,
    divt,
    timestep,
    Kvf,
    Region,
    Fa_E_top,
    Fa_N_top,
    Fa_E_down,
    Fa_N_down,
    Fa_Vert,
    E,
    P,
    W_top,
    W_down,
    Sa_track_top_last,
    Sa_track_down_last,
    Sa_time_top_last,
    Sa_time_down_last,
):

    # make P_region matrix
    Region3D = np.tile(
        np.reshape(Region, [1, len(latitude), len(longitude)]), [len(P[:, 0, 0]), 1, 1]
    )
    P_region = Region3D * P

    # Total moisture in the column
    W = W_top + W_down

    # separate the direction of the vertical flux and make it absolute
    Fa_upward = np.zeros(np.shape(Fa_Vert))
    Fa_upward[Fa_Vert <= 0] = Fa_Vert[Fa_Vert <= 0]
    Fa_downward = np.zeros(np.shape(Fa_Vert))
    Fa_downward[Fa_Vert >= 0] = Fa_Vert[Fa_Vert >= 0]
    Fa_upward = np.abs(Fa_upward)

    # include the vertical dispersion
    if Kvf == 0:
        pass
        # do nothing
    else:
        Fa_upward = (1.0 + Kvf) * Fa_upward
        Fa_upward[Fa_Vert >= 0] = Fa_Vert[Fa_Vert >= 0] * Kvf
        Fa_downward = (1.0 + Kvf) * Fa_downward
        Fa_downward[Fa_Vert <= 0] = np.abs(Fa_Vert[Fa_Vert <= 0]) * Kvf

    # define the horizontal fluxes over the boundaries
    # fluxes over the eastern boundary
    Fa_E_top_boundary = np.zeros(np.shape(Fa_E_top))
    Fa_E_top_boundary[:, :, :-1] = 0.5 * (Fa_E_top[:, :, :-1] + Fa_E_top[:, :, 1:])
    if isglobal == 1:
        Fa_E_top_boundary[:, :, -1] = 0.5 * (Fa_E_top[:, :, -1] + Fa_E_top[:, :, 0])
    Fa_E_down_boundary = np.zeros(np.shape(Fa_E_down))
    Fa_E_down_boundary[:, :, :-1] = 0.5 * (Fa_E_down[:, :, :-1] + Fa_E_down[:, :, 1:])
    if isglobal == 1:
        Fa_E_down_boundary[:, :, -1] = 0.5 * (Fa_E_down[:, :, -1] + Fa_E_down[:, :, 0])

    # find out where the positive and negative fluxes are
    Fa_E_top_pos = np.ones(np.shape(Fa_E_top))
    Fa_E_down_pos = np.ones(np.shape(Fa_E_down))
    Fa_E_top_pos[Fa_E_top_boundary < 0] = 0
    Fa_E_down_pos[Fa_E_down_boundary < 0] = 0
    Fa_E_top_neg = Fa_E_top_pos - 1
    Fa_E_down_neg = Fa_E_down_pos - 1

    # separate directions west-east (all positive numbers)
    Fa_E_top_WE = Fa_E_top_boundary * Fa_E_top_pos
    Fa_E_top_EW = Fa_E_top_boundary * Fa_E_top_neg
    Fa_E_down_WE = Fa_E_down_boundary * Fa_E_down_pos
    Fa_E_down_EW = Fa_E_down_boundary * Fa_E_down_neg

    # fluxes over the western boundary
    Fa_W_top_WE = np.nan * np.zeros(np.shape(P))
    Fa_W_top_WE[:, :, 1:] = Fa_E_top_WE[:, :, :-1]
    Fa_W_top_WE[:, :, 0] = Fa_E_top_WE[:, :, -1]
    Fa_W_top_EW = np.nan * np.zeros(np.shape(P))
    Fa_W_top_EW[:, :, 1:] = Fa_E_top_EW[:, :, :-1]
    Fa_W_top_EW[:, :, 0] = Fa_E_top_EW[:, :, -1]
    Fa_W_down_WE = np.nan * np.zeros(np.shape(P))
    Fa_W_down_WE[:, :, 1:] = Fa_E_down_WE[:, :, :-1]
    Fa_W_down_WE[:, :, 0] = Fa_E_down_WE[:, :, -1]
    Fa_W_down_EW = np.nan * np.zeros(np.shape(P))
    Fa_W_down_EW[:, :, 1:] = Fa_E_down_EW[:, :, :-1]
    Fa_W_down_EW[:, :, 0] = Fa_E_down_EW[:, :, -1]

    # fluxes over the northern boundary
    Fa_N_top_boundary = np.nan * np.zeros(np.shape(Fa_N_top))
    Fa_N_top_boundary[:, 1:, :] = 0.5 * (Fa_N_top[:, :-1, :] + Fa_N_top[:, 1:, :])
    Fa_N_down_boundary = np.nan * np.zeros(np.shape(Fa_N_down))
    Fa_N_down_boundary[:, 1:, :] = 0.5 * (Fa_N_down[:, :-1, :] + Fa_N_down[:, 1:, :])

    # find out where the positive and negative fluxes are
    Fa_N_top_pos = np.ones(np.shape(Fa_N_top))
    Fa_N_down_pos = np.ones(np.shape(Fa_N_down))
    Fa_N_top_pos[Fa_N_top_boundary < 0] = 0
    Fa_N_down_pos[Fa_N_down_boundary < 0] = 0
    Fa_N_top_neg = Fa_N_top_pos - 1
    Fa_N_down_neg = Fa_N_down_pos - 1

    # separate directions south-north (all positive numbers)
    Fa_N_top_SN = Fa_N_top_boundary * Fa_N_top_pos
    Fa_N_top_NS = Fa_N_top_boundary * Fa_N_top_neg
    Fa_N_down_SN = Fa_N_down_boundary * Fa_N_down_pos
    Fa_N_down_NS = Fa_N_down_boundary * Fa_N_down_neg

    # fluxes over the southern boundary
    Fa_S_top_SN = np.nan * np.zeros(np.shape(P))
    Fa_S_top_SN[:, :-1, :] = Fa_N_top_SN[:, 1:, :]
    Fa_S_top_NS = np.nan * np.zeros(np.shape(P))
    Fa_S_top_NS[:, :-1, :] = Fa_N_top_NS[:, 1:, :]
    Fa_S_down_SN = np.nan * np.zeros(np.shape(P))
    Fa_S_down_SN[:, :-1, :] = Fa_N_down_SN[:, 1:, :]
    Fa_S_down_NS = np.nan * np.zeros(np.shape(P))
    Fa_S_down_NS[:, :-1, :] = Fa_N_down_NS[:, 1:, :]

    # defining size of output
    Sa_track_down = np.zeros(np.shape(W_down))
    Sa_track_top = np.zeros(np.shape(W_top))
    Sa_time_down = np.zeros(np.shape(W_down))
    Sa_time_top = np.zeros(np.shape(W_top))

    # assign begin values of output == last (but first index) values of the previous time slot
    Sa_track_down[-1, :, :] = Sa_track_down_last
    Sa_track_top[-1, :, :] = Sa_track_top_last
    Sa_time_down[-1, :, :] = Sa_time_down_last
    Sa_time_top[-1, :, :] = Sa_time_top_last

    # defining sizes of tracked moisture
    Sa_track_after_Fa_down = np.zeros(np.shape(Sa_track_down_last))
    Sa_track_after_Fa_P_E_down = np.zeros(np.shape(Sa_track_down_last))
    Sa_track_E_down = np.zeros(np.shape(Sa_track_down_last))
    Sa_track_W_down = np.zeros(np.shape(Sa_track_down_last))
    Sa_track_N_down = np.zeros(np.shape(Sa_track_down_last))
    Sa_track_S_down = np.zeros(np.shape(Sa_track_down_last))
    Sa_track_after_Fa_top = np.zeros(np.shape(Sa_track_top_last))
    Sa_track_after_Fa_P_E_top = np.zeros(np.shape(Sa_track_top_last))
    Sa_track_E_top = np.zeros(np.shape(Sa_track_top_last))
    Sa_track_W_top = np.zeros(np.shape(Sa_track_top_last))
    Sa_track_N_top = np.zeros(np.shape(Sa_track_top_last))
    Sa_track_S_top = np.zeros(np.shape(Sa_track_top_last))

    # define sizes of total moisture
    Sa_E_down = np.zeros(np.shape(Sa_track_down_last))
    Sa_W_down = np.zeros(np.shape(Sa_track_down_last))
    Sa_N_down = np.zeros(np.shape(Sa_track_down_last))
    Sa_S_down = np.zeros(np.shape(Sa_track_down_last))
    Sa_E_top = np.zeros(np.shape(Sa_track_top_last))
    Sa_W_top = np.zeros(np.shape(Sa_track_top_last))
    Sa_N_top = np.zeros(np.shape(Sa_track_top_last))
    Sa_S_top = np.zeros(np.shape(Sa_track_top_last))

    # define variables that find out what happens to the water
    north_loss = np.zeros((int(count_time * divt), 1, len(longitude)))
    south_loss = np.zeros((int(count_time * divt), 1, len(longitude)))
    down_to_top = np.zeros(np.shape(P))
    top_to_down = np.zeros(np.shape(P))
    water_lost = np.zeros(np.shape(P))
    water_lost_down = np.zeros(np.shape(P))
    water_lost_top = np.zeros(np.shape(P))

    # Sa calculation backward in time
    for t in np.arange(int(count_time * divt), 0, -1):
        # down: define values of total moisture
        Sa_E_down[0, :, :-1] = W_down[
            t, :, 1:
        ]  # Atmospheric storage of the cell to the east [m3]
        # to make dependent on isglobal but for now kept to avoid division by zero errors
        Sa_E_down[0, :, -1] = W_down[
            t, :, 0
        ]  # Atmospheric storage of the cell to the east [m3]
        Sa_W_down[0, :, 1:] = W_down[
            t, :, :-1
        ]  # Atmospheric storage of the cell to the west [m3]
        # to make dependent on isglobal but for now kept to avoid division by zero errors
        Sa_W_down[0, :, 0] = W_down[
            t, :, -1
        ]  # Atmospheric storage of the cell to the west [m3]
        Sa_N_down[0, 1:, :] = W_down[
            t, 0:-1, :
        ]  # Atmospheric storage of the cell to the north [m3]
        Sa_S_down[0, :-1, :] = W_down[
            t, 1:, :
        ]  # Atmospheric storage of the cell to the south [m3]

        # top: define values of total moisture
        Sa_E_top[0, :, :-1] = W_top[
            t, :, 1:
        ]  # Atmospheric storage of the cell to the east [m3]
        # to make dependent on isglobal but for now kept to avoid division by zero errors
        Sa_E_top[0, :, -1] = W_top[
            t, :, 0
        ]  # Atmospheric storage of the cell to the east [m3]
        Sa_W_top[0, :, 1:] = W_top[
            t, :, :-1
        ]  # Atmospheric storage of the cell to the west [m3]
        # to make dependent on isglobal but for now kept to avoid division by zero errors
        Sa_W_top[0, :, 0] = W_top[
            t, :, -1
        ]  # Atmospheric storage of the cell to the west [m3]
        Sa_N_top[0, 1:, :] = W_top[
            t, :-1, :
        ]  # Atmospheric storage of the cell to the north [m3]
        Sa_S_top[0, :-1, :] = W_top[
            t, 1:, :
        ]  # Atmospheric storage of the cell to the south [m3]

        # down: define values of tracked moisture of neighbouring grid cells
        Sa_track_E_down[0, :, :-1] = Sa_track_down[
            t, :, 1:
        ]  # Atmospheric tracked storage of the cell to the east [m3]
        if isglobal == 1:
            Sa_track_E_down[0, :, -1] = Sa_track_down[
                t, :, 0
            ]  # Atmospheric tracked storage of the cell to the east [m3]
        Sa_track_W_down[0, :, 1:] = Sa_track_down[
            t, :, :-1
        ]  # Atmospheric storage of the cell to the west [m3]
        if isglobal == 1:
            Sa_track_W_down[0, :, 0] = Sa_track_down[
                t, :, -1
            ]  # Atmospheric storage of the cell to the west [m3]
        Sa_track_N_down[0, 1:, :] = Sa_track_down[
            t, :-1, :
        ]  # Atmospheric storage of the cell to the north [m3]
        Sa_track_S_down[0, :-1, :] = Sa_track_down[
            t, 1:, :
        ]  # Atmospheric storage of the cell to the south [m3]

        # down: calculate with moisture fluxes
        Sa_track_after_Fa_down[0, 1:-1, :] = (
            Sa_track_down[t, 1:-1, :]
            + Fa_E_down_WE[t - 1, 1:-1, :]
            * (Sa_track_E_down[0, 1:-1, :] / Sa_E_down[0, 1:-1, :])
            - Fa_E_down_EW[t - 1, 1:-1, :]
            * (Sa_track_down[t, 1:-1, :] / W_down[t, 1:-1, :])
            - Fa_W_down_WE[t - 1, 1:-1, :]
            * (Sa_track_down[t, 1:-1, :] / W_down[t, 1:-1, :])
            + Fa_W_down_EW[t - 1, 1:-1, :]
            * (Sa_track_W_down[0, 1:-1, :] / Sa_W_down[0, 1:-1, :])
            + Fa_N_down_SN[t - 1, 1:-1, :]
            * (Sa_track_N_down[0, 1:-1, :] / Sa_N_down[0, 1:-1, :])
            - Fa_N_down_NS[t - 1, 1:-1, :]
            * (Sa_track_down[t, 1:-1, :] / W_down[t, 1:-1, :])
            - Fa_S_down_SN[t - 1, 1:-1, :]
            * (Sa_track_down[t, 1:-1, :] / W_down[t, 1:-1, :])
            + Fa_S_down_NS[t - 1, 1:-1, :]
            * (Sa_track_S_down[0, 1:-1, :] / Sa_S_down[0, 1:-1, :])
            - Fa_downward[t - 1, 1:-1, :]
            * (Sa_track_down[t, 1:-1, :] / W_down[t, 1:-1, :])
            + Fa_upward[t - 1, 1:-1, :] * (Sa_track_top[t, 1:-1, :] / W_top[t, 1:-1, :])
        )

        # top: define values of tracked moisture of neighbouring grid cells
        Sa_track_E_top[0, :, :-1] = Sa_track_top[
            t, :, 1:
        ]  # Atmospheric tracked storage of the cell to the east [m3]
        if isglobal == 1:
            Sa_track_E_top[0, :, -1] = Sa_track_top[
                t, :, 0
            ]  # Atmospheric tracked storage of the cell to the east [m3]
        Sa_track_W_top[0, :, 1:] = Sa_track_top[
            t, :, :-1
        ]  # Atmospheric tracked storage of the cell to the west [m3]
        if isglobal == 1:
            Sa_track_W_top[0, :, 0] = Sa_track_top[
                t, :, -1
            ]  # Atmospheric tracked storage of the cell to the west [m3]
        Sa_track_N_top[0, 1:, :] = Sa_track_top[
            t, :-1, :
        ]  # Atmospheric tracked storage of the cell to the north [m3]
        Sa_track_S_top[0, :-1, :] = Sa_track_top[
            t, 1:, :
        ]  # Atmospheric tracked storage of the cell to the south [m3]

        # top: calculate with moisture fluxes
        Sa_track_after_Fa_top[0, 1:-1, :] = (
            Sa_track_top[t, 1:-1, :]
            + Fa_E_top_WE[t - 1, 1:-1, :]
            * (Sa_track_E_top[0, 1:-1, :] / Sa_E_top[0, 1:-1, :])
            - Fa_E_top_EW[t - 1, 1:-1, :]
            * (Sa_track_top[t, 1:-1, :] / W_top[t, 1:-1, :])
            - Fa_W_top_WE[t - 1, 1:-1, :]
            * (Sa_track_top[t, 1:-1, :] / W_top[t, 1:-1, :])
            + Fa_W_top_EW[t - 1, 1:-1, :]
            * (Sa_track_W_top[0, 1:-1, :] / Sa_W_top[0, 1:-1, :])
            + Fa_N_top_SN[t - 1, 1:-1, :]
            * (Sa_track_N_top[0, 1:-1, :] / Sa_N_top[0, 1:-1, :])
            - Fa_N_top_NS[t - 1, 1:-1, :]
            * (Sa_track_top[t, 1:-1, :] / W_top[t, 1:-1, :])
            - Fa_S_top_SN[t - 1, 1:-1, :]
            * (Sa_track_top[t, 1:-1, :] / W_top[t, 1:-1, :])
            + Fa_S_top_NS[t - 1, 1:-1, :]
            * (Sa_track_S_top[0, 1:-1, :] / Sa_S_top[0, 1:-1, :])
            + Fa_downward[t - 1, 1:-1, :]
            * (Sa_track_down[t, 1:-1, :] / W_down[t, 1:-1, :])
            - Fa_upward[t - 1, 1:-1, :] * (Sa_track_top[t, 1:-1, :] / W_top[t, 1:-1, :])
        )

        # losses to the north and south
        north_loss[t - 1, 0, :] = Fa_N_top_NS[t - 1, 1, :] * (
            Sa_track_top[t, 1, :] / W_top[t, 1, :]
        ) + Fa_N_down_NS[t - 1, 1, :] * (Sa_track_down[t, 1, :] / W_down[t, 1, :])
        south_loss[t - 1, 0, :] = Fa_S_top_SN[t - 1, -2, :] * (
            Sa_track_top[t, -2, :] / W_top[t, -2, :]
        ) + Fa_S_down_SN[t - 1, -2, :] * (Sa_track_down[t, -2, :] / W_down[t, -2, :])

        # down: add precipitation and subtract evaporation
        Sa_track_after_Fa_P_E_down[0, 1:-1, :] = (
            Sa_track_after_Fa_down[0, 1:-1, :]
            + P_region[t - 1, 1:-1, :] * (W_down[t, 1:-1, :] / W[t, 1:-1, :])
            - E[t - 1, 1:-1, :] * (Sa_track_down[t, 1:-1, :] / W_down[t, 1:-1, :])
        )

        # top: add precipitation
        Sa_track_after_Fa_P_E_top[0, 1:-1, :] = Sa_track_after_Fa_top[
            0, 1:-1, :
        ] + P_region[t - 1, 1:-1, :] * (W_top[t, 1:-1, :] / W[t, 1:-1, :])

        # down and top: redistribute unaccounted water that is otherwise lost from the sytem
        down_to_top[t - 1, :, :] = np.reshape(
            np.maximum(
                0,
                np.reshape(
                    Sa_track_after_Fa_P_E_down, (np.size(Sa_track_after_Fa_P_E_down))
                )
                - np.reshape(W_down[t - 1, :, :], (np.size(W_down[t - 1, :, :]))),
            ),
            (len(latitude), len(longitude)),
        )
        top_to_down[t - 1, :, :] = np.reshape(
            np.maximum(
                0,
                np.reshape(
                    Sa_track_after_Fa_P_E_top, (np.size(Sa_track_after_Fa_P_E_top))
                )
                - np.reshape(W_top[t - 1, :, :], (np.size(W_top[t - 1, :, :]))),
            ),
            (len(latitude), len(longitude)),
        )
        Sa_track_after_all_down = (
            Sa_track_after_Fa_P_E_down
            - down_to_top[t - 1, :, :]
            + top_to_down[t - 1, :, :]
        )
        Sa_track_after_all_top = (
            Sa_track_after_Fa_P_E_top
            - top_to_down[t - 1, :, :]
            + down_to_top[t - 1, :, :]
        )

        # down and top: water lost to the system:
        water_lost_down[t - 1, :, :] = np.reshape(
            np.maximum(
                0,
                np.reshape(Sa_track_after_all_down, (np.size(Sa_track_after_all_down)))
                - np.reshape(W_down[t - 1, :, :], (np.size(W_down[t - 1, :, :]))),
            ),
            (len(latitude), len(longitude)),
        )
        water_lost_top[t - 1, :, :] = np.reshape(
            np.maximum(
                0,
                np.reshape(Sa_track_after_all_top, (np.size(Sa_track_after_all_top)))
                - np.reshape(W_top[t - 1, :, :], (np.size(W_top[t - 1, :, :]))),
            ),
            (len(latitude), len(longitude)),
        )
        water_lost = water_lost_down + water_lost_top

        # down: determine Sa_region of this next timestep 100% stable
        Sa_track_down[t - 1, 1:-1, :] = np.reshape(
            np.maximum(
                0,
                np.minimum(
                    np.reshape(W_down[t - 1, 1:-1, :], np.size(W_down[t - 1, 1:-1, :])),
                    np.reshape(
                        Sa_track_after_all_down[0, 1:-1, :],
                        np.size(Sa_track_after_all_down[0, 1:-1, :]),
                    ),
                ),
            ),
            (len(latitude[1:-1]), len(longitude)),
        )
        # top: determine Sa_region of this next timestep 100% stable
        Sa_track_top[t - 1, 1:-1, :] = np.reshape(
            np.maximum(
                0,
                np.minimum(
                    np.reshape(W_top[t - 1, 1:-1, :], np.size(W_top[t - 1, 1:-1, :])),
                    np.reshape(
                        Sa_track_after_all_top[0, 1:-1, :],
                        np.size(Sa_track_after_all_top[0, 1:-1, :]),
                    ),
                ),
            ),
            (len(latitude[1:-1]), len(longitude)),
        )

        #############################################################
        # timetracking start

        # defining sizes of timed moisture
        Sa_time_after_Fa_down = np.zeros(np.shape(Sa_time_down_last))
        Sa_time_after_Fa_P_E_down = np.zeros(np.shape(Sa_time_down_last))
        Sa_time_E_down = np.zeros(np.shape(Sa_time_down_last))
        Sa_time_W_down = np.zeros(np.shape(Sa_time_down_last))
        Sa_time_N_down = np.zeros(np.shape(Sa_time_down_last))
        Sa_time_S_down = np.zeros(np.shape(Sa_time_down_last))
        Sa_time_after_Fa_top = np.zeros(np.shape(Sa_time_top_last))
        Sa_time_after_Fa_P_E_top = np.zeros(np.shape(Sa_time_top_last))
        Sa_time_E_top = np.zeros(np.shape(Sa_time_top_last))
        Sa_time_W_top = np.zeros(np.shape(Sa_time_top_last))
        Sa_time_N_top = np.zeros(np.shape(Sa_time_top_last))
        Sa_time_S_top = np.zeros(np.shape(Sa_time_top_last))

        # time increase
        ti = timestep / divt

        # down: define values of timeed moisture of neighbouring grid cells
        Sa_time_E_down[0, :, :-1] = Sa_time_down[
            t, :, 1:
        ]  # Atmospheric timeed storage of the cell to the east [s]
        if isglobal == 1:
            Sa_time_E_down[0, :, -1] = Sa_time_down[
                t, :, 0
            ]  # Atmospheric timeed storage of the cell to the east [s]
        Sa_time_W_down[0, :, 1:] = Sa_time_down[
            t, :, :-1
        ]  # Atmospheric timeed storage of the cell to the west [s]
        if isglobal == 1:
            Sa_time_W_down[0, :, 0] = Sa_time_down[
                t, :, -1
            ]  # Atmospheric timeed storage of the cell to the west [s]
        Sa_time_N_down[0, 1:, :] = Sa_time_down[
            t, :-1, :
        ]  # Atmospheric timeed storage of the cell to the north [s]
        Sa_time_S_down[0, :-1, :] = Sa_time_down[
            t, 1:, :
        ]  # Atmospheric timeed storage of the cell to the south [s]

        # down: calculate with moisture fluxes
        Sa_time_after_Fa_down[0, 1:-1, :] = (
            Sa_track_down[t, 1:-1, :] * (ti + Sa_time_down[t, 1:-1, :])
            + Fa_E_down_WE[t - 1, 1:-1, :]
            * (ti + Sa_time_E_down[0, 1:-1, :])
            * (Sa_track_E_down[0, 1:-1, :] / Sa_E_down[0, 1:-1, :])
            - Fa_E_down_EW[t - 1, 1:-1, :]
            * (ti + Sa_time_down[t, 1:-1, :])
            * (Sa_track_down[t, 1:-1, :] / W_down[t, 1:-1, :])
            - Fa_W_down_WE[t - 1, 1:-1, :]
            * (ti + Sa_time_down[t, 1:-1, :])
            * (Sa_track_down[t, 1:-1, :] / W_down[t, 1:-1, :])
            + Fa_W_down_EW[t - 1, 1:-1, :]
            * (ti + Sa_time_W_down[0, 1:-1, :])
            * (Sa_track_W_down[0, 1:-1, :] / Sa_W_down[0, 1:-1, :])
            + Fa_N_down_SN[t - 1, 1:-1, :]
            * (ti + Sa_time_N_down[0, 1:-1, :])
            * (Sa_track_N_down[0, 1:-1, :] / Sa_N_down[0, 1:-1, :])
            - Fa_N_down_NS[t - 1, 1:-1, :]
            * (ti + Sa_time_down[t, 1:-1, :])
            * (Sa_track_down[t, 1:-1, :] / W_down[t, 1:-1, :])
            - Fa_S_down_SN[t - 1, 1:-1, :]
            * (ti + Sa_time_down[t, 1:-1, :])
            * (Sa_track_down[t, 1:-1, :] / W_down[t, 1:-1, :])
            + Fa_S_down_NS[t - 1, 1:-1, :]
            * (ti + Sa_time_S_down[0, 1:-1, :])
            * (Sa_track_S_down[0, 1:-1, :] / Sa_S_down[0, 1:-1, :])
            - Fa_downward[t - 1, 1:-1, :]
            * (ti + Sa_time_down[t, 1:-1, :])
            * (Sa_track_down[t, 1:-1, :] / W_down[t, 1:-1, :])
            + Fa_upward[t - 1, 1:-1, :]
            * (ti + Sa_time_top[t, 1:-1, :])
            * (Sa_track_top[t, 1:-1, :] / W_top[t, 1:-1, :])
        ) / Sa_track_after_Fa_down[0, 1:-1, :]

        where_are_NaNs = np.isnan(Sa_time_after_Fa_down)
        Sa_time_after_Fa_down[where_are_NaNs] = 0

        # top: define values of timeed moisture of neighbouring grid cells
        Sa_time_E_top[0, :, :-1] = Sa_time_top[
            t, :, 1:
        ]  # Atmospheric storage of the cell to the east [m3]
        if isglobal == 1:
            Sa_time_E_top[0, :, -1] = Sa_time_top[
                t, :, 0
            ]  # Atmospheric storage of the cell to the east [m3]
        Sa_time_W_top[0, :, 1:] = Sa_time_top[
            t, :, :-1
        ]  # Atmospheric storage of the cell to the west [m3]
        if isglobal == 1:
            Sa_time_W_top[0, :, 0] = Sa_time_top[
                t, :, -1
            ]  # Atmospheric storage of the cell to the west [m3]
        Sa_time_N_top[0, 1:, :] = Sa_time_top[
            t, :-1, :
        ]  # Atmospheric storage of the cell to the north [m3]
        Sa_time_S_top[0, :-1, :] = Sa_time_top[
            t, 1:, :
        ]  # Atmospheric storage of the cell to the south [m3]

        # top: calculate with moisture fluxes
        Sa_time_after_Fa_top[0, 1:-1, :] = (
            Sa_track_top[t, 1:-1, :] * (ti + Sa_time_top[t, 1:-1, :])
            + Fa_E_top_WE[t - 1, 1:-1, :]
            * (ti + Sa_time_E_top[0, 1:-1, :])
            * (Sa_track_E_top[0, 1:-1, :] / Sa_E_top[0, 1:-1, :])
            - Fa_E_top_EW[t - 1, 1:-1, :]
            * (ti + Sa_time_top[t, 1:-1, :])
            * (Sa_track_top[t, 1:-1, :] / W_top[t, 1:-1, :])
            - Fa_W_top_WE[t - 1, 1:-1, :]
            * (ti + Sa_time_top[t, 1:-1, :])
            * (Sa_track_top[t, 1:-1, :] / W_top[t, 1:-1, :])
            + Fa_W_top_EW[t - 1, 1:-1, :]
            * (ti + Sa_time_W_top[0, 1:-1, :])
            * (Sa_track_W_top[0, 1:-1, :] / Sa_W_top[0, 1:-1, :])
            + Fa_N_top_SN[t - 1, 1:-1, :]
            * (ti + Sa_time_N_top[0, 1:-1, :])
            * (Sa_track_N_top[0, 1:-1, :] / Sa_N_top[0, 1:-1, :])
            - Fa_N_top_NS[t - 1, 1:-1, :]
            * (ti + Sa_time_top[t, 1:-1, :])
            * (Sa_track_top[t, 1:-1, :] / W_top[t, 1:-1, :])
            - Fa_S_top_SN[t - 1, 1:-1, :]
            * (ti + Sa_time_top[t, 1:-1, :])
            * (Sa_track_top[t, 1:-1, :] / W_top[t, 1:-1, :])
            + Fa_S_top_NS[t - 1, 1:-1, :]
            * (ti + Sa_time_S_top[0, 1:-1, :])
            * (Sa_track_S_top[0, 1:-1, :] / Sa_S_top[0, 1:-1, :])
            + Fa_downward[t - 1, 1:-1, :]
            * (ti + Sa_time_down[t, 1:-1, :])
            * (Sa_track_down[t, 1:-1, :] / W_down[t, 1:-1, :])
            - Fa_upward[t - 1, 1:-1, :]
            * (ti + Sa_time_top[t, 1:-1, :])
            * (Sa_track_top[t, 1:-1, :] / W_top[t, 1:-1, :])
        ) / Sa_track_after_Fa_top[0, 1:-1, :]

        where_are_NaNs = np.isnan(Sa_time_after_Fa_top)
        Sa_time_after_Fa_top[where_are_NaNs] = 0

        # down: add precipitation and substract evaporation
        Sa_time_after_Fa_P_E_down[0, 1:-1, :] = (
            Sa_track_after_Fa_down[0, 1:-1, :] * Sa_time_after_Fa_down[0, 1:-1, :]
            + P_region[t - 1, 1:-1, :] * ti / 2.0 * (W_down[t, 1:-1, :] / W[t, 1:-1, :])
            - E[t - 1, 1:-1, :]
            * (ti + Sa_time_down[t, 1:-1, :])
            * (Sa_track_down[t, 1:-1, :] / W_down[t, 1:-1, :])
        ) / Sa_track_after_Fa_P_E_down[0, 1:-1, :]

        where_are_NaNs = np.isnan(Sa_time_after_Fa_P_E_down)
        Sa_time_after_Fa_P_E_down[where_are_NaNs] = 0

        # top: add precipitation
        Sa_time_after_Fa_P_E_top[0, 1:-1, :] = (
            Sa_track_after_Fa_top[0, 1:-1, :] * Sa_time_after_Fa_top[0, 1:-1, :]
            + P_region[t - 1, 1:-1, :] * ti / 2 * (W_top[t, 1:-1, :] / W[t, 1:-1, :])
        ) / Sa_track_after_Fa_P_E_top[0, 1:-1, :]

        where_are_NaNs = np.isnan(Sa_time_after_Fa_P_E_top)
        Sa_time_after_Fa_P_E_top[where_are_NaNs] = 0

        # down: redistribute water
        Sa_time_after_all_down = (
            Sa_track_after_Fa_P_E_down * Sa_time_after_Fa_P_E_down
            - down_to_top[t - 1, :, :] * Sa_time_after_Fa_P_E_down
            + top_to_down[t - 1, :, :] * Sa_time_after_Fa_P_E_top
        ) / Sa_track_after_all_down

        where_are_NaNs = np.isnan(Sa_time_after_all_down)
        Sa_time_after_all_down[where_are_NaNs] = 0

        # top: redistribute water
        Sa_time_after_all_top = (
            Sa_track_after_Fa_P_E_top * Sa_time_after_Fa_P_E_top
            - top_to_down[t - 1, :, :] * Sa_time_after_Fa_P_E_top
            + down_to_top[t - 1, :, :] * Sa_time_after_Fa_P_E_down
        ) / Sa_track_after_all_top

        where_are_NaNs = np.isnan(Sa_time_after_all_top)
        Sa_time_after_all_top[where_are_NaNs] = 0

        # down: determine Sa_region of this next timestep 100% stable
        Sa_time_down[t - 1, 1:-1, :] = Sa_time_after_all_down[0, 1:-1, :]

        # top: determine Sa_region of this next timestep 100% stable
        Sa_time_top[t - 1, 1:-1, :] = Sa_time_after_all_top[0, 1:-1, :]
        #############################################################

    return (
        Sa_time_top,
        Sa_time_down,
        Sa_track_top,
        Sa_track_down,
        north_loss,
        south_loss,
        down_to_top,
        top_to_down,
        water_lost,
    )


#%% create empty array for track and time


def create_empty_array(count_time, divt, latitude, longitude, year, doy):
    """Create empty array using specified year and doy (dummy array for backtracking)."""
    Sa_time_top = np.zeros((int(count_time * divt) + 1, len(latitude), len(longitude)))
    Sa_time_down = np.zeros((int(count_time * divt) + 1, len(latitude), len(longitude)))
    Sa_track_top = np.zeros((int(count_time * divt) + 1, len(latitude), len(longitude)))
    Sa_track_down = np.zeros(
        (int(count_time * divt) + 1, len(latitude), len(longitude))
    )
    sio.savemat(
        os.path.join(sub_interdata_folder, f"{year}-{doy}Sa_track.mat"),
        {
            "Sa_track_top": Sa_track_top,
            "Sa_track_down": Sa_track_down,
        },
        do_compression=True,
    )
    return


######################################### start of tracking ####################################################

if __name__ == "__main__":

    #### Read parameters #####
    parser = argparse.ArgumentParser()
    parser.add_argument("--params_fp", dest="params_path")
    parser.add_argument("--year1", dest="year1", type=int)
    parser.add_argument("--year2", dest="year2", type=int)
    parser.add_argument("--veryfirstrun", dest="veryfirstrun", type=int)
    parser.add_argument("--region_fp", dest="region_fp")
    parser.add_argument("--list_of_days_fp", dest="list_of_days_fp", default="none")
    parser.add_argument("--Kvf", dest="Kvf", default=3.0, type=float)
    args = parser.parse_args()

    params_path = args.params_path
    veryfirstrun = args.veryfirstrun
    region_fp = args.region_fp
    list_of_days_fp = args.list_of_days_fp
    years = np.arange(args.year1, args.year2 + 1)
    Kvf = args.Kvf

    ###  Import parameters used for flux computations  ##############
    spec = importlib.util.spec_from_file_location("module.name", params_path)
    params = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(params)

    ######### Parameters (shared for all experiments) ###################
    latitude = params.latitude
    longitude = params.longitude
    lsm_path = params.lsm_path
    divt = params.divt
    count_time = params.count_time
    isglobal = params.isglobal
    interdata_folder = params.interdata_folder
    sub_interdata_folder = os.path.join(params.output_folder, f"kvf_{int(Kvf)}", "back")
    lake_mask = np.nan  # dummy argument for compatibility with other stuff...

    timetracking = 0  # 0 for not tracking time and 1 for tracking time

    # Check if interdata folder exists:
    assert os.path.isdir(
        interdata_folder
    ), "Please create the interdata_folder before running the script"
    # Check if sub interdata folder exists otherwise create it:
    if os.path.isdir(sub_interdata_folder):
        pass
    else:
        os.makedirs(sub_interdata_folder)

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

    ########### Load list of days to track moisture for, if path is specified ##################
    if list_of_days_fp == "none":
        extreme_days = None
    else:
        with open(list_of_days_fp, "rb") as f:
            extreme_days = pickle.load(f)

            # if in array format, convert to list of tuples
            if type(extreme_days) is np.ndarray:
                extreme_days = list(zip(extreme_days[:, 0], extreme_days[:, 1]))

    # Create mask for region (load vertices of path from specified file)
    region_path = Path(np.load(region_fp))
    Region = lsm * makeMask(region_path, latitude, longitude)

    #%% Runtime & Results
    start1 = timer()

    # loop through the years (start with most recent first)
    yearnumber = years[0]  # this version of the script only takes 1 year

    # ## Get indices of valid days
    # if yearnumber==2021:
    #     last_day = 90
    # elif calendar.isleap(yearnumber):
    #     last_day = 366
    # else:
    #     last_day = 365
    #
    # thisyearpart = np.arange(last_day)[::-1]
    thisyearpart = np.arange(91)[::-1]
    # doy_start = 90 if yearnumber==2021 else 0
    doy_start = 91
    # year_start= yearnumber if (yearnumber==2021) else yearnumber+1
    year_start = yearnumber

    # # only use for short experiment
    # thisyearpart = np.arange(183,214)[::-1]
    # doy_start = 214
    # year_start = 1993

    create_empty_array(
        count_time - 1, divt, latitude, longitude, year_start, doy_start
    )  # creates empty arrays for first day run

    for a in thisyearpart:
        # last day has one fewer timestep (no data for April 1 2021)
        count_time_ = count_time - 1 if (a == thisyearpart[0]) else count_time

        start = timer()

        if a == (364 + calendar.isleap(yearnumber)):  # a == 31 December
            previous_data_to_load = str(yearnumber + 1) + "-0"
        else:  # a != 31 December
            previous_data_to_load = str(yearnumber) + "-" + str(a + 1)
        datapath = data_path(previous_data_to_load, yearnumber, a)

        loading_ST = sio.loadmat(datapath[0], verify_compressed_data_integrity=False)
        Sa_track_top = loading_ST["Sa_track_top"]
        Sa_track_down = loading_ST["Sa_track_down"]
        Sa_track_top_last_scheef = Sa_track_top[0, :, :]
        Sa_track_down_last_scheef = Sa_track_down[0, :, :]
        Sa_track_top_last = np.reshape(
            Sa_track_top_last_scheef, (1, len(latitude), len(longitude))
        )
        Sa_track_down_last = np.reshape(
            Sa_track_down_last_scheef, (1, len(latitude), len(longitude))
        )

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

        # Check for extremes: if criterion not satisfied, zero out the precip
        if extreme_days is not None:
            if (yearnumber, a) not in extreme_days:
                P = np.zeros_like(P)

        # call the backward tracking function
        if timetracking == 0:
            (
                Sa_track_top,
                Sa_track_down,
                north_loss,
                south_loss,
                down_to_top,
                top_to_down,
                water_lost,
            ) = get_Sa_track_backward(
                latitude,
                longitude,
                count_time_,
                divt,
                Kvf,
                Region,
                Fa_E_top,
                Fa_N_top,
                Fa_E_down,
                Fa_N_down,
                Fa_Vert,
                E,
                P,
                W_top,
                W_down,
                Sa_track_top_last,
                Sa_track_down_last,
            )
        elif timetracking == 1:
            loading_STT = sio.loadmat(
                datapath[2], verify_compressed_data_integrity=False
            )
            Sa_time_top = loading_STT["Sa_time_top"]  # [seconds]
            Sa_time_down = loading_STT["Sa_time_down"]
            Sa_time_top_last_scheef = Sa_time_top[0, :, :]
            Sa_time_down_last_scheef = Sa_time_down[0, :, :]
            Sa_time_top_last = np.reshape(
                Sa_time_top_last_scheef, (1, len(latitude), len(longitude))
            )
            Sa_time_down_last = np.reshape(
                Sa_time_down_last_scheef, (1, len(latitude), len(longitude))
            )

            (
                Sa_time_top,
                Sa_time_down,
                Sa_track_top,
                Sa_track_down,
                north_loss,
                south_loss,
                down_to_top,
                top_to_down,
                water_lost,
            ) = get_Sa_track_backward_TIME(
                latitude,
                longitude,
                count_time_,
                divt,
                timestep,
                Kvf,
                Region,
                Fa_E_top,
                Fa_N_top,
                Fa_E_down,
                Fa_N_down,
                Fa_Vert,
                E,
                P,
                W_top,
                W_down,
                Sa_track_top_last,
                Sa_track_down_last,
                Sa_time_top_last,
                Sa_time_down_last,
            )

        # save this data
        sio.savemat(
            datapath[3],
            {
                "Sa_track_top": Sa_track_top,
                "Sa_track_down": Sa_track_down,
                "north_loss": north_loss,
                "south_loss": south_loss,
                "down_to_top": down_to_top,
                "top_to_down": top_to_down,
                "water_lost": water_lost,
            },
            do_compression=True,
        )
        if timetracking == 1:
            sio.savemat(
                datapath[4],
                {"Sa_time_top": Sa_time_top, "Sa_time_down": Sa_time_down},
                do_compression=True,
            )

        end = timer()
        print(
            "Runtime Sa_track for day "
            + str(a + 1)
            + " in year "
            + str(yearnumber)
            + " is",
            (end - start),
            " seconds.",
        )

    end1 = timer()
    print("The total runtime of precipitation_shed is", (end1 - start1), " seconds.")
