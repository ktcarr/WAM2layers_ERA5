"""backtrack.py: Track moisture back in time, using pre-computed fluxes
   This script takes two arguments: the path to the flux and tracking parameter files"""

import numpy as np
import scipy.io as sio
import os


def create_empty_array(freq, divt, latitude, longitude, fp):
    """Create empty array at specified filepath  (dummy array for backtracking)."""

    ## dimensions for empty array
    dims = (int(freq * divt) + 1, len(latitude), len(longitude))
    Sa_track_top = np.zeros(dims)
    Sa_track_down = np.zeros(dims)

    sio.savemat(
        fp,
        {
            "Sa_track_top": Sa_track_top,
            "Sa_track_down": Sa_track_down,
        },
        do_compression=True,
    )
    return


def data_path(previous_data_to_load, year, doy_idx, fluxes_fp, tracked_moisture_fp):
    load_Sa_track = os.path.join(
        tracked_moisture_fp, previous_data_to_load + "Sa_track.mat"
    )
    load_fluxes_and_storages = os.path.join(
        fluxes_fp, str(year) + "-" + str(doy_idx) + "fluxes_storages.mat"
    )

    save_path_track = os.path.join(
        tracked_moisture_fp, str(year) + "-" + str(doy_idx) + "Sa_track.mat"
    )

    return (
        load_Sa_track,
        load_fluxes_and_storages,
        save_path_track,
    )


#%% Code (no need to look at this for running)
def get_Sa_track_backward(
    latitude,
    longitude,
    freq,
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
    is_global,
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
    if is_global == 1:
        Fa_E_top_boundary[:, :, -1] = 0.5 * (Fa_E_top[:, :, -1] + Fa_E_top[:, :, 0])
    Fa_E_down_boundary = np.zeros(np.shape(Fa_E_down))
    Fa_E_down_boundary[:, :, :-1] = 0.5 * (Fa_E_down[:, :, :-1] + Fa_E_down[:, :, 1:])
    if is_global == 1:
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
    north_loss = np.zeros((int(freq * divt), 1, len(longitude)))
    south_loss = np.zeros((int(freq * divt), 1, len(longitude)))
    down_to_top = np.zeros(np.shape(P))
    top_to_down = np.zeros(np.shape(P))
    water_lost = np.zeros(np.shape(P))
    water_lost_down = np.zeros(np.shape(P))
    water_lost_top = np.zeros(np.shape(P))

    # Sa calculation backward in time
    for t in np.arange(int(freq * divt), 0, -1):
        # down: define values of total moisture
        Sa_E_down[0, :, :-1] = W_down[
            t, :, 1:
        ]  # Atmospheric storage of the cell to the east [m3]
        # to make dependent on is_global but for now kept to avoid division by zero errors
        Sa_E_down[0, :, -1] = W_down[
            t, :, 0
        ]  # Atmospheric storage of the cell to the east [m3]
        Sa_W_down[0, :, 1:] = W_down[
            t, :, :-1
        ]  # Atmospheric storage of the cell to the west [m3]
        # to make dependent on is_global but for now kept to avoid division by zero errors
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
        # to make dependent on is_global but for now kept to avoid division by zero errors
        Sa_E_top[0, :, -1] = W_top[
            t, :, 0
        ]  # Atmospheric storage of the cell to the east [m3]
        Sa_W_top[0, :, 1:] = W_top[
            t, :, :-1
        ]  # Atmospheric storage of the cell to the west [m3]
        # to make dependent on is_global but for now kept to avoid division by zero errors
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
        if is_global == 1:
            Sa_track_E_down[0, :, -1] = Sa_track_down[
                t, :, 0
            ]  # Atmospheric tracked storage of the cell to the east [m3]
        Sa_track_W_down[0, :, 1:] = Sa_track_down[
            t, :, :-1
        ]  # Atmospheric storage of the cell to the west [m3]
        if is_global == 1:
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
        if is_global == 1:
            Sa_track_E_top[0, :, -1] = Sa_track_top[
                t, :, 0
            ]  # Atmospheric tracked storage of the cell to the east [m3]
        Sa_track_W_top[0, :, 1:] = Sa_track_top[
            t, :, :-1
        ]  # Atmospheric tracked storage of the cell to the west [m3]
        if is_global == 1:
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

    #### EXPLANATION OF WATER LOST ####
    # # down and top: water lost to the system:
    # # Can't have more tracked moisture than total moisture
    # water_lost_down[t-1,:,:] = np.maximum(0, Sa_track - W_down)

    # # Stabilize: if more tracked moisture than total moisture,
    # # set tracked moisture to total moisture
    # # otherwise, take tracked moisture
    # Sa_track_down[t-1,1:-1,:] = np.maximum(0, np.minimum(W_down, Sa_track))


if __name__ == "__main__":

    import xarray as xr
    import argparse
    import utils
    from timeit import default_timer as timer

    start1 = timer()

    #### Read parameters #####
    parser = argparse.ArgumentParser()

    ## Specify region to track moisture from
    parser.add_argument("--input_fp", dest="input_fp")
    parser.add_argument("--region_fp", dest="region_fp")
    parser.add_argument("--list_of_days_fp", dest="list_of_days_fp")

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

    ## Numerical parameters
    parser.add_argument("--divt", dest="divt", default=45, type=int)
    parser.add_argument("--kvf", dest="kvf", default=2.0, type=float)
    parser.add_argument("--freq", dest="freq", default=8, type=int)
    parser.add_argument("--is_global", dest="is_global", default=1, type=int)

    ## Data folders
    parser.add_argument("--fluxes_fp", dest="fluxes_fp", type=str)
    parser.add_argument("--tracked_moisture_fp", dest="tracked_moisture_fp", type=str)

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

    # Check if interdata folder exists:
    assert os.path.isdir(args.fluxes_fp), "fluxes_fp does not exist"

    ########### Load list of days to track moisture for, if path is specified ##################
    if args.list_of_days_fp is None:
        extreme_days = None
    else:
        extreme_days = utils.load_doy_list(list_of_days_fp)

    # Create mask for region (load vertices of path from specified file)
    outline = np.loadtxt(args.region_fp, delimiter=",")
    lsm = utils.load_lsm(lsm_fp=os.path.join(args.input_fp, "lsm.nc"))
    Region = utils.makeMask(outline, latitude, longitude, lsm).values

    #%% Runtime & Results
    start1 = timer()

    # Get days of year
    doy_indices = utils.get_doy_indices(args.doy_start, args.doy_end, args.year)
    doy_indices = doy_indices[::-1]  # put in descending order

    for doy_idx in doy_indices: 

        start = timer()

        ## Get filename for data loading
        if utils.is_last_doy_idx(args.year, doy_idx):
            previous_data_to_load = str(args.year + 1) + "-0"
        else:
            previous_data_to_load = str(args.year) + "-" + str(doy_idx + 1)

        ## Get filepaths
        datapath = data_path(
            previous_data_to_load,
            args.year,
            doy_idx,
            args.fluxes_fp,
            args.tracked_moisture_fp,
        )

        ## Create 'dummy' array for first backtracking step
        if doy_idx == doy_indices[0]:
            create_empty_array( 
                args.freq, args.divt, latitude, longitude, datapath[0]
            )

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
            if (args.year, doy_idx) not in extreme_days:
                P = np.zeros_like(P)

        # call the backward tracking function
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
            args.freq,
            args.divt,
            args.kvf,
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
            args.is_global,
        )
        # save this data
        sio.savemat(
            datapath[2],
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

        end = timer()
        print(f"backtrack runtime for day {doy_idx+1}: {end-start:.1f} seconds")

    end1 = timer()
    print(f"Total runtime for backtracking is {end1 - start1:.1f} seconds.")
