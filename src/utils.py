import numpy as np
import xarray as xr
import calendar
import pickle
from matplotlib.path import Path


def load_doy_list(fp):
    """load list of doy indices for moisture tracking"""

    with open(fp, "rb") as f:
        doy_indices = pickle.load(f)

        # if in array format, convert to list of tuples
        if type(doy_indices) is np.ndarray:
            doy_indices = list(zip(doy_indices[:, 0], doy_indices[:, 1]))

    return doy_indices


def load_lsm(lsm_fp):
    """load LSM from filepath"""
    lsm = xr.open_dataarray(lsm_fp)
    lsm = lsm.isel(time=0, drop=True)
    return lsm.astype("bool")


def is_last_doy_idx(year, doy_idx):
    """check if doy_idx is the last for the given year"""
    return doy_idx == 364 + calendar.isleap(year)


def get_next_doy_idx(year, doy_idx):
    """get next year, day of year index"""

    if is_last_doy_idx(year, doy_idx):
        next_year = year + 1
        next_doy_idx = 1
    else:
        next_year = year
        next_doy_idx = doy_idx + 1

    return next_year, next_doy_idx


def makeMask(vertices, lat, lon, lsm=None):
    """Function returns a mask with with 1s representing the area inside of vertices
    'vertices' is an Nx2 array, representing boundaries of a region"""

    # create 2-D grid from lat/lon coords
    lon_lat_grid = np.meshgrid(lon, lat)

    # next, get pairs of lon/lat
    t = zip(lon_lat_grid[0].flatten(), lon_lat_grid[1].flatten())

    # convert back to array
    t = np.array(list(t))  # convert to array

    # convert vertices into a matplotlib Path
    path = Path(vertices)
    mask = path.contains_points(t).reshape(len(lat), len(lon))  # create mask

    # if LSM is supplied, mask out the ocean
    if lsm is not None:
        mask *= lsm.sel(latitude=lat, longitude=lon)

    return mask


## get indices for day of year
def get_doy_indices(doy_start, doy_end, year):
    """Get indices for days of year in given year,
    based on specifed start and end doy"""

    if doy_start is None:
        is_leap_year = int(calendar.isleap(year))
        doy = np.arange(1, 366 + is_leap_year)
    else:
        doy = np.arange(doy_start, doy_end + 1)
    doy_indices = doy - 1  # indices start with 0

    return doy_indices


def get_longitude(lon_min, lon_max, dlon):
    """Get longitude"""
    return np.arange(lon_min, lon_max + 1, dlon)


def get_latitude(lat_min, lat_max, dlat):
    """Get latitude (in descending order)"""
    latitude = np.arange(lat_min, lat_max + 1, dlat)
    return latitude[::-1]


def get_gridcell_dims(latitude, dlon, dlat):
    """Get longitude/latitude range, and gridcell dimensions,
    based on user inputs"""

    ## constants
    dg = 111089.56  # [m] length of 1 degree latitude
    Erad = 6.371e6  # [m] Earth radius

    # new area size calculation:
    lat_n_bound = np.minimum(90.0, latitude + 0.5 * dlat)
    lat_s_bound = np.maximum(-90.0, latitude - 0.5 * dlat)

    A_gridcell = np.zeros([len(latitude), 1])
    A_gridcell[:, 0] = (
        (np.pi / 180.0)
        * Erad**2
        * dlat
        * abs(np.sin(lat_s_bound * np.pi / 180.0) - np.sin(lat_n_bound * np.pi / 180.0))
    )

    # Get length of N/S gridcell boundaries (in meters)
    L_N_gridcell = dlon * np.cos((latitude + dlat / 2.0) * np.pi / 180.0) * dg
    L_S_gridcell = dlon * np.cos((latitude - dlat / 2.0) * np.pi / 180.0) * dg

    # Get length of east/west bondary
    L_EW_gridcell = dlat * dg

    return A_gridcell, L_N_gridcell, L_S_gridcell, L_EW_gridcell


def get_constants(lon_min, lon_max, dlon, lat_min, lat_max, dlat):
    """Wrapper function to get constants based on specified inputs.
    Returns dictionary with constants"""
    latitude = get_latitude(lat_min, lat_max, dlat)
    longitude = get_longitude(lon_min, lon_max, dlon)
    A_gridcell, L_N_gridcell, L_S_gridcell, L_EW_gridcell = get_gridcell_dims(
        latitude=latitude, dlon=dlon, dlat=dlat
    )

    ## Put results in dictionary
    constants = {
        "g": 9.80665,  # [m/s2] from ERA-interim archive
        "density_water": 1000,  # [kg/m3]
        "longitude": longitude,
        "latitude": latitude,
        "A_gridcell": A_gridcell,
        "L_N_gridcell": L_N_gridcell,
        "L_S_gridcell": L_S_gridcell,
        "L_EW_gridcell": L_EW_gridcell,
    }

    return constants


def get_constants_from_args(args):
    """Wrapper function for get_constants, which takes in args structure.
    Useful for parsing user inputs from argparse"""
    constants = get_constants(
        lon_min=args.lon_min,
        lon_max=args.lon_max,
        dlon=args.dlon,
        lat_min=args.lat_min,
        lat_max=args.lat_max,
        dlat=args.dlat,
    )

    return constants


# def get_lsm(lsm_path):
#     '''Tentative function to get LSM'''
#     # Create land-sea-mask (in this model lakes are considered part of the land)
#     # 0 = sea, 1 = land
#     lsm = xr.open_dataset(lsm_path)['lsm'].isel(time=0, drop=True)
#     lsm = lsm.sel(latitude=latitude, longitude=longitude) # determine lon/lat values to use
#     lsm       = lsm.values
#
#     lsm[0,:] = 0 # the northern boundary is always oceanic = 0
#     lsm[-1,:] = 0 # the southern boundary is always oceanic = 0
