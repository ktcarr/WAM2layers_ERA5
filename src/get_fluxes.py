import numpy as np
import pandas as pd
import xarray as xr
from scipy.interpolate import interp1d
from multiprocessing import Pool
import numpy.ma as ma
from os.path import join

eps = np.finfo(float).eps  # small number to avoid divide-by-zero

####### Functions #######
def remove_negative_longitudes(lon):
    """Convert negative longitude values to positive"""
    lon[lon < 0] = lon[lon < 0] + 360
    return lon


def fix_lonlat(data):
    """Fix longitude/latitude coordinates on xr.dataarray/xr.dataset"""

    # make sure name is "longitude" and "latitude", not "lon"/"lat"
    if "lon" in data.coords:
        data = data.rename({"lon": "longitude", "lat": "latitude"})

    # make sure longitude values are in [0,360) and not [-180,180]
    data["longitude"] = remove_negative_longitudes(data["longitude"].values)

    return data

def safe_open(fp, varname, longitude, latitude, times):
    """Load data from specified fp. Get the specified varname,
    longitude, latitude, and times. 
    Make sure lon/lat dimensions are named consistently,
    and that longitude is in range [0,360)"""
    data = fix_lonlat(xr.open_dataset(fp))
    data = data[varname].sel(
        longitude=longitude,
        latitude=latitude,
        time=times)
    
    return data


def get_fps(name, year, input_fp):
    """Get filepaths for given variable ('name') and year"""
    sep = "-"  # if name=='uvq' else '_'
    fp0 = join(input_fp, f"{year}{sep}{name}.nc")
    fp1 = join(input_fp, f"{year+1}{sep}{name}.nc")
    return fp0, fp1


# other scripts use exactly this sequence, do not change it unless you change it also in the scripts
def get_datapath(date, input_fp, fluxes_fp):
    """Get filepaths for data loading and saving"""

    # Get year, DOY, and DOY-index
    year = date.year
    doy = date.dayofyear
    doy_idx = doy-1

    sp_data, sp_eoy_data = get_fps("sp", year, input_fp)
    tcw_data, tcw_eoy_data = get_fps("tcw", year, input_fp)

    ewvf_data, ewvf_eoy_data = get_fps("vert-int", year, input_fp)
    nwvf_data, nwvf_eoy_data = get_fps("vert-int", year, input_fp)
    eclwf_data, eclwf_eoy_data = get_fps("vert-int", year, input_fp)
    nclwf_data, nclwf_eoy_data = get_fps("vert-int", year, input_fp)
    ecfwf_data, ecfwf_eoy_data = get_fps("vert-int", year, input_fp)
    ncfwf_data, ncfwf_eoy_data = get_fps("vert-int", year, input_fp)

    u_f_data, u_f_eoy_data = get_fps("uvq", year, input_fp)
    v_f_data, v_f_eoy_data = get_fps("uvq", year, input_fp)
    q_f_data, q_f_eoy_data = get_fps("uvq", year, input_fp)

    evaporation_data = get_fps("evaporation", year, input_fp)[0]
    precipitation_data = get_fps("total_precipitation", year, input_fp)[0]

    save_path = join(fluxes_fp, str(year) + "-" + str(doy_idx) + "fluxes_storages.mat")

    return (
        sp_data,
        sp_eoy_data,
        q_f_data,
        q_f_eoy_data,
        tcw_data,
        tcw_eoy_data,
        u_f_data,
        u_f_eoy_data,
        v_f_data,
        v_f_eoy_data,
        ewvf_data,
        ewvf_eoy_data,
        nwvf_data,
        nwvf_eoy_data,
        eclwf_data,
        eclwf_eoy_data,
        nclwf_data,
        nclwf_eoy_data,
        ecfwf_data,
        ecfwf_eoy_data,
        ncfwf_data,
        ncfwf_eoy_data,
        evaporation_data,
        precipitation_data,
        save_path,
    )

def load_data(
    varname,
    fp,
    fp_next,
    latitude,
    longitude,
    date,
    freq,
    get_next_day=True,
    ):
    """Load specified data based, including handling end-of-year stuff.
    'get_next_day' is a boolean specifying whether to get first timestep
    of data from the next day"""

    ## Get timesteps to load
    n_periods = (freq + 1) if get_next_day else freq
    times = pd.date_range(start=date, periods=n_periods, freq=f"{24/freq}H")

    # Check if this is the last day in the year
    is_end_of_year = (times[0].year != times[-1].year)

    if is_end_of_year: 
        data = xr.concat([
            safe_open(fp, varname, longitude, latitude, times[:-1]),
            safe_open(fp_next, varname, longitude, latitude, times[-1])
        ], dim="time",
        )
    else: 
        data = safe_open(fp, varname, longitude, latitude, times)
    
    ## make sure dimension ordering is correct
    data = data.transpose("time",...)

    return data 

def getUVQ(latitude, longitude, date, freq):
    """Load u,v,q data to memory"""
    
    def load(varname, fp, fp_next):
        return load_data(
            varname=varname, 
            fp=fp, 
            fp_next=fp_next, 
            latitude=latitude,
            longitude=longitude,
            date=date,
            freq=freq,
        ) 

    q_f = load("q", datapath[2], datapath[3])
    u_f = load("u", datapath[6], datapath[7])
    v_f = load("v", datapath[8], datapath[9])
    uvq = np.stack([u_f.values, v_f.values, q_f.values], axis=0)

    return uvq


def getPres(
    latitude,
    longitude,
    date,
    freq,
    top_level=100,
    n_levels=37,
):
    """Get pressure levels, after adjusting for surface pressure. Also get pressure difference between levels"""
    # Get surface pressure
    sp = load_data( 
        varname="sp",
        fp=datapath[0],
        fp_next=datapath[1],
        latitude=latitude,
        longitude=longitude,
        date=date,
        freq=freq,
    ).values

    # Get original pressure levels and convert from hPa to Pa
    levels = xr.open_dataset(datapath[2])["level"].values * 100

    # Get new pressure levels (which we'll interpolate to)
    dp = (sp - top_level) / (n_levels - 1)  # get pressure difference between levels
    p = np.stack(
        [sp - n * dp for n in np.arange(n_levels)[::-1]], axis=1
    )  # compute pressure based on dp

    return p, dp, levels


def interp_timeslice(args):
    ti, name = args
    """Interpolate u, v, and q to match surface pressure"""

    new = np.zeros_like(Pres[ti, ...])

    for ilon in range(len(longitude)):
        for ilat in range(len(latitude)):
            data = UVQ_RAW[name, ti, :, ilat, ilon]
            bad_idx = LEVELS > data
            x_masked = ma.MaskedArray(LEVELS, mask=bad_idx)
            y_masked = ma.MaskedArray(data, mask=bad_idx)
            fit = interp1d(x_masked, y_masked, fill_value="extrapolate", kind="linear")
            new[:, ilat, ilon] = fit(Pres[ti, :, ilat, ilon])
    return new


def interp_uvq(freq):
    """Interpolate the UVQ data to match the surface pressure"""

    # Interpolate in parallel (across variables and timesteps
    time_idx = np.arange(freq + 1)
    var_idx = np.arange(3)
    with Pool(len(time_idx)) as p:
        new = p.map(interp_timeslice, [(t, i) for t in time_idx for i in var_idx])

    u_new = np.stack(new[0::3], axis=0)
    v_new = np.stack(new[1::3], axis=0)
    q_new = np.stack(new[2::3], axis=0)

    return u_new, v_new, q_new


def getW(
    q,
    dp,
    latitude,
    longitude,
    date,
    freq,
    density_water,
    g,
    A_gridcell,
    boundary,
):
    """Load W data"""
    tcw_ERA = load_data(
        varname="tcw",
        fp=datapath[4],
        fp_next=datapath[5],
        latitude=latitude,
        longitude=longitude,
        date=date,
        freq=freq,
    ).values
     
    # Compute column water vapor based on defined pressure levels
    q_ = (
        q[:, 1:] + q[:, :-1]
    ) / 2  # estimate average specific humidity between pressure levels
    cwv = q_ * dp[:, None, ...] / g
    tcwv = np.sum(
        cwv, axis=1
    )  # total column water vapor, cwv is summed over the vertical [kg/m2]

    # Compute difference versus actual
    diff = (tcwv - tcw_ERA) / (tcw_ERA + eps)  # avoid divide-by-zero
    diff = np.mean(np.abs(diff)).item()

    # make cw vector (column water based on ratio of estimated to actual tcwv)
    ratio = tcw_ERA / (tcwv + eps)  # avoid divide-by-zero
    cw = ratio[:, None, ...] * cwv

    # put A_gridcell on a 3D grid
    A_gridcell2D = np.tile(A_gridcell, [1, len(longitude)])
    A_gridcell_1_2D = np.reshape(A_gridcell2D, [1, len(latitude), len(longitude)])
    A_gridcell_plus3D = np.tile(A_gridcell_1_2D, [freq + 1, 1, 1])

    # water volumes - might be better to advect interpolated vapor totals here, instead of the ERA totals - not sure how accurate ERA totals are
    vapor_top = np.squeeze(np.sum(cwv[:, :boundary, ...], axis=1))
    vapor_down = np.squeeze(np.sum(cwv[:, boundary:, ...], axis=1))
    vapor = (
        vapor_top + vapor_down + eps
    )  # add small number to avoid divide-by-zero below
    W_top = tcw_ERA * (vapor_top / vapor) * A_gridcell_plus3D / density_water  # m3
    W_down = tcw_ERA * (vapor_down / vapor) * A_gridcell_plus3D / density_water  # m3

    return cw, W_top, W_down


#%% Code
def getFa(
    latitude,
    longitude,
    boundary,
    cw,
    U,
    V,
    freq,
    date,
):
    """Get vertically integrated variables"""

    def load(varname, fp, fp_next):
        return load_data(
            varname=varname,
            fp=fp,
            fp_next=fp_next,
            latitude=latitude,
            longitude=longitude,
            date=date,
            freq=freq,
        )

    # load different variables
    ewvf = load("p71.162", datapath[10], datapath[11]).values
    nwvf = load("p72.162", datapath[12], datapath[13]).values
    eclwf = load("p88.162", datapath[14], datapath[15]).values
    nclwf = load("p89.162", datapath[16], datapath[17]).values
    ecfwf = load("p90.162", datapath[18], datapath[19]).values
    ncfwf = load("p91.162", datapath[20], datapath[21]).values

    ewf = ewvf + eclwf + ecfwf  # kg*m-1*s-1
    nwf = nwvf + nclwf + ncfwf  # kg*m-1*s-1

    # eastward and northward fluxes
    U_ = (
        U[:, 1:] + U[:, :-1]
    ) / 2  # interpolate velocities to match column water estimates
    V_ = (V[:, 1:] + V[:, :-1]) / 2
    Fa_E_p = U_ * cw
    Fa_N_p = V_ * cw

    # uncorrected down and top fluxes
    Fa_E_down_uncorr = np.squeeze(np.sum(Fa_E_p[:, boundary:, :, :], 1))  # kg*m-1*s-1
    Fa_N_down_uncorr = np.squeeze(np.sum(Fa_N_p[:, boundary:, :, :], 1))  # kg*m-1*s-1
    Fa_E_top_uncorr = np.squeeze(np.sum(Fa_E_p[:, 0:boundary, :, :], 1))  # kg*m-1*s-1
    Fa_N_top_uncorr = np.squeeze(np.sum(Fa_N_p[:, 0:boundary, :, :], 1))  # kg*m-1*s-1

    # correct down and top fluxes
    corr_E = np.zeros([freq + 1, len(latitude), len(longitude)])
    corr_N = np.zeros([freq + 1, len(latitude), len(longitude)])
    for i in range(len(longitude)):
        for j in range(len(latitude)):
            for k in range(freq + 1):
                corr_E[k, j, i] = min(
                    2,
                    max(
                        0,
                        ewf[k, j, i]
                        / (Fa_E_down_uncorr[k, j, i] + Fa_E_top_uncorr[k, j, i] + eps),
                    ),
                )
                corr_N[k, j, i] = min(
                    2,
                    max(
                        0,
                        nwf[k, j, i]
                        / (Fa_N_down_uncorr[k, j, i] + Fa_N_top_uncorr[k, j, i] + eps),
                    ),
                )
    Fa_E_down = corr_E * Fa_E_down_uncorr  # kg*m-1*s-1
    Fa_N_down = corr_N * Fa_N_down_uncorr  # kg*m-1*s-1
    Fa_E_top = corr_E * Fa_E_top_uncorr  # kg*m-1*s-1
    Fa_N_top = corr_N * Fa_N_top_uncorr  # kg*m-1*s-1

    # make the fluxes during the timestep
    Fa_E_down = 0.5 * (Fa_E_down[0:-1, :, :] + Fa_E_down[1:, :, :])
    Fa_N_down = 0.5 * (Fa_N_down[0:-1, :, :] + Fa_N_down[1:, :, :])
    Fa_E_top = 0.5 * (Fa_E_top[0:-1, :, :] + Fa_E_top[1:, :, :])
    Fa_N_top = 0.5 * (Fa_N_top[0:-1, :, :] + Fa_N_top[1:, :, :])

    return Fa_E_top, Fa_N_top, Fa_E_down, Fa_N_down


#%% Code
def getEP(latitude, longitude, date, freq_ep, A_gridcell):
    """Load evaporation and precipitation data into memory.
    Note ERA5 precip is hourly, unlike ERA-Interim, which was accumulated in 12-hour periods.
    See https://confluence.ecmwf.int/pages/viewpage.action?pageId=197702790 for details"""

    def load(varname, fp):
        return load_data(
            varname=varname, 
            fp=fp, 
            fp_next=None, 
            latitude=latitude,
            longitude=longitude,
            date=date,
            freq=freq_ep, 
            get_next_day=False,
        )

    evaporation = load("e", datapath[22]).values
    precipitation = load("tp", datapath[23]).values

    # delete and transfer negative values, change sign convention to all positive
    precipitation = np.reshape(
        np.maximum(
            np.reshape(precipitation, (np.size(precipitation)))
            + np.maximum(np.reshape(evaporation, (np.size(evaporation))), 0.0),
            0.0,
        ),
        (int(freq_ep), len(latitude), len(longitude)),
    )
    evaporation = np.reshape(
        np.abs(np.minimum(np.reshape(evaporation, (np.size(evaporation))), 0.0)),
        (int(freq_ep), len(latitude), len(longitude)),
    )

    # calculate volumes
    A_gridcell2D = np.tile(A_gridcell, [1, len(longitude)])
    A_gridcell_1_2D = np.reshape(A_gridcell2D, [1, len(latitude), len(longitude)])
    A_gridcell_max3D = np.tile(A_gridcell_1_2D, [freq_ep, 1, 1])

    E = evaporation * A_gridcell_max3D
    P = precipitation * A_gridcell_max3D

    return E, P


#%% Code
def getrefined(
    Fa_E_top,
    Fa_N_top,
    Fa_E_down,
    Fa_N_down,
    W_top,
    W_down,
    E,
    P,
    divt,
    freq,
    freq_ep,
    latitude,
    longitude,
):

    # for hourly data (i.e., evap/precip)
    # divt2 = divt / 3.0
    divt_ep = divt * freq / freq_ep
    oddvector2 = np.zeros((1, int(freq_ep * divt_ep)))
    partvector2 = np.zeros((1, int(freq_ep * divt_ep)))
    da = np.arange(1, divt_ep)

    for o in np.arange(0, int(freq_ep * divt_ep), int(divt_ep)):
        for i in range(len(da)):
            oddvector2[0, o + i] = (divt_ep - da[i]) / divt_ep
            partvector2[0, o + i + 1] = da[i] / divt_ep

    E_small = np.nan * np.zeros(
        (int(freq_ep * divt_ep), len(latitude), len(longitude))
    )
    for t in range(1, int(freq_ep * divt_ep) + 1):
        E_small[t - 1] = (1.0 / divt_ep) * E[int(t / divt_ep + oddvector2[0, t - 1] - 1)]
    E = E_small

    P_small = np.nan * np.zeros(
        (int(freq_ep * divt_ep), len(latitude), len(longitude))
    )
    for t in range(1, int(freq_ep * divt_ep) + 1):
        P_small[t - 1] = (1.0 / divt_ep) * P[int(t / divt_ep + oddvector2[0, t - 1] - 1)]
    P = P_small

    # for 3 hourly info
    oddvector = np.zeros((1, int(freq * divt)))
    partvector = np.zeros((1, int(freq * divt)))
    da = np.arange(1, divt)
    divt = float(divt)
    for o in np.arange(0, int(freq * divt), int(divt)):
        for i in range(len(da)):
            oddvector[0, o + i] = (divt - da[i]) / divt
            partvector[0, o + i + 1] = da[i] / divt

    W_top_small = np.nan * np.zeros(
        (int(freq * divt + 1), len(latitude), len(longitude))
    )
    for t in range(1, int(freq * divt) + 1):
        W_top_small[t - 1] = W_top[
            int(t / divt + oddvector[0, t - 1] - 1)
        ] + partvector[0, t - 1] * (
            W_top[int(t / divt + oddvector[0, t - 1])]
            - W_top[int(t / divt + oddvector[0, t - 1] - 1)]
        )
        W_top_small[-1] = W_top[-1]
    W_top = W_top_small

    W_down_small = np.nan * np.zeros(
        (int(freq * divt + 1), len(latitude), len(longitude))
    )
    for t in range(1, int(freq * divt) + 1):
        W_down_small[t - 1] = W_down[
            int(t / divt + oddvector[0, t - 1] - 1)
        ] + partvector[0, t - 1] * (
            W_down[int(t / divt + oddvector[0, t - 1])]
            - W_down[int(t / divt + oddvector[0, t - 1] - 1)]
        )
        W_down_small[-1] = W_down[-1]
    W_down = W_down_small

    Fa_E_down_small = np.nan * np.zeros(
        (int(freq * divt), len(latitude), len(longitude))
    )
    Fa_N_down_small = np.nan * np.zeros(
        (int(freq * divt), len(latitude), len(longitude))
    )
    Fa_E_top_small = np.nan * np.zeros(
        (int(freq * divt), len(latitude), len(longitude))
    )
    Fa_N_top_small = np.nan * np.zeros(
        (int(freq * divt), len(latitude), len(longitude))
    )
    for t in range(1, int(freq * divt) + 1):
        Fa_E_down_small[t - 1] = Fa_E_down[int(t / divt + oddvector[0, t - 1] - 1)]
        Fa_N_down_small[t - 1] = Fa_N_down[int(t / divt + oddvector[0, t - 1] - 1)]
        Fa_E_top_small[t - 1] = Fa_E_top[int(t / divt + oddvector[0, t - 1] - 1)]
        Fa_N_top_small[t - 1] = Fa_N_top[int(t / divt + oddvector[0, t - 1] - 1)]

    Fa_E_down = Fa_E_down_small
    Fa_N_down = Fa_N_down_small
    Fa_E_top = Fa_E_top_small
    Fa_N_top = Fa_N_top_small

    return Fa_E_top, Fa_N_top, Fa_E_down, Fa_N_down, E, P, W_top, W_down


#%% Code
def get_stablefluxes(
    W_top,
    W_down,
    Fa_E_top_1,
    Fa_E_down_1,
    Fa_N_top_1,
    Fa_N_down_1,
    divt,
    L_EW_gridcell,
    density_water,
    L_N_gridcell,
    L_S_gridcell,
    latitude,
    freq,
):

    # redefine according to units
    Fa_E_top_kgpmps = Fa_E_top_1
    Fa_E_down_kgpmps = Fa_E_down_1
    Fa_N_top_kgpmps = Fa_N_top_1
    Fa_N_down_kgpmps = Fa_N_down_1

    # compute number of seconds in timestep
    timestep = 1/freq * 24 * 3600 # units: (days) * (hours/day) * (sec/hour) = sec

    # convert to m3
    Fa_E_top = (
        Fa_E_top_kgpmps * timestep / float(divt) * L_EW_gridcell / density_water
    )  # [kg*m^-1*s^-1*s*m*kg^-1*m^3]=[m3]
    Fa_E_down = (
        Fa_E_down_kgpmps * timestep / float(divt) * L_EW_gridcell / density_water
    )  # [s*m*kg*m^-1*s^-1*kg^-1*m^3]=[m3]

    Fa_N_top_swap = np.zeros(
        (len(latitude), int(freq * float(divt)), len(longitude))
    )
    Fa_N_down_swap = np.zeros(
        (len(latitude), int(freq * float(divt)), len(longitude))
    )
    Fa_N_top_kgpmps_swap = np.swapaxes(Fa_N_top_kgpmps, 0, 1)
    Fa_N_down_kgpmps_swap = np.swapaxes(Fa_N_down_kgpmps, 0, 1)
    for c in range(len(latitude)):
        Fa_N_top_swap[c] = (
            Fa_N_top_kgpmps_swap[c]
            * timestep
            / float(divt)
            * 0.5
            * (L_N_gridcell[c] + L_S_gridcell[c])
            / density_water
        )  # [s*m*kg*m^-1*s^-1*kg^-1*m^3]=[m3]
        Fa_N_down_swap[c] = (
            Fa_N_down_kgpmps_swap[c]
            * timestep
            / float(divt)
            * 0.5
            * (L_N_gridcell[c] + L_S_gridcell[c])
            / density_water
        )  # [s*m*kg*m^-1*s^-1*kg^-1*m^3]=[m3]

    Fa_N_top = np.swapaxes(Fa_N_top_swap, 0, 1)
    Fa_N_down = np.swapaxes(Fa_N_down_swap, 0, 1)

    # find out where the negative fluxes are
    Fa_E_top_posneg = np.ones(np.shape(Fa_E_top))
    Fa_E_top_posneg[Fa_E_top < 0] = -1
    Fa_N_top_posneg = np.ones(np.shape(Fa_E_top))
    Fa_N_top_posneg[Fa_N_top < 0] = -1
    Fa_E_down_posneg = np.ones(np.shape(Fa_E_top))
    Fa_E_down_posneg[Fa_E_down < 0] = -1
    Fa_N_down_posneg = np.ones(np.shape(Fa_E_top))
    Fa_N_down_posneg[Fa_N_down < 0] = -1

    # make everything absolute
    Fa_E_top_abs = np.abs(Fa_E_top)
    Fa_E_down_abs = np.abs(Fa_E_down)
    Fa_N_top_abs = np.abs(Fa_N_top)
    Fa_N_down_abs = np.abs(Fa_N_down)

    # stabilize the outfluxes / influxes
    stab = (
        1.0 / 2.0
    )  # during the reduced timestep the water cannot move further than 1/x * the gridcell,
    # in other words at least x * the reduced timestep is needed to cross a gridcell
    Fa_E_top_stable = np.reshape(
        np.minimum(
            np.reshape(Fa_E_top_abs, (np.size(Fa_E_top_abs))),
            (
                np.reshape(Fa_E_top_abs, (np.size(Fa_E_top_abs)))
                / (
                    eps
                    + np.reshape(Fa_E_top_abs, (np.size(Fa_E_top_abs)))
                    + np.reshape(Fa_N_top_abs, (np.size(Fa_N_top_abs)))
                )
            )
            * stab
            * np.reshape(W_top[:-1, :, :], (np.size(W_top[:-1, :, :]))),
        ),
        (int(freq * float(divt)), len(latitude), len(longitude)),
    )
    Fa_N_top_stable = np.reshape(
        np.minimum(
            np.reshape(Fa_N_top_abs, (np.size(Fa_N_top_abs))),
            (
                np.reshape(Fa_N_top_abs, (np.size(Fa_N_top_abs)))
                / (
                    eps
                    + np.reshape(Fa_E_top_abs, (np.size(Fa_E_top_abs)))
                    + np.reshape(Fa_N_top_abs, (np.size(Fa_N_top_abs)))
                )
            )
            * stab
            * np.reshape(W_top[:-1, :, :], (np.size(W_top[:-1, :, :]))),
        ),
        (int(freq * float(divt)), len(latitude), len(longitude)),
    )
    Fa_E_down_stable = np.reshape(
        np.minimum(
            np.reshape(Fa_E_down_abs, (np.size(Fa_E_down_abs))),
            (
                np.reshape(Fa_E_down_abs, (np.size(Fa_E_down_abs)))
                / (
                    eps
                    + np.reshape(Fa_E_down_abs, (np.size(Fa_E_down_abs)))
                    + np.reshape(Fa_N_down_abs, (np.size(Fa_N_down_abs)))
                )
            )
            * stab
            * np.reshape(W_down[:-1, :, :], (np.size(W_down[:-1, :, :]))),
        ),
        (int(freq * float(divt)), len(latitude), len(longitude)),
    )
    Fa_N_down_stable = np.reshape(
        np.minimum(
            np.reshape(Fa_N_down_abs, (np.size(Fa_N_down_abs))),
            (
                np.reshape(Fa_N_down_abs, (np.size(Fa_N_down_abs)))
                / (
                    eps
                    + np.reshape(Fa_E_down_abs, (np.size(Fa_E_down_abs)))
                    + np.reshape(Fa_N_down_abs, (np.size(Fa_N_down_abs)))
                )
            )
            * stab
            * np.reshape(W_down[:-1, :, :], (np.size(W_down[:-1, :, :]))),
        ),
        (int(freq * float(divt)), len(latitude), len(longitude)),
    )

    # get rid of the nan values
    Fa_E_top_stable[np.isnan(Fa_E_top_stable)] = 0
    Fa_N_top_stable[np.isnan(Fa_N_top_stable)] = 0
    Fa_E_down_stable[np.isnan(Fa_E_down_stable)] = 0
    Fa_N_down_stable[np.isnan(Fa_N_down_stable)] = 0

    # redefine
    Fa_E_top = Fa_E_top_stable * Fa_E_top_posneg
    Fa_N_top = Fa_N_top_stable * Fa_N_top_posneg
    Fa_E_down = Fa_E_down_stable * Fa_E_down_posneg
    Fa_N_down = Fa_N_down_stable * Fa_N_down_posneg

    return Fa_E_top, Fa_E_down, Fa_N_top, Fa_N_down


#%% Code
def getFa_Vert(
    Fa_E_top,
    Fa_E_down,
    Fa_N_top,
    Fa_N_down,
    E,
    P,
    W_top,
    W_down,
    divt,
    freq,
    latitude,
    longitude,
    is_global,
):

    # total moisture in the column
    W = W_top + W_down

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

    # check the water balance
    Sa_after_Fa_down = np.zeros([1, len(latitude), len(longitude)])
    Sa_after_Fa_top = np.zeros([1, len(latitude), len(longitude)])
    Sa_after_all_down = np.zeros([1, len(latitude), len(longitude)])
    Sa_after_all_top = np.zeros([1, len(latitude), len(longitude)])
    residual_down = np.zeros(np.shape(P))  # residual factor [m3]
    residual_top = np.zeros(np.shape(P))  # residual factor [m3]

    for t in range(int(freq * divt)):
        # down: calculate with moisture fluxes:
        Sa_after_Fa_down[0, 1:-1, :] = (
            W_down[t, 1:-1, :]
            - Fa_E_down_WE[t, 1:-1, :]
            + Fa_E_down_EW[t, 1:-1, :]
            + Fa_W_down_WE[t, 1:-1, :]
            - Fa_W_down_EW[t, 1:-1, :]
            - Fa_N_down_SN[t, 1:-1, :]
            + Fa_N_down_NS[t, 1:-1, :]
            + Fa_S_down_SN[t, 1:-1, :]
            - Fa_S_down_NS[t, 1:-1, :]
        )

        # top: calculate with moisture fluxes:
        Sa_after_Fa_top[0, 1:-1, :] = (
            W_top[t, 1:-1, :]
            - Fa_E_top_WE[t, 1:-1, :]
            + Fa_E_top_EW[t, 1:-1, :]
            + Fa_W_top_WE[t, 1:-1, :]
            - Fa_W_top_EW[t, 1:-1, :]
            - Fa_N_top_SN[t, 1:-1, :]
            + Fa_N_top_NS[t, 1:-1, :]
            + Fa_S_top_SN[t, 1:-1, :]
            - Fa_S_top_NS[t, 1:-1, :]
        )

        # down: substract precipitation and add evaporation
        Sa_after_all_down[0, 1:-1, :] = (
            Sa_after_Fa_down[0, 1:-1, :]
            - P[t, 1:-1, :] * (W_down[t, 1:-1, :] / W[t, 1:-1, :])
            + E[t, 1:-1, :]
        )

        # top: substract precipitation
        Sa_after_all_top[0, 1:-1, :] = Sa_after_Fa_top[0, 1:-1, :] - P[t, 1:-1, :] * (
            W_top[t, 1:-1, :] / W[t, 1:-1, :]
        )

        # down: calculate the residual
        residual_down[t, 1:-1, :] = (
            W_down[t + 1, 1:-1, :] - Sa_after_all_down[0, 1:-1, :]
        )

        # top: calculate the residual
        residual_top[t, 1:-1, :] = W_top[t + 1, 1:-1, :] - Sa_after_all_top[0, 1:-1, :]

    # compute the resulting vertical moisture flux
    Fa_Vert_raw = (
        W_down[1:, :, :] / W[1:, :, :] * (residual_down + residual_top) - residual_down
    )  # the vertical velocity so that the new residual_down/W_down =  residual_top/W_top (positive downward)

    # find out where the negative vertical flux is
    Fa_Vert_posneg = np.ones(np.shape(Fa_Vert_raw))
    Fa_Vert_posneg[Fa_Vert_raw < 0] = -1

    # make the vertical flux absolute
    Fa_Vert_abs = np.abs(Fa_Vert_raw)

    # stabilize the outfluxes / influxes
    stab = (
        1.0 / 4.0
    )  # during the reduced timestep the vertical flux can maximally empty/fill 1/x of the top or down storage

    Fa_Vert_stable = np.reshape(
        np.minimum(
            np.reshape(Fa_Vert_abs, (np.size(Fa_Vert_abs))),
            np.minimum(
                stab * np.reshape(W_top[1:, :, :], (np.size(W_top[1:, :, :]))),
                stab * np.reshape(W_down[1:, :, :], (np.size(W_down[1:, :, :]))),
            ),
        ),
        (int(freq * float(divt)), len(latitude), len(longitude)),
    )

    # redefine the vertical flux
    Fa_Vert = Fa_Vert_stable * Fa_Vert_posneg

    return Fa_Vert_raw, Fa_Vert


# #### End of functions


#%% Runtime & Results

if __name__ == "__main__":

    import argparse
    import src.utils
    from timeit import default_timer as timer
    import scipy.io as sio

    start1 = timer()

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

    ## Numerical parameters
    parser.add_argument("--divt", dest="divt", type=int, default=45)
    parser.add_argument("--boundary", dest="boundary", type=int, default=29)
    parser.add_argument("--freq", dest="freq", type=int, default=8)
    parser.add_argument("--freq_ep", dest="freq_ep", type=int, default=24)
    parser.add_argument("--is_global", dest="is_global", type=int, default=1)

    ## Data folders
    parser.add_argument("--input_fp", dest="input_fp", type=str)
    parser.add_argument("--fluxes_fp", dest="fluxes_fp", type=str)

    args = parser.parse_args()

    # Get lon/lat, and gridcell dimensions
    constants = src.utils.get_constants_from_args(args)

    # Parse constants
    g = constants["g"]
    density_water = constants["density_water"]
    longitude = constants["longitude"]
    latitude = constants["latitude"]
    A_gridcell = constants["A_gridcell"]
    L_N_gridcell = constants["L_N_gridcell"]
    L_S_gridcell = constants["L_S_gridcell"]
    L_EW_gridcell = constants["L_EW_gridcell"]

    # Get days of year
    # doy_indices = src.utils.get_doy_indices(args.doy_start, args.doy_end, args.year)
    dates = src.utils.get_dates(args.doy_start, args.doy_end, args.year)

    #### Loop through days
    for date in dates:

        start = timer()

        # define save path
        datapath = get_datapath(
            date=date,
            input_fp=args.input_fp,
            fluxes_fp=args.fluxes_fp,
        )  # global variable
 
        # 1 Interpolate U,V,Q data to match surface pressure
        UVQ_RAW = getUVQ(
            latitude=latitude, 
            longitude=longitude, 
            date=date,
            freq=args.freq,
        )

        Pres, DP, LEVELS = getPres(
            latitude=latitude,
            longitude=longitude,
            date=date, 
            freq=args.freq,
            top_level=100,
            n_levels=37,
        )
        U, V, Q = interp_uvq(args.freq)

        # 2 integrate specific humidity to get the (total) column water (vapor)
        cw, W_top, W_down = getW(
            q=Q,
            dp=DP,
            latitude=latitude,
            longitude=longitude,
            date=date, 
            freq=args.freq,
            density_water=density_water,
            g=g,
            A_gridcell=A_gridcell,
            boundary=args.boundary,
        )

        # 3 calculate horizontal moisture fluxes
        Fa_E_top, Fa_N_top, Fa_E_down, Fa_N_down = getFa(
            latitude=latitude,
            longitude=longitude,
            boundary=args.boundary,
            cw=cw,
            U=U,
            V=V,
            freq=args.freq,
            date=date,
        )

        # 4 evaporation and precipitation
        E, P = getEP(
            latitude=latitude, 
            longitude=longitude, 
            date=date,
            freq_ep=args.freq_ep, 
            A_gridcell=A_gridcell,
        )

        # 5 put data on a smaller time step
        (
            Fa_E_top_1,
            Fa_N_top_1,
            Fa_E_down_1,
            Fa_N_down_1,
            E,
            P,
            W_top,
            W_down,
        ) = getrefined(
            Fa_E_top=Fa_E_top,
            Fa_N_top=Fa_N_top,
            Fa_E_down=Fa_E_down,
            Fa_N_down=Fa_N_down,
            W_top=W_top,
            W_down=W_down,
            E=E,
            P=P,
            divt=args.divt,
            freq=args.freq,
            freq_ep=args.freq_ep,
            latitude=latitude,
            longitude=longitude,
        )

        # 6 stabilize horizontal fluxes and get everything in (m3 per smaller timestep)
        Fa_E_top, Fa_E_down, Fa_N_top, Fa_N_down = get_stablefluxes(
            W_top=W_top,
            W_down=W_down,
            Fa_E_top_1=Fa_E_top_1,
            Fa_E_down_1=Fa_E_down_1,
            Fa_N_top_1=Fa_N_top_1,
            Fa_N_down_1=Fa_N_down_1,
            divt=args.divt,
            L_EW_gridcell=L_EW_gridcell,
            density_water=density_water,
            L_N_gridcell=L_N_gridcell,
            L_S_gridcell=L_S_gridcell,
            latitude=latitude,
            freq=args.freq,
        )

        # 7 determine the vertical moisture flux
        Fa_Vert_raw, Fa_Vert = getFa_Vert(
            Fa_E_top=Fa_E_top,
            Fa_E_down=Fa_E_down,
            Fa_N_top=Fa_N_top,
            Fa_N_down=Fa_N_down,
            E=E,
            P=P,
            W_top=W_top,
            W_down=W_down,
            divt=args.divt,
            freq=args.freq,
            latitude=latitude,
            longitude=longitude,
            is_global=args.is_global,
        )

        sio.savemat(
            datapath[24],
            {
                "Fa_E_top": Fa_E_top,
                "Fa_N_top": Fa_N_top,
                "Fa_E_down": Fa_E_down,
                "Fa_N_down": Fa_N_down,
                "Fa_Vert": Fa_Vert,
                "E": E,
                "P": P,
                "W_top": W_top,
                "W_down": W_down,
            },
            do_compression=True,
        )

        end = timer()
        print(f"get_fluxes runtime for {date.year}-{date.month:02d}-{date.day:02d}: {end - start : .1f} sec.")

    end1 = timer()
    print(f"Total get_fluxes runtime for {args.year}: {end1 - start1 : .1f} seconds.")
