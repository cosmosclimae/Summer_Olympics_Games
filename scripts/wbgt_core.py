"""
wbgt_core.py
------------
Shared WBGT computation for the Summer Olympics thermal-constraints paper.

Physics (identical to the FIFA / PLOS / heat-stress papers):
  - RH from t2m and d2m
  - Tw : Stull (2011) psychrometric approximation
  - Tg : Hajizadeh et al. (2017) empirical regression  (uses SSRD in W/m^2)
  - WBGT_outdoor = 0.7*Tw + 0.2*Tg + 0.1*Ta        (Yaglou & Minard 1957 weights)
  - WBGT_indoor  : ISO 7243 shade case -> SSRD set to 0 in Tg

Units in:  t2m, d2m in K (ERA5-Land);  ssrd in W/m^2 (already de-accumulated)
Units out: degrees C

Author: CosmosClimae
"""

import numpy as np
import xarray as xr


# ----------------------------------------------------------------------
# Relative humidity from temperature and dewpoint (both in Kelvin)
# ----------------------------------------------------------------------
def relative_humidity(t2m_K, d2m_K):
    """RH (%) from 2 m temperature and dewpoint, Magnus/Tetens over water."""
    Ta = t2m_K - 273.15
    Td = d2m_K - 273.15
    # saturation vapour pressure (hPa), Tetens
    es = 6.112 * np.exp(17.67 * Ta / (Ta + 243.5))
    e  = 6.112 * np.exp(17.67 * Td / (Td + 243.5))
    rh = 100.0 * (e / es)
    return rh.clip(1.0, 100.0)


# ----------------------------------------------------------------------
# Stull (2011) natural wet-bulb temperature
#   Ta in degrees C, RH in %
# ----------------------------------------------------------------------
def wet_bulb_stull(Ta_C, RH):
    """Natural wet-bulb temperature (deg C), Stull 2011."""
    Tw = (Ta_C * np.arctan(0.151977 * np.sqrt(RH + 8.313659))
          + np.arctan(Ta_C + RH)
          - np.arctan(RH - 1.676331)
          + 0.00391838 * np.power(RH, 1.5) * np.arctan(0.023101 * RH)
          - 4.686035)
    return Tw


# ----------------------------------------------------------------------
# Hajizadeh et al. (2017) black-globe temperature
#   ssrd_Wm2 in W/m^2, Ta in deg C, RH in %
# ----------------------------------------------------------------------
def globe_hajizadeh(ssrd_Wm2, Ta_C, RH):
    """Black-globe temperature (deg C), Hajizadeh 2017 regression."""
    Tg = 0.01498 * ssrd_Wm2 + 1.184 * Ta_C - 0.0789 * RH - 2.739
    return Tg


# ----------------------------------------------------------------------
# WBGT outdoor and indoor
# ----------------------------------------------------------------------
def wbgt_outdoor(Ta_C, Tw, Tg):
    return 0.7 * Tw + 0.2 * Tg + 0.1 * Ta_C


def wbgt_indoor(Ta_C, Tw, Tg_noSun):
    """ISO 7243 shade: 0.7 Tw + 0.3 Tg, with Tg computed at SSRD=0."""
    return 0.7 * Tw + 0.3 * Tg_noSun


# ----------------------------------------------------------------------
# Full WBGT from raw ERA5-Land hourly fields
# Returns (wbgt_out, wbgt_in) in deg C
# ----------------------------------------------------------------------
def compute_wbgt(t2m_K, d2m_K, ssrd_Wm2):
    Ta_C = t2m_K - 273.15
    RH   = relative_humidity(t2m_K, d2m_K)
    Tw   = wet_bulb_stull(Ta_C, RH)

    # outdoor: full solar load
    Tg_sun   = globe_hajizadeh(ssrd_Wm2, Ta_C, RH)
    # indoor / shade: no direct solar -> SSRD = 0
    Tg_shade = globe_hajizadeh(xr.zeros_like(ssrd_Wm2), Ta_C, RH) \
               if isinstance(ssrd_Wm2, xr.DataArray) \
               else globe_hajizadeh(0.0, Ta_C, RH)

    wb_out = wbgt_outdoor(Ta_C, Tw, Tg_sun)
    wb_in  = wbgt_indoor(Ta_C, Tw, Tg_shade)
    return wb_out, wb_in


# ----------------------------------------------------------------------
# SSRD de-accumulation
# ERA5-Land accumulates SSRD (J/m^2) from 00 UTC each day.
# Instantaneous hourly mean flux (W/m^2) = (value - previous value)/3600,
# with the daily reset handled: the 01 UTC value is itself the first-hour
# accumulation, so its flux = value/3600.
# ----------------------------------------------------------------------
def deaccumulate_ssrd(ssrd_accum_Jm2, time_dim="valid_time"):
    """
    Convert accumulated SSRD (J/m^2) to instantaneous hourly mean flux (W/m^2).

    ERA5-Land accumulation convention (verified on the data):
      - the accumulation counter resets AT 01 UTC each day;
      - the value at 01 UTC is the first within-day accumulation, i.e. it is
        itself the first hour's flux (00->01 UTC);
      - the value increases through the day and the value stored AT 00 UTC
        (next day) is the FULL 24 h total of the day that just ended;
      - therefore every hour EXCEPT 01 UTC is a simple hourly difference,
        including 00 UTC whose flux = value(00) - value(23) (a real,
        non-zero late-afternoon flux in western-hemisphere longitudes).

    Example (Los Angeles, 2010-07): value(00 UTC, next day)=3.063e7,
    value(23 UTC)=2.825e7 -> flux = 2.38e6 J/m^2 -> ~661 W/m^2 (17h local).
    """
    da = ssrd_accum_Jm2
    hour = da[time_dim].dt.hour

    # hourly difference, first element aligned
    diff = da.diff(time_dim)
    diff = xr.concat([da.isel({time_dim: 0}), diff], dim=time_dim)

    # 01 UTC is the first within-day accumulation -> use the raw value;
    # every other hour (including 00 UTC) uses the hourly difference.
    flux_Jm2 = xr.where(hour == 1, da, diff)

    # Edge case: the very first timestep of the file. If it is 00 UTC it holds
    # the total of a previous day not present in the file, so its "diff" is the
    # raw day-total. Rather than zeroing it (which would drop a real late-
    # afternoon flux in western-hemisphere longitudes, e.g. 1 Aug 00 UTC = 17h
    # local in Los Angeles), we approximate it by the following hour's flux.
    # The residual error is small and is further diluted by the 0.2 globe
    # weight in the outdoor WBGT, so its effect on threshold-day counts is
    # negligible.
    first_is_midnight = bool(hour.isel({time_dim: 0}) == 0)
    if first_is_midnight and da.sizes[time_dim] > 1:
        second = flux_Jm2.isel({time_dim: 1})
        flux_Jm2 = xr.concat(
            [second.expand_dims({time_dim: [da[time_dim].values[0]]}),
             flux_Jm2.isel({time_dim: slice(1, None)})],
            dim=time_dim,
        )

    # guard against tiny negatives from float noise
    flux_Jm2 = flux_Jm2.clip(min=0.0)

    # J/m^2 over one hour -> mean W/m^2
    return flux_Jm2 / 3600.0
