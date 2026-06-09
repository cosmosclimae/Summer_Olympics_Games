"""
pipeline1_global_maps.py
------------------------
Pipeline 1 -- GLOBAL GRID for the Figure 1 maps.

For each (year, month) in 2006-2025, July & August:
  - load t2m, d2m, u10, v10, ssrd hourly over the full grid
  - de-accumulate ssrd -> W/m^2
  - compute hourly WBGT outdoor and indoor
  - flag each hour exceeding 28 C
  - per day: at least one hour > 28 C  (boolean per grid cell)
  - sum constrained days over the season
Then average the per-season constrained-day count over the 20 years.

Output: NetCDF with two fields (outdoor, indoor): mean Jul-Aug days with
>=1 h WBGT > 28 C, on the full 0.1 deg grid.

NOTE: heavy. Uses dask chunking. Wind (u10,v10) is loaded for completeness
/ provenance but the Hajizadeh globe does not use wind, so it is not required
for the WBGT itself.
"""

import os
import glob
import argparse
import numpy as np
import xarray as xr
import wbgt_core as w

# ----------------------------------------------------------------------
# Config (overridable by CLI)
# ----------------------------------------------------------------------
MONTHS   = [7, 8]                      # July, August
THRESHOLDS = [28.0, 32.0]              # high-risk and extreme-risk
TIME_DIM = "valid_time"
# Process the global grid in latitude BANDS, in float32, to bound peak memory.
# Each band loads t2m/d2m/ssrd (744 x LAT_BAND x 3600) plus ~6 intermediates.
# float32 + small band keeps the working set ~1-1.5 GB regardless of grid size.
LAT_BAND = 44
FT = np.float32


def _deaccum_np(ssrd, hours):
    """De-accumulate SSRD (numpy float32, time-major) -> W/m^2."""
    diff = np.empty_like(ssrd)
    diff[0] = ssrd[0]
    diff[1:] = ssrd[1:] - ssrd[:-1]
    h = hours
    flux = np.where(h[:, None, None] == 1, ssrd, diff)
    flux = np.where(h[:, None, None] == 0, diff, flux)
    if h[0] == 0:
        flux[0] = flux[1]
    np.clip(flux, 0.0, None, out=flux)
    return (flux / FT(3600.0)).astype(FT)


def season_constrained_days(year, data_dir, domain):
    """Per-grid-cell counts of Jul+Aug days with >=1 h WBGT > threshold,
    for each (formulation, threshold). Memory bounded via float32 lat bands."""
    acc = {}

    for m in MONTHS:
        fname = os.path.join(data_dir, f"era5land_{domain}_{year}-{m:02d}.nc")
        if not os.path.exists(fname):
            raise FileNotFoundError(fname)

        ds = xr.open_dataset(fname)            # lazy header only
        nlat = ds.sizes["latitude"]
        hours = ds[TIME_DIM].dt.hour.values
        day_idx = (ds[TIME_DIM].dt.floor("D")
                   - ds[TIME_DIM].dt.floor("D")[0]).values.astype("timedelta64[D]").astype(int)
        ndays = int(day_idx.max()) + 1

        for lat0 in range(0, nlat, LAT_BAND):
            lat1 = min(lat0 + LAT_BAND, nlat)
            sub = ds.isel(latitude=slice(lat0, lat1))

            t2m = sub["t2m"].values.astype(FT)
            d2m = sub["d2m"].values.astype(FT)
            ssrd = sub["ssrd"].values.astype(FT)

            flux = _deaccum_np(ssrd, hours)
            del ssrd

            Ta = t2m - FT(273.15)
            es = FT(6.112) * np.exp(FT(17.67) * Ta / (Ta + FT(243.5)))
            Td = d2m - FT(273.15)
            e  = FT(6.112) * np.exp(FT(17.67) * Td / (Td + FT(243.5)))
            RH = np.clip(FT(100.0) * e / es, FT(1.0), FT(100.0)).astype(FT)
            del d2m, e, es, Td

            Tw = (Ta * np.arctan(FT(0.151977) * np.sqrt(RH + FT(8.313659)))
                  + np.arctan(Ta + RH) - np.arctan(RH - FT(1.676331))
                  + FT(0.00391838) * RH**FT(1.5) * np.arctan(FT(0.023101) * RH)
                  - FT(4.686035)).astype(FT)
            Tg_sun = (FT(0.01498) * flux + FT(1.184) * Ta
                      - FT(0.0789) * RH - FT(2.739)).astype(FT)
            Tg_sh  = (FT(1.184) * Ta - FT(0.0789) * RH - FT(2.739)).astype(FT)
            del flux

            wb_out = (FT(0.7) * Tw + FT(0.2) * Tg_sun + FT(0.1) * Ta).astype(FT)
            wb_in  = (FT(0.7) * Tw + FT(0.3) * Tg_sh).astype(FT)
            del Ta, RH, Tw, Tg_sun, Tg_sh, t2m

            for thr in THRESHOLDS:
                for label, wb in (("outdoor", wb_out), ("indoor", wb_in)):
                    exc = wb > FT(thr)
                    band_days = np.zeros((lat1 - lat0, exc.shape[2]),
                                         dtype=np.int16)
                    for d in range(ndays):
                        band_days += exc[day_idx == d].any(axis=0)
                    key = f"days_gt{int(thr)}_{label}"
                    if key not in acc:
                        acc[key] = np.zeros((nlat, exc.shape[2]),
                                            dtype=np.float64)
                    acc[key][lat0:lat1] += band_days
                    del exc, band_days
            del wb_out, wb_in

        ds.close()

    # wrap back into DataArrays with coords from a reference file
    ref = xr.open_dataset(
        os.path.join(data_dir, f"era5land_{domain}_{year}-{MONTHS[0]:02d}.nc"))
    coords = {"latitude": ref["latitude"].values,
              "longitude": ref["longitude"].values}
    ref.close()
    out = {k: xr.DataArray(v, dims=("latitude", "longitude"), coords=coords)
           for k, v in acc.items()}
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", required=True)
    ap.add_argument("--domain", default="world")
    ap.add_argument("--years", nargs=2, type=int, default=[2006, 2025],
                    metavar=("START", "END"))
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    years = range(args.years[0], args.years[1] + 1)

    accum = {}
    n = 0
    for yr in years:
        print(f"  processing {yr} ...", flush=True)
        fields = season_constrained_days(yr, args.data_dir, args.domain)
        for k, v in fields.items():
            v = v.compute()
            accum[k] = v if k not in accum else accum[k] + v
        n += 1

    data_vars = {k: (v / n) for k, v in accum.items()}

    result = xr.Dataset(
        data_vars,
        attrs={
            "description": "Mean Jul-Aug days with >=1 h WBGT > threshold",
            "wbgt_outdoor": "0.7 Tw(Stull2011) + 0.2 Tg(Hajizadeh2017) + 0.1 Ta",
            "wbgt_indoor":  "0.7 Tw + 0.3 Tg(SSRD=0), ISO 7243 shade",
            "thresholds_C": ",".join(str(int(t)) for t in THRESHOLDS),
            "years": f"{args.years[0]}-{args.years[1]}",
        },
    )
    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    result.to_netcdf(args.output)
    print(f"Saved {args.output} with vars: {list(data_vars)}")


if __name__ == "__main__":
    main()
