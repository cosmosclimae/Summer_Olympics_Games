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
THRESH   = 28.0                        # deg C
TIME_DIM = "valid_time"
CHUNKS   = {TIME_DIM: 24, "latitude": 363, "longitude": 900}  # tune to RAM


def season_constrained_days(year, data_dir, domain):
    """Return (days_out, days_in): per-grid-cell count of Jul+Aug days
    with at least one hour WBGT > THRESH, for one year."""
    days_out_list, days_in_list = [], []

    for m in MONTHS:
        fname = os.path.join(data_dir, f"era5land_{domain}_{year}-{m:02d}.nc")
        if not os.path.exists(fname):
            raise FileNotFoundError(fname)

        ds = xr.open_dataset(fname, chunks=CHUNKS)

        ssrd_flux = w.deaccumulate_ssrd(ds["ssrd"], time_dim=TIME_DIM)
        wb_out, wb_in = w.compute_wbgt(ds["t2m"], ds["d2m"], ssrd_flux)

        # hourly exceedance
        exc_out = wb_out > THRESH
        exc_in  = wb_in  > THRESH

        # group hours by calendar day, "any hour exceeds" per day
        day = ds[TIME_DIM].dt.floor("D")
        exc_out = exc_out.assign_coords(day=day)
        exc_in  = exc_in.assign_coords(day=day)

        any_out = exc_out.groupby("day").any(dim=TIME_DIM)
        any_in  = exc_in.groupby("day").any(dim=TIME_DIM)

        # count constrained days in this month
        days_out_list.append(any_out.sum(dim="day"))
        days_in_list.append(any_in.sum(dim="day"))

        ds.close()

    days_out = sum(days_out_list)   # Jul + Aug
    days_in  = sum(days_in_list)
    return days_out, days_in


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", required=True)
    ap.add_argument("--domain", default="world")
    ap.add_argument("--years", nargs=2, type=int, default=[2006, 2025],
                    metavar=("START", "END"))
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    years = range(args.years[0], args.years[1] + 1)

    out_accum, in_accum = None, None
    n = 0
    for yr in years:
        print(f"  processing {yr} ...", flush=True)
        d_out, d_in = season_constrained_days(yr, args.data_dir, args.domain)
        d_out = d_out.compute()
        d_in  = d_in.compute()
        out_accum = d_out if out_accum is None else out_accum + d_out
        in_accum  = d_in  if in_accum  is None else in_accum  + d_in
        n += 1

    mean_out = out_accum / n
    mean_in  = in_accum  / n

    result = xr.Dataset(
        {
            "days_gt28_outdoor": mean_out,
            "days_gt28_indoor":  mean_in,
        },
        attrs={
            "description": "Mean Jul-Aug days with >=1 h WBGT > 28 C",
            "wbgt_outdoor": "0.7 Tw(Stull2011) + 0.2 Tg(Hajizadeh2017) + 0.1 Ta",
            "wbgt_indoor":  "0.7 Tw + 0.3 Tg(SSRD=0), ISO 7243 shade",
            "threshold_C": THRESH,
            "years": f"{args.years[0]}-{args.years[1]}",
        },
    )
    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    result.to_netcdf(args.output)
    print(f"Saved {args.output}")


if __name__ == "__main__":
    main()
