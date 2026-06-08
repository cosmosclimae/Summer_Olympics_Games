"""
diag_landmask.py
----------------
Inspect the ERA5-Land land mask around the problem cities and report the
nearest VALID (non-NaN) land grid cell, with its distance from the requested
point. Run this to choose clean coordinates for coastal cities.

Usage:
  python diag_landmask.py --data-dir /mnt/f/climae_CIO/raw/world/era5land --domain world
"""
import os
import glob
import argparse
import numpy as np
import xarray as xr

TIME_DIM = "valid_time"

CITIES = {
    "Rio de Janeiro": (-22.91, -43.20),
    "Tokyo":          ( 35.68, 139.69),
    "Paris":          ( 48.86,   2.35),
    "Los Angeles":    ( 34.05, -118.24),
    "Brisbane":       (-27.47, 153.03),
    "Ahmedabad":      ( 23.03,  72.58),
    "Doha":           ( 25.29,  51.53),
    "Cape Town":      (-33.92,  18.42),
}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", required=True)
    ap.add_argument("--domain", default="world")
    ap.add_argument("--window", type=float, default=0.6,
                    help="half-width of search window in degrees")
    args = ap.parse_args()

    f = sorted(glob.glob(os.path.join(args.data_dir,
                                      f"era5land_{args.domain}_*.nc")))[0]
    print(f"[diag] file: {f}")
    ds = xr.open_dataset(f)
    land0 = ds["t2m"].isel({TIME_DIM: 0})  # first hour, NaN = sea

    print(f"{'city':18s} {'nearest_NaN?':12s} {'land_lat':>9s} {'land_lon':>9s} {'dist_km':>8s}")
    for city, (lat, lon) in CITIES.items():
        near = land0.sel(latitude=lat, longitude=lon, method="nearest")
        is_nan = bool(np.isnan(near))

        w = args.window
        win = land0.sel(latitude=slice(lat + w, lat - w),     # lat descending
                        longitude=slice(lon - w, lon + w))
        valid = win.notnull()
        if bool(valid.any()):
            la = win["latitude"].values
            lo = win["longitude"].values
            LO, LA = np.meshgrid(lo, la)
            mask = valid.values
            dist = np.sqrt((LA - lat) ** 2 + (LO - lon) ** 2)
            dist[~mask] = np.inf
            j, i = np.unravel_index(np.argmin(dist), dist.shape)
            blat, blon = la[j], lo[i]
            dkm = np.sqrt((blat - lat) ** 2 + (blon - lon) ** 2) * 111.0
            print(f"{city:18s} {str(is_nan):12s} {blat:9.3f} {blon:9.3f} {dkm:8.1f}")
        else:
            print(f"{city:18s} {str(is_nan):12s}  NO LAND within +/-{w} deg")


if __name__ == "__main__":
    main()
