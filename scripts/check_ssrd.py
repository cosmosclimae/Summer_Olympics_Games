"""
check_ssrd.py
-------------
Diagnostic: verify that SSRD de-accumulation produces a clean diurnal curve
on ONE real ERA5-Land file, BEFORE running the heavy pipelines.

Picks the first available file, extracts a few cities, de-accumulates SSRD,
and plots the first 72 hours. A correct result shows clean daytime bumps
(zero at night, peak at local noon), no negatives, no saw-teeth.
"""
import os
import glob
import argparse
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import wbgt_core as w

TIME_DIM = "valid_time"

# a few cities at very different longitudes -> peaks should shift in UTC
CHECK_CITIES = {
    "Tokyo":  (35.68, 139.69),   # peak ~03 UTC
    "Paris":  (48.86,   2.35),   # peak ~12 UTC
    "Los Angeles": (34.05, -118.24),  # peak ~20 UTC
}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", required=True)
    ap.add_argument("--domain", default="world")
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    files = sorted(glob.glob(os.path.join(args.data_dir,
                                          f"era5land_{args.domain}_*.nc")))
    if not files:
        raise FileNotFoundError(f"No ERA5-Land files in {args.data_dir}")
    fname = files[0]
    print(f"[check] using {fname}")

    ds = xr.open_dataset(fname)
    print(f"[check] dims: {dict(ds.dims)}")
    print(f"[check] ssrd units: {ds['ssrd'].attrs.get('units')}, "
          f"stepType: {ds['ssrd'].attrs.get('GRIB_stepType')}")

    fig, ax = plt.subplots(figsize=(10, 5))
    for city, (lat, lon) in CHECK_CITIES.items():
        pt = ds.sel(latitude=lat, longitude=lon, method="nearest")
        flux = w.deaccumulate_ssrd(pt["ssrd"], time_dim=TIME_DIM)
        f72 = flux.isel({TIME_DIM: slice(0, 72)})
        ax.plot(np.arange(len(f72)), f72.values, label=city, marker=".")
        print(f"[check] {city}: min={float(flux.min()):.1f} "
              f"max={float(flux.max()):.1f} W/m2")

    ax.set_xlabel("hours from file start (UTC)")
    ax.set_ylabel("de-accumulated SSRD (W/m$^2$)")
    ax.set_title("SSRD de-accumulation check — expect clean diurnal bumps, no negatives")
    ax.axhline(0, color="k", lw=0.5)
    ax.legend()
    ax.grid(alpha=0.3)
    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    fig.tight_layout()
    fig.savefig(args.output, dpi=120)
    print(f"[check] saved {args.output}")
    print("[check] INSPECT THE PLOT: night=0, peak at local noon, no negatives, "
          "no saw-teeth. Tokyo peak should be ~03 UTC, Paris ~12, LA ~20.")


if __name__ == "__main__":
    main()
