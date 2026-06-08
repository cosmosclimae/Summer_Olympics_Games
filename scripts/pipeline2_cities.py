"""
pipeline2_cities.py
-------------------
Pipeline 2 -- CITY SCALE for Figures 2-5.

Extracts hourly WBGT (outdoor + indoor) at 8 cities over Jul-Aug 2006-2025
and builds:
  Fig 2 : per-city days with >=1h >28 and >=1h >32  (outdoor + indoor)
  Fig 3 : diurnal profiles, % of local hours exceeding 28 / 32
  Fig 4 : decomposition of peak-hour outdoor WBGT (T / humidity / radiation share)
  Fig 5 : per-year (per-August) constrained-day count + linear trend

Local time: UTC -> mean local solar time via longitude offset (hours = lon/15).

Output: a tidy NetCDF / CSVs with all city-level series.
"""

import os
import argparse
import numpy as np
import pandas as pd
import xarray as xr
import wbgt_core as w

MONTHS   = [7, 8]
TIME_DIM = "valid_time"

# ----------------------------------------------------------------------
# Eight study cities  [VERIFY coordinates: city-centre vs venue cluster]
# lon in -180..180 to match the grid
# ----------------------------------------------------------------------
CITIES = {
    # Rio & Cape Town nudged to nearest ERA5-Land LAND cell (0.1deg grid masks
    # the exact city centre as sea); land cells verified via diag_landmask.py.
    "Rio de Janeiro": (-22.90, -43.10),   # was -22.91,-43.20 (sea in ERA5-Land)
    "Tokyo":          ( 35.68, 139.69),
    "Paris":          ( 48.86,   2.35),
    "Los Angeles":    ( 34.05, -118.24),
    "Brisbane":       (-27.47, 153.03),
    "Ahmedabad":      ( 23.03,  72.58),
    "Doha":           ( 25.29,  51.53),
    "Cape Town":      (-33.90,  18.50),   # was -33.92,18.42 (sea in ERA5-Land)
}

# Olympic session windows (local solar hour) for Fig 3 shading
SESSIONS = {"morning": (9, 12), "evening": (19, 22)}


def local_hour_offset(lon):
    """Mean solar time offset in (fractional) hours from UTC."""
    return lon / 15.0


def extract_city(lat, lon, data_dir, domain, years):
    """Concatenated hourly series (all Jul-Aug, all years) for one city,
    with WBGT outdoor/indoor and the components needed for decomposition."""
    recs = []
    nan_report = []
    for yr in years:
        for m in MONTHS:
            fname = os.path.join(data_dir, f"era5land_{domain}_{yr}-{m:02d}.nc")
            ds = xr.open_dataset(fname)
            pt = ds.sel(latitude=lat, longitude=lon, method="nearest")

            # ERA5-Land masks sea cells as NaN (the mask is time-invariant, so
            # test the first hour). If the nearest cell is sea, search a window
            # for the closest valid land cell.
            is_sea = bool(np.isnan(pt["t2m"].isel({TIME_DIM: 0})))
            if is_sea:
                wdeg = 0.6
                win = ds.sel(
                    latitude=slice(lat + wdeg, lat - wdeg),   # lat descending
                    longitude=slice(lon - wdeg, lon + wdeg),
                )
                land = win["t2m"].isel({TIME_DIM: 0}).notnull()
                if bool(land.any()):
                    la = win["latitude"].values
                    lo = win["longitude"].values
                    LO, LA = np.meshgrid(lo, la)
                    dist = np.sqrt((LA - lat) ** 2 + (LO - lon) ** 2)
                    dist[~land.values] = np.inf
                    j, i = np.unravel_index(np.argmin(dist), dist.shape)
                    pt = win.sel(latitude=la[j], longitude=lo[i],
                                 method="nearest")

            frac_nan = float(np.isnan(pt["t2m"]).mean())
            nan_report.append(frac_nan)

            ssrd_flux = w.deaccumulate_ssrd(pt["ssrd"], time_dim=TIME_DIM)
            Ta_C = pt["t2m"] - 273.15
            RH   = w.relative_humidity(pt["t2m"], pt["d2m"])
            Tw   = w.wet_bulb_stull(Ta_C, RH)
            Tg_sun   = w.globe_hajizadeh(ssrd_flux, Ta_C, RH)
            Tg_shade = w.globe_hajizadeh(xr.zeros_like(ssrd_flux), Ta_C, RH)
            wb_out = 0.7 * Tw + 0.2 * Tg_sun   + 0.1 * Ta_C
            wb_in  = 0.7 * Tw + 0.3 * Tg_shade

            t = pt[TIME_DIM].values
            df = pd.DataFrame({
                "time_utc": t,
                "Ta": Ta_C.values,
                "RH": RH.values,
                "ssrd": ssrd_flux.values,
                "Tw": Tw.values,
                "Tg_sun": Tg_sun.values,
                "wbgt_out": wb_out.values,
                "wbgt_in": wb_in.values,
            })
            recs.append(df)
            ds.close()

    out = pd.concat(recs, ignore_index=True)
    out["time_utc"] = pd.to_datetime(out["time_utc"])
    # local solar time
    off = local_hour_offset(lon)
    out["time_local"] = out["time_utc"] + pd.to_timedelta(off, unit="h")
    out["local_hour"] = out["time_local"].dt.hour
    out["date_local"] = out["time_local"].dt.floor("D")
    out["year"] = out["time_local"].dt.year

    mean_nan = float(np.mean(nan_report)) if nan_report else 0.0
    if mean_nan > 0.05:
        print(f"    [WARN] mean NaN fraction at this city = {mean_nan:.1%} "
              f"-- check coordinates (sea cell?)", flush=True)
    return out


def indicators_for_city(df):
    """Compute the per-city indicators from the hourly dataframe."""
    res = {}

    # Fig 2: days with >=1h above thresholds (outdoor & indoor), mean per year
    def days_with_exc(col, thr):
        daily = df.groupby("date_local")[col].apply(lambda s: (s > thr).any())
        n_years = df["year"].nunique()
        return daily.sum() / n_years
    res["days_out_28"] = days_with_exc("wbgt_out", 28)
    res["days_out_32"] = days_with_exc("wbgt_out", 32)
    res["days_in_28"]  = days_with_exc("wbgt_in", 28)
    res["days_in_32"]  = days_with_exc("wbgt_in", 32)

    # Fig 3: diurnal profile, % of hours exceeding 28 / 32 by local hour
    prof = []
    for h in range(24):
        sub = df[df["local_hour"] == h]
        prof.append({
            "local_hour": h,
            "pct_out_28": 100 * (sub["wbgt_out"] > 28).mean(),
            "pct_out_32": 100 * (sub["wbgt_out"] > 32).mean(),
            "pct_in_28":  100 * (sub["wbgt_in"] > 28).mean(),
            "pct_in_32":  100 * (sub["wbgt_in"] > 32).mean(),
        })
    res["diurnal"] = pd.DataFrame(prof)

    # Fig 4: decomposition of peak-hour outdoor WBGT into term contributions
    # contributions: 0.7*Tw (humidity/temp), 0.2*Tg_sun (radiation+temp), 0.1*Ta
    # Drop hours/days where wbgt_out is NaN (ERA5-Land sea-mask cells) so that
    # idxmax does not hit all-NaN groups.
    df_valid = df.dropna(subset=["wbgt_out"])
    if len(df_valid) == 0:
        res["decomp"] = {k: np.nan for k in
                         ["term_Tw", "term_Tg", "term_Ta", "mean_peak_wbgt",
                          "mean_RH_at_peak", "mean_ssrd_at_peak"]}
    else:
        idx = df_valid.groupby("date_local")["wbgt_out"].idxmax()
        peak = df_valid.loc[idx]
        res["decomp"] = {
            "term_Tw":  (0.7 * peak["Tw"]).mean(),
            "term_Tg":  (0.2 * peak["Tg_sun"]).mean(),
            "term_Ta":  (0.1 * peak["Ta"]).mean(),
            "mean_peak_wbgt": peak["wbgt_out"].mean(),
            "mean_RH_at_peak": peak["RH"].mean(),
            "mean_ssrd_at_peak": peak["ssrd"].mean(),
        }

    # Fig 5: per-year constrained-day count (>=1h >28 outdoor) + linear trend
    yearly = (df.groupby(["year", "date_local"])["wbgt_out"]
                .apply(lambda s: (s > 28).any())
                .groupby("year").sum())
    years = yearly.index.values.astype(float)
    vals  = yearly.values.astype(float)
    if len(years) > 1:
        slope, intercept = np.polyfit(years, vals, 1)
        # 95% CI on slope
        n = len(years)
        yhat = slope * years + intercept
        resid = vals - yhat
        se = np.sqrt(np.sum(resid**2) / (n - 2)) / np.sqrt(np.sum((years - years.mean())**2))
        ci95 = 1.96 * se
    else:
        slope, intercept, ci95 = np.nan, np.nan, np.nan
    res["trend"] = {"slope_days_per_yr": slope, "intercept": intercept,
                    "ci95": ci95}
    res["yearly"] = yearly.reset_index().rename(columns={"wbgt_out": "constrained_days"})

    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", required=True)
    ap.add_argument("--domain", default="world")
    ap.add_argument("--years", nargs=2, type=int, default=[2006, 2025],
                    metavar=("START", "END"))
    ap.add_argument("--out-dir", required=True)
    args = ap.parse_args()

    years = range(args.years[0], args.years[1] + 1)
    out_dir = args.out_dir
    os.makedirs(out_dir, exist_ok=True)

    summary_rows = []
    for city, (lat, lon) in CITIES.items():
        print(f"  extracting {city} ...", flush=True)
        df = extract_city(lat, lon, args.data_dir, args.domain, years)
        df.to_parquet(os.path.join(out_dir, f"hourly_{city.replace(' ','_')}.parquet"))

        ind = indicators_for_city(df)
        ind["diurnal"].to_csv(
            os.path.join(out_dir, f"diurnal_{city.replace(' ','_')}.csv"), index=False)
        ind["yearly"].to_csv(
            os.path.join(out_dir, f"yearly_{city.replace(' ','_')}.csv"), index=False)

        summary_rows.append({
            "city": city, "lat": lat, "lon": lon,
            **{k: ind[k] for k in ["days_out_28","days_out_32","days_in_28","days_in_32"]},
            **{f"decomp_{k}": v for k, v in ind["decomp"].items()},
            **{f"trend_{k}": v for k, v in ind["trend"].items()},
        })

    pd.DataFrame(summary_rows).to_csv(
        os.path.join(out_dir, "city_summary.csv"), index=False)
    print(f"Saved city summary to {out_dir}/city_summary.csv")


if __name__ == "__main__":
    main()
