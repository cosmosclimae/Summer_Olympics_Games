import argparse
import calendar
from pathlib import Path

import cdsapi

# [N, W, S, E]
AREAS = {
    # Domaine plus raisonnable pour FIFA 2026 / USA + sud Canada + Mexique.
    "north_america_fifa": [60.0, -140.0, 15.0, -50.0],
    # Domaine large si tu veux garder Alaska / nord Canada / contexte complet.
    "north_america_large": [75.0, -170.0, 10.0, -50.0],
    "world": [85.0, -180.0, -60.0, 180.0],
}

VARS_ERA5LAND = [
    "2m_temperature",
    "2m_dewpoint_temperature",
    "10m_u_component_of_wind",
    "10m_v_component_of_wind",
    "surface_solar_radiation_downwards",
]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", choices=sorted(AREAS.keys()), required=True)
    ap.add_argument("--year", required=True)
    ap.add_argument("--month", required=True)
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)

    ndays = calendar.monthrange(int(args.year), int(args.month))[1]

    request = {
        "variable": VARS_ERA5LAND,
        "year": args.year,
        "month": args.month,
        "day": [f"{d:02d}" for d in range(1, ndays + 1)],
        "time": [f"{h:02d}:00" for h in range(24)],
        "area": AREAS[args.domain],
        "data_format": "netcdf",
        "download_format": "unarchived",
    }

    print(
        f"[CDS] ERA5-Land download {args.domain} {args.year}-{args.month} -> {out}"
    )
    client = cdsapi.Client()
    client.retrieve("reanalysis-era5-land", request, str(out))


if __name__ == "__main__":
    main()
