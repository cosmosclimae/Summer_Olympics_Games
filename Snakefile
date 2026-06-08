# -*- coding: utf-8 -*-
"""
Snakefile — FIFA North America / ERA5-Land + ERA5
            + WBGT post-processing (CIO Olympics thermal-constraints paper)

Téléchargements horaires + calcul WBGT (outdoor/indoor) pour le papier JO.

Téléchargement (inchangé):
  snakemake -s Snakefile -j 1 --config target=era5land domain=world

Post-traitement WBGT (nouveau):
  # cartes mondiales (Fig 1) — lourd
  snakemake -s Snakefile -j 1 --config target=wbgt_global domain=world
  # 8 villes (Fig 2-5) — léger
  snakemake -s Snakefile -j 1 --config target=wbgt_cities domain=world
  # test diagnostic SSRD sur un fichier
  snakemake -s Snakefile -j 1 --config target=ssrd_check domain=world

Conseil: ne pas dépasser -j 1 ou -j 2 avec CDS, sinon tu vas te faire bloquer ou ralentir.
"""
from pathlib import Path

# =========================
# CONFIG
# =========================

# [N, W, S, E]
AREAS = {
    "north_america_fifa": [60.0, -140.0, 15.0, -50.0],
    "north_america_large": [75.0, -170.0, 10.0, -50.0],
    "world": [85.0, -180.0, -60.0, 180.0],
}

DATA_ROOT = Path(config.get("data_root", "/mnt/f/climae_CIO/")).expanduser()
VENV_ACTIVATE = str(config.get("venv_activate", "~/venv_cds/bin/activate"))
# all | era5land | era5 | wbgt_global | wbgt_cities | ssrd_check
TARGET = str(config.get("target", "all")).lower()
DOMAIN = str(config.get("domain", "north_america_fifa")).lower()

DEFAULT_YEARS = [str(y) for y in range(2006, 2026)]
DEFAULT_MONTHS = ["07", "08"]

YEARS = [str(y) for y in config.get("years", DEFAULT_YEARS)]
MONTHS = [str(m).zfill(2) for m in config.get("months", DEFAULT_MONTHS)]

# Plain scalars for use inside shell blocks (Snakemake cannot index lists there)
YEAR_START = YEARS[0]
YEAR_END = YEARS[-1]

VALID_TARGETS = {"all", "era5land", "era5", "wbgt_global", "wbgt_cities", "ssrd_check"}
if TARGET not in VALID_TARGETS:
    raise ValueError(f"config target must be one of: {' | '.join(sorted(VALID_TARGETS))}")

if DOMAIN not in AREAS:
    raise ValueError(f"config domain must be one of: {' | '.join(AREAS.keys())}")

AREA = AREAS[DOMAIN]

print(f"[CONFIG] TARGET={TARGET}")
print(f"[CONFIG] DOMAIN={DOMAIN}")
print(f"[CONFIG] AREA={AREA}")
print(f"[CONFIG] YEARS={YEARS[0]}->{YEARS[-1]} | MONTHS={MONTHS}")

SCRIPT_ERA5LAND = Path(workflow.basedir) / "scripts" / "download_era5land_fifa_month.py"
SCRIPT_ERA5 = Path(workflow.basedir) / "scripts" / "download_era5_fifa_month.py"

# WBGT post-processing scripts
SCRIPT_WBGT_GLOBAL = Path(workflow.basedir) / "scripts" / "pipeline1_global_maps.py"
SCRIPT_WBGT_CITIES = Path(workflow.basedir) / "scripts" / "pipeline2_cities.py"
SCRIPT_SSRD_CHECK  = Path(workflow.basedir) / "scripts" / "check_ssrd.py"

for script in [SCRIPT_ERA5LAND]:
    if not script.exists():
        raise FileNotFoundError(f"Missing downloader script: {script}")

# =========================
# PATHS
# =========================

RAW = DATA_ROOT / "raw" / DOMAIN
ERA5LAND_DIR = RAW / "era5land"
ERA5_DIR = RAW / "era5"

PROCESSED = DATA_ROOT / "processed" / DOMAIN
PROC_GLOBAL_DIR = PROCESSED / "global"
PROC_CITIES_DIR = PROCESSED / "cities"

ERA5LAND_BASENAME = f"era5land_{DOMAIN}" + "_{year}-{month}.nc"
ERA5_BASENAME = f"era5_{DOMAIN}" + "_{year}-{month}.nc"

ERA5LAND_FILES = expand(
    str(ERA5LAND_DIR / ERA5LAND_BASENAME),
    year=YEARS,
    month=MONTHS,
)

ERA5_FILES = expand(
    str(ERA5_DIR / ERA5_BASENAME),
    year=YEARS,
    month=MONTHS,
)

# WBGT outputs
GLOBAL_MAP_OUT = str(PROC_GLOBAL_DIR / f"global_days_gt28_{YEAR_START}-{YEAR_END}.nc")
CITIES_SUMMARY_OUT = str(PROC_CITIES_DIR / "city_summary.csv")
SSRD_CHECK_OUT = str(PROC_GLOBAL_DIR / "ssrd_check.png")


def targets():
    if TARGET == "era5land":
        return ERA5LAND_FILES
    if TARGET == "era5":
        return ERA5_FILES
    if TARGET == "wbgt_global":
        return [GLOBAL_MAP_OUT]
    if TARGET == "wbgt_cities":
        return [CITIES_SUMMARY_OUT]
    if TARGET == "ssrd_check":
        return [SSRD_CHECK_OUT]
    return ERA5LAND_FILES + ERA5_FILES


rule all:
    input:
        targets()


# =========================
# DOWNLOAD RULES (unchanged)
# =========================

rule download_era5land_month:
    output:
        ERA5LAND_DIR / ERA5LAND_BASENAME
    shell:
        r"""
        bash -lc "set -euo pipefail; \
          mkdir -p '{ERA5LAND_DIR}' && \
          source {VENV_ACTIVATE} && \
          python '{SCRIPT_ERA5LAND}' \
            --domain '{DOMAIN}' \
            --year '{wildcards.year}' \
            --month '{wildcards.month}' \
            --output '{output}'"
        """


# =========================
# WBGT POST-PROCESSING RULES (new)
# =========================

# Diagnostic: de-accumulated SSRD curve on one file (RUN THIS FIRST)
rule ssrd_check:
    input:
        ERA5LAND_FILES
    output:
        SSRD_CHECK_OUT
    shell:
        r"""
        bash -lc "set -euo pipefail; \
          mkdir -p '{PROC_GLOBAL_DIR}' && \
          source {VENV_ACTIVATE} && \
          python '{SCRIPT_SSRD_CHECK}' \
            --data-dir '{ERA5LAND_DIR}' \
            --domain '{DOMAIN}' \
            --output '{output}'"
        """

# Pipeline 1: global maps (Fig 1) — depends on ALL Jul/Aug ERA5-Land files
rule wbgt_global:
    input:
        ERA5LAND_FILES
    output:
        GLOBAL_MAP_OUT
    shell:
        r"""
        bash -lc "set -euo pipefail; \
          mkdir -p '{PROC_GLOBAL_DIR}' && \
          source {VENV_ACTIVATE} && \
          python '{SCRIPT_WBGT_GLOBAL}' \
            --data-dir '{ERA5LAND_DIR}' \
            --domain '{DOMAIN}' \
            --years '{YEAR_START}' '{YEAR_END}' \
            --output '{output}'"
        """

# Pipeline 2: 8 cities (Fig 2-5) — depends on the same files
rule wbgt_cities:
    input:
        ERA5LAND_FILES
    output:
        CITIES_SUMMARY_OUT
    shell:
        r"""
        bash -lc "set -euo pipefail; \
          mkdir -p '{PROC_CITIES_DIR}' && \
          source {VENV_ACTIVATE} && \
          python '{SCRIPT_WBGT_CITIES}' \
            --data-dir '{ERA5LAND_DIR}' \
            --domain '{DOMAIN}' \
            --years '{YEAR_START}' '{YEAR_END}' \
            --out-dir '{PROC_CITIES_DIR}'"
        """
