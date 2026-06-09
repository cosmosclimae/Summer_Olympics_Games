# Olympic Summer Heat Risk Pipeline

This repository contains the data-processing workflow used to assess heat-related climate risks for potential Summer Olympic Games host cities. The pipeline computes sport-relevant thermal exposure indicators from hourly ERA5-Land reanalysis data and extracts them for selected candidate cities, venues, or urban locations during the typical Summer Olympic period.

The workflow is designed to support the study of how present-day summer climate conditions constrain the feasibility, safety, and scheduling of outdoor Olympic events. It focuses on the interaction between heat, humidity, time of day, and event location, with particular attention to the implications for athletes, spectators, staff, and event operations.

## Scientific objective

The objective of this pipeline is to quantify historical heat exposure during the July–August Olympic window and to compare climatic constraints across potential host cities. The workflow allows the calculation of hourly and daily thermal indicators, their aggregation over relevant time windows, and their extraction at Olympic-relevant locations.

The pipeline supports analyses such as:

* comparison of heat exposure between candidate host cities;
* identification of high-risk hours for outdoor events;
* estimation of the frequency of hazardous thermal conditions during the Olympic period;
* assessment of whether morning or evening scheduling reduces thermal risk;
* production of maps, tables, and figures for scientific publication.

## Input data

The workflow uses hourly ERA5-Land reanalysis data for July and August over the study period. Required variables include, depending on the indicator computed:

* 2 m air temperature;
* 2 m dew point temperature or relative humidity;
* 10 m wind components;
* surface solar radiation downward;

Input files are expected to be provided as NetCDF files organised by year and month. The original ERA5-Land data are not redistributed in this repository.

## Main indicators

The pipeline computes heat-stress indicators relevant to outdoor sport and event planning, including:

* Wet-Bulb Globe Temperature (WBGT), estimated from hourly meteorological variables;
* air temperature thresholds;
* humidity-related heat exposure metrics;
* hourly exceedance counts above sport-relevant thresholds;
* aggregated exposure statistics for the Olympic period.

Additional indicators can be added modularly depending on the needs of the analysis.

## Workflow

The pipeline is structured to process the data in reproducible steps:

1. prepare and harmonise hourly ERA5-Land input files;
2. compute derived meteorological variables where needed;
3. calculate heat-stress indicators at hourly resolution;
4. aggregate indicators by day, hour, city, and Olympic time window;
5. extract values for selected candidate cities or venues;
6. generate summary tables and figure-ready outputs.

The workflow is intended to be run year by year to reduce memory usage and facilitate reproducibility on local machines or HPC environments.

## Repository structure

```text
.
├── scripts/             # Processing and analysis scripts
├── SNAKEworkflow_file            # Snakemake workflow files, if applicable
├── results/ on Zenodo             # Derived outputs, tables and figure-ready files
└── README.md            # Repository documentation
```

## Reproducibility

The repository provides the scripts and workflow required to reproduce the processed indicators, summary tables, and figures used in the associated manuscript. Large input climate files are not included and must be downloaded separately from the original ERA5-Land source.

Users should update the configuration files to specify:

* local paths to ERA5-Land input files;
* the study period;
* the list of candidate cities or venues;
* the target indicators and thresholds;
* the output directory.

## Associated manuscript

This pipeline is associated with a manuscript assessing heat-related climate constraints for Summer Olympic Games host cities. The repository is intended to provide transparency on the data-processing workflow and to support reproducibility of the reported results.

## Citation

If you use or adapt this workflow, please cite the associated manuscript and the archived software release once available.

## License

The code in this repository is released under an open-source license. The ERA5-Land input data remain subject to the terms and conditions of the original data provider.

