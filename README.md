# pollock-ddc

> Code is always in development  
>  
> This document is for informational purposes only and does not necessarily represent the views or official position of the Department of Commerce, the National Oceanic and Atmospheric Administration, or the National Marine Fisheries Service.

---

## Maintainers

This code is primarily maintained by:

- **Sophia Wassermann**

Alaska Fisheries Science Center  
National Marine Fisheries Service  
National Oceanic and Atmospheric Administration  
Seattle, WA 98195

---

## Overview

This repository contains R code to apply a density-dependent correction (DDC) to pollock survey data based on the Kotwicki et al. (2014) method. The workflow combines survey haul data, catch/length/specimen information, and model-based output to estimate corrected CPUE, biomass, and age/length-composition summaries for eastern Bering Sea pollock.

The analysis focuses on correcting survey indices for density dependence in the relationship between backscatter and fish density, and then producing downstream outputs for assessment support and diagnostic plotting.

---

## Application Structure

The project is organized around a small set of R scripts and supporting directories:

### Core analysis scripts

- `R/Density_Dep_estimates_MAIN.R`  
  Runs the main density-dependent correction workflow, including data processing, corrected CPUE/biomass estimation, and output generation.

- `R/plot_ddc_effect.R`  
  Produces plots illustrating the density-dependent correction effect and related diagnostics.

- `R/plot_db_index.R`  
  Generates plots for design-based index summaries and comparisons.

- `R/plot_comps.R`  
  Produces composition and size-structure diagnostic figures.

- `R/drive_upload.R`  
  Provides utilities for uploading generated outputs to Google Drive.

### Supporting project files

- `R/functions/`  
  Contains the helper functions used throughout the workflow for data retrieval, CPUE calculations, population estimates, bootstrap routines, and output preparation.

- `from_stan/`  
  Stores external Stan-derived results used in the correction process, including `abc_exp.csv`.

- `output/`  
  Folder for generated model outputs, summaries, plots, and supporting tables.

- `pollock-ddc.Rproj`  
  RStudio project file for working within the repository environment.

### General workflow

- Prepare survey haul, catch, length, and specimen data  
- Apply the density-dependent correction to CPUE and biomass estimates  
- Derive corrected length and age composition summaries  
- Generate diagnostics and comparison plots  
- Save outputs for assessment use and reproducibility

---

# NOAA README

This repository is a scientific product and is not official communication of the National Oceanic and Atmospheric Administration, or the United States Department of Commerce. All NOAA GitHub project content is provided on an 'as is' basis and the user assumes responsibility for its use. Any reference to specific commercial products, processes, or services by service mark, trademark, manufacturer, or otherwise, does not constitute or imply their endorsement, recommendation, or favoring by the United States Government, Department of Commerce, or NOAA.

---

# NOAA License

Software code created by U.S. Government employees is not subject to copyright in the United States (17 U.S.C. §105).

The United States/Department of Commerce reserve all rights to seek and obtain copyright protection in countries other than the United States for Software authored in its entirety by the Department of Commerce employees.

---

<img src="https://raw.githubusercontent.com/nmfs-general-modeling-tools/nmfspalette/main/man/figures/noaa-fisheries-rgb-2line-horizontal-small.png" height="75" alt="NOAA Fisheries">

[U.S. Department of Commerce](https://www.commerce.gov/) |  
[National Oceanographic and Atmospheric Administration](https://www.noaa.gov) |  
[NOAA Fisheries](https://www.fisheries.noaa.gov/)
