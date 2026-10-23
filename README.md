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

This repository contains R code to apply the Kotwicki et al. (2014) density-dependent correction to pollock survey data. The workflow estimates corrected CPUE, biomass, and age/length compositions for eastern Bering Sea pollock using survey haul, catch, length, and specimen data.

---

## Application Structure

The project is organized around a small set of R scripts and supporting directories:

### Core analysis scripts

- `R/Density_Dep_estimates_MAIN.R`  
  Runs the main density-dependent correction workflow.

- `R/plot_ddc_effect.R`  
  Produces diagnostic plots for the density-dependent correction.

- `R/plot_db_index.R`  
  Generates design-based index plots.

- `R/plot_comps.R`  
  Produces composition and size-structure figures.

### Supporting project files

- `R/functions/`  
  Contains helper functions for data processing, CPUE, biomass, age/length composition, and diagnostics.

- `from_stan/`  
  Stores Stan-derived outputs used in the correction process.

- `output/`  
  Stores generated results, tables, and figures.

### General workflow

- Prepare haul, catch, length, and specimen data  
- Apply the density-dependent correction  
- Estimate corrected CPUE and biomass  
- Generate age/length composition summaries  
- Save diagnostic outputs and tables

---

# NOAA README

This repository is a scientific product and is not official communication of the National Oceanic and Atmospheric Administration, or the United States Department of Commerce. All NOAA GitHub project content is provided on an "as is" basis. Any reference to specific commercial products, processes, or services by service mark, trademark, manufacturer, or otherwise, does not constitute or imply their endorsement, recommendation, or favoring by the United States Government, Department of Commerce, or NOAA.

---

# NOAA License

Software code created by U.S. Government employees is not subject to copyright in the United States (17 U.S.C. §105).

The United States/Department of Commerce reserve all rights to seek and obtain copyright protection in countries other than the United States for software authored in its entirety by Department of Commerce employees.

---

<img src="https://raw.githubusercontent.com/nmfs-general-modeling-tools/nmfspalette/main/man/figures/noaa-fisheries-rgb-2line-horizontal-small.png" height="75" alt="NOAA Fisheries">

[U.S. Department of Commerce](https://www.commerce.gov/) |  
[National Oceanographic and Atmospheric Administration](https://www.noaa.gov) |  
[NOAA Fisheries](https://www.fisheries.noaa.gov/)
