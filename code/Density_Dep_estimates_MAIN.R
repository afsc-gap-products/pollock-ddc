# This code produces density-dependent estimates for the Bering Sea Pollock assessment
# This code is converted from SQL files by Stan Kotwicki
# Author: Caitlin I. Allen Akselrud
# Maintained by: Sophia N. Wassermann
# Contact: sophia.wassermann@noaa.gov
# Date created: 2020.08.24
# Date updated: 2024.12.04

# notes -------------------------------------------------------------------

# Legacy code information: Stan's origianl SQL and R code:
#   (fpc = fishing power correction; invalid method)
#   tables generated with fpc have ending '_c'; DO NOT USE
# tables for index corrected density dependence have ending '_n'

# there are sections in this code that ask for a user input. 
#   # If you run code from the top, it will run default settings
#   # Please PAY ATTENTION to ensure you are using the setting you need

# age-length keys (alk) methods may need updating in 2023 assessment
#     # Correa et al 2020 CJFAS
#     # GitHub: https://github.com/gmoroncorrea/STageCompsEstimation 

# ASSESSMENT TABLES -------------------------------------------------------

#   Tables Jim Ianelli needs (all EBS + NBS):
#     - CPUE (number of fish) - uncorrected
#     - CPUE (kg/hectare)     - uncorrected
#     - CPUE (number of fish) - dd corrected
#     - CPUE (kg/hectare)     - dd corrected
#     - CPUE (number of fish, by length) - uncorrected
#     - CPUE (kg/hectare, by length)     - uncorrected
# 
#     - Biomass (thousands of tons) - uncorrected
#     - Biomass (thousands of tons) - dd corrected
#     - Biomass (thousands of tons) just EBS - dd corrected
#     - Biomass (thousands of tons) just NBS - dd corrected
# 
#     - Co-variance matrix (on biomass, thousands of tons) - dd corrected
# 
#     - Length compositions - dd corrected

#     PREP for VAST:
#      - biomass for EBS + NBS together  - dd correction
#      - biomass for just EBS            - dd correction
#      - biomass for just NBS            - dd correction
#      - age comps for EBS + NBS together- dd correction

#     OTHER INFO:
#       - number of stations (by year)
#       - number of fish measured (by year)
#       - number of fish aged (by year)


#   -> TABLES FROM VAST:
#      - VAST biomass for EBS + NBS together  - dd correction
#      - VAST biomass for just EBS            - dd correction
#      - VAST biomass for just NBS            - dd correction
#      - VAST covariance matrix for EBS + NBS together - dd correction
#      - VAST covariance matrix for just EBS  - dd correction
#      - VAST covariance matrix for just NBS  - dd correction
#      - Age compositions for EBS + NBS together - dd correction


# Dependencies ------------------------------------------------------------
library(here)
library(tidyverse)
library(RODBC)
library(janitor)
library(mgcv)

# library(magrittr)
# library(lubridate)
# library(modelr)
# library(googledrive)

functions <- list.files(here("functions"))
walk(functions, ~ source(here("functions", .x)))

# Connect to database -----------------------------------------------------
# Connection established either through saved username and password, or by entering directly
# make sure you are connected to VPN first (and have an AFSC Oracle login)
if (file.exists("Z:/Projects/ConnectToOracle.R")) {
  source("Z:/Projects/ConnectToOracle.R")
} else {
  # For those without a ConnectToOracle file
  channel <- odbcConnect(dsn = "AFSC", 
                         uid = rstudioapi::showPrompt(title = "Username", 
                                                      message = "Oracle Username", 
                                                      default = ""), 
                         pwd = rstudioapi::askForPassword("Enter Password"),
                         believeNRows = FALSE)
}

## checks to see if connection has been established
odbcGetInfo(channel)

# Season-specific fixed inputs --------------------------------------------
# May need to change the vessels and whether the NBS should be included, here.

current_year <- year(Sys.Date())
set_inputs <- function(vessel1 = 162, vessel2 = 134, include_NBS = FALSE) {
  if(include_NBS == TRUE) {
    cruise <- paste0(current_year, "01", ",", current_year, "02")
  }
  
  if(include_NBS == FALSE) {
    cruise <- paste0(current_year, "01")
  }

  vessel_code <- paste0(vessel1, ",", vessel2)  # List each vessel, separated by commas
  vessel_nums <- c(vessel1, vessel2)
  
  # get cruise id from here
  query_command <- paste0(" select * from race_data.v_cruises where cruise in (", cruise,");")
  cruise_info <- sqlQuery(channel, query_command) %>% 
    as_tibble() %>% 
    clean_names()
  
  cruise_id_nums <- cruise_info %>% filter(vessel_id %in% vessel_nums) 
  
  # create cruise ID object - when no NBS, avoid NAs from non-existent cruise
  if(include_NBS == TRUE) {
    cruise_id <- paste0(cruise_id_nums$cruise_id[1], "," , cruise_id_nums$cruise_id[2],
                        "," , cruise_id_nums$cruise_id[3], "," , cruise_id_nums$cruise_id[4])
  }
  
  if(include_NBS == FALSE) {
    cruise_id <- paste0(cruise_id_nums$cruise_id[1], "," , cruise_id_nums$cruise_id[2])
  }
  
  return(list(cruise = cruise, vessel_code = vessel_code, cruise_id = cruise_id))
}

inputs <- set_inputs()

# Define inputs required below
cruise <- inputs$cruise
vessel_code <- inputs$vessel_code
cruise_id <- inputs$cruise_id

# NBS subarea stratum-- this shouldn't change too much, but is a fixed input
NBS_subarea <- c(81, 70, 71, 99) # NBS stratum numbers; added 99 to indicate 2018 NBS emergency survey; diff survey methods

# Strata metadata year; 2022 is the latest update (use for current assessments)
strat_meta_year <- 2022

# Set output - model- or design-based
data_type <- "db"

# Set up folder 
dir_thisyr <- paste0(current_year,"_", data_type, "_data_", strat_meta_year, "_strata")
dir.create(here("output",dir_thisyr))

# data --------------------------------------------------------------------
process_data <- function(first_run = FALSE, estimate_ages = FALSE, save_data = FALSE) {
  # Don't include slope survey for VAST
  slope_survey <- slope_survey_d()
  
  #' Once per year, when the new survey data is in, you need to update the 
  #' "valid hauls". Jason Conner put code in the Oracle safe folder to do this. 
  #' If you do not run this, you will not get the most current year of survey 
  #' data, but you only have to run it once after the survey data are finalized 
  #' (but it doesn't hurt anything if you run it again).
  if(first_run == TRUE) {
    RODBC::sqlQuery(channel, "BEGIN safe.UPDATE_SURVEY; END;")
  }
  
  ## Haul data ----------------------------------------------------------------
  #' Get the haul data. You pull different hauls depending on whether you're 
  #' running the model-based or design-based indices/comps. Design-based indices
  #' do not include the NBS prior to 2010 and do not include the non-standard
  #' 2018 NBS survey; model-based indices include all valid NBS sations 
  #' throughout the time series, including 2018.
  get_hauls <- haul_data_d(data_selection = data_type, 
                           nbs_subarea = NBS_subarea, 
                           slope_info = slope_survey)
  
  if(data_type == "db") {
    hauls_survey_ebs <- get_hauls$good_hauls_DB_EBS
    hauls_survey_nbs <- get_hauls$good_hauls_DB_NBS
    hauls_survey_bad_ebs <- get_hauls$bad_hauls_DB_ebs
    hauls_survey_bad_nbs <- get_hauls$bad_hauls_DB_nbs
    all_hauljoins <- c(hauls_survey_ebs$hauljoin, hauls_survey_bad_ebs$hauljoin)
    all_hauljoins_nbs <- c(hauls_survey_nbs$hauljoin, hauls_survey_bad_nbs$hauljoin)
    valid_hauljoins_nbs <- hauls_survey_nbs$hauljoin
    hauls_survey <- hauls_survey_ebs %>% bind_rows(hauls_survey_nbs)
    hauls_survey_bad <- hauls_survey_bad_ebs %>% bind_rows(hauls_survey_bad_nbs)
  }
  if(data_type == "mb") {
    hauls_survey <- get_hauls$good_hauls_MB
    hauls_survey_bad <- get_hauls$bad_hauls
    all_hauljoins <- c(hauls_survey$hauljoin, hauls_survey_bad$hauljoin)
  } 
  
  valid_hauljoins <- hauls_survey$hauljoin
  
  ## Specimen data ------------------------------------------------------------
  pollock_specimen <- specimen_data_d(hauls_survey_dat = hauls_survey)
  
  # some bad hauls have useful age and length info that we use in the ALK, but not the age comps
  pollock_specimen_bad <- specimen_data_d(hauls_survey_dat = hauls_survey_bad)
  
  pollock_specimen_all <- pollock_specimen %>% 
    bind_rows(pollock_specimen_bad) 
  
  #' Remove NAs if not estimating ages for the current year, using the ALK. 
  #' Because we usually have pollock ages for the current year by the time we 
  #' run the comp code, we don't need to estimate ages, so this is the default.
  if(estimate_ages == FALSE) {
    pollock_specimen_all <- pollock_specimen_all %>% dplyr::filter(!is.na(age))
  }
  
  # Separate EBS and NBS specimen data for design-based age comps
  if(data_type == "db") {
    pollock_specimen_nbs <- specimen_data_d(hauls_survey_dat = hauls_survey_nbs)
    pollock_specimen_bad_nbs <- specimen_data_d(hauls_survey_dat = hauls_survey_bad_nbs)
    pollock_specimen_all_nbs <- pollock_specimen_nbs %>%
      bind_rows(pollock_specimen_bad_nbs) %>%
      filter(!is.na(age))
    
    pollock_specimen_ebs <- pollock_specimen %>% filter(!stratum %in% NBS_subarea)
  }
  
  ## Catch data ---------------------------------------------------------------
  pollock_catch <- catch_data_d(hauljoins = valid_hauljoins)
  
  ## Length data --------------------------------------------------------------
  pollock_length_info <- length_data_d(valid_hauljoins)
  
  pollock_length <- pollock_length_info$pollock_length
  pollock_raw_length <- pollock_length_info$pollock_raw_length
  
  # Separate EBS and NBS length data for design-based age comps
  if(data_type == "db") {
    pollock_length_nbs <- length_data_d(hauljoins = valid_hauljoins_nbs)
    pollock_raw_length_nbs <- pollock_length_nbs$pollock_raw_length
    pollock_raw_length_ebs <- left_join(pollock_raw_length, hauls_survey) %>% 
      filter(!stratum %in% NBS_subarea)
  }
  
  ## Metadata (aka update cruise) ---------------------------------------------
  strata_metadata <- metadata_d(cruise_id = cruise, 
                                vessel_code_id = vessel_code, 
                                pollock_specimen_data = pollock_specimen, 
                                meta_select = strat_meta_year)
  
  # Save raw data -------------------------------------------------------------
  if(save_data == TRUE) {
    raw_data_save <- save_raw_data(hauls_survey = hauls_survey, 
                                   pollock_specimen = pollock_specimen, 
                                   pollock_catch = pollock_catch,
                                   pollock_length = pollock_length,
                                   # ebs_shelf_cruise = ebs_shelf_cruise, 
                                   metadata = strata_metadata) #, 
                                   # strata_metadata_raw = strata_metadata_raw)
  }
  
  ## Data Processing ----------------------------------------------------------
  # Calculate CPUE 
  cpue_info <- get_cpue(hauls = hauls_survey, 
                        catch = pollock_catch, 
                        strata_filter = strata_metadata)
  
  # Convert CPUE to population (numbers) and biomass using stratum areas
  all_strata <- get_stratum_proportions(ebs_strata = strata_metadata)
  pop_info <- population(ebs_strata = all_strata, 
                         avg_cpue = cpue_info$avg_cpue_yr_strat)
  
  # Generate length comps
  length_comps <- length_comp_f(length = pollock_length, 
                                hauls = hauls_survey, 
                                catch = pollock_catch, 
                                stratum = pop_info$ebs_strata, 
                                biomass_kg = pop_info$pollock_biomass_kg_ha, 
                                biomass_mt = pop_info$pollock_biomass_MT_ha)
  
  # Put length comps into a table
  cpue_length_table <- full_join(length_comps$sizecomp_cpue_stn, hauls_survey) %>%  #, by = c("vessel", "year", "haul"))
    select(cruisejoin, vessel, year, haul, start_latitude, start_longitude, 
           stratum, stationid, sex, length, cpue_length_num_ha)
  
  #' NB: the age length key and age comp names were reversed in Stan's code, 
  #' which is reflected in the function names here.
  age_comps <- get_agecomps(specimen = pollock_specimen_all, 
                            pop_lengths = length_comps$pollock_length_comp)
  
  # Extract the full key for use
  age_comp_full_key <- bind_rows(age_comps$aal_filled_key) 
  
  # CIA: pollock_specimen vs pollock_specimen_all here?? (think about whether bad hauls should be included)
  al_key <- get_al_key(pollock_specimen_all, length_comps$pollock_length_comp)
  # CIA: NOTE- consider updating alk methods for 2022
  
  ### Convert new 1 -----------------------------------------------------------
  # this is the conversion of stan's "new1.R' which starts the density-dependent conversion process
  # requres: pollock_cpue_length = cpue_length_table
  # pollock_cpue = cpue_info$avg_cpue_yr_strat
  # z csxzaazaserraqwdfdhghggf tteszzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz <- ELI IS HELPING
  equiv_bs_length <- backscatter_est(cpue_length_table)
  
  ### Convert new 2 -----------------------------------------------------------
  #' this is the conversion of Stan's "new2.R' which continues the DDC process by getting the backscatter correction for CPUE
  # read in abc_exp.csv table
  param_ests <-  read_csv(here("from_stan", "abc_exp.csv"))
  # asp, slope, and intercept = 3 parameters of the correction fxn; generated in other code; 
  # saved as is until new data (future)
  backscatter_cpue_corr <- ddc_cpue_bs(equiv_bs_length, param_ests)
  
  ### Convert new 3 -----------------------------------------------------------
  # this is the conversion of Stan's 'new3.R' which applies the correction and builds the ddc cpue table
  ddc_table <- ddc_fxn(backscatter_cpue_corr, cpue_info$p_cpue)
  
  check_ddc <- ddc_table %>% 
    select(vessel, year, haul, cruise, stratum, start_latitude, start_longitude, 
           species_code, ddc_cpue_kg_ha, ddc_cpue_num_ha, area_fished_ha) %>% 
    rename(latitude = start_latitude, longitude = start_longitude) %>% 
    select(- latitude, -longitude, -species_code)
  
  ## Return tables for next steps ---------------------------------------------
  out <- list(slope_survey = slope_survey,
              hauls_survey = hauls_survey, 
              pollock_specimen = pollock_specimen,
              pollock_catch = pollock_catch,
              pollock_length = pollock_length,
              all_strata = all_strata,
              pop_info = pop_info,
              strata_metadata = strata_metadata,
              cpue_length_table = cpue_length_table,
              age_comp_full_key = age_comp_full_key,
              ddc_table = ddc_table, 
              cpue_info = cpue_info)
  
  if(data_type == "mb") {return(out)}
  
  # Include additional tables generated in the design-based section
  if(data_type == "db") {
    extra <- list(hauls_survey_ebs = hauls_survey_ebs, 
                  hauls_survey_nbs = hauls_survey_nbs,
                  pollock_specimen_ebs = pollock_specimen_ebs,
                  pollock_specimen_nbs = pollock_specimen_nbs,
                  pollock_raw_length_ebs = pollock_raw_length_ebs,
                  pollock_raw_length_nbs = pollock_raw_length_nbs)
    out2 <- c(out, extra)
    return(out2)
  }
}

tables <- process_data()

# Define inputs for next function - save to environment variables explicitly
hauls_survey <- tables$hauls_survey
pollock_specimen <- tables$pollock_specimen
pollock_catch <- tables$pollock_catch
pollock_length <- tables$pollock_length
all_strata <- tables$all_strata
strata_metadata <- tables$strata_metadata
ddc_table <- tables$ddc_table


# Density-dependent correction ------------------------------------------------
ddc_conversion <- function() {
  # CPUE
  ddc_cpue <- ddc_table %>%
    select(-cpue_num_ha, -cpue_kg_ha) %>% 
    rename(cpue_num_ha = ddc_cpue_num_ha,
           cpue_kg_ha = ddc_cpue_kg_ha) %>% 
    ungroup() %>% 
    group_by(year, stratum) %>% 
    summarize(number_hauls_no = length(cpue_num_ha[!is.na(cpue_num_ha)]), 
              number_hauls_kg = length(cpue_kg_ha[!is.na(cpue_kg_ha )]),
              avg_cpue_no_ha = mean(cpue_num_ha, na.rm = TRUE), 
              var_cpue_no_ha = var(cpue_num_ha,  na.rm = TRUE), 
              avg_cpue_kg_ha = mean(cpue_kg_ha,  na.rm = TRUE), 
              var_cpue_kg_ha = var(cpue_kg_ha,   na.rm = TRUE), 
              number_hauls = length(haul),
              tot_area_fished = sum(area_fished_ha, na.rm = TRUE))
  
  # Population
  ddc_pop_ests <- population(ebs_strata = all_strata, avg_cpue = ddc_cpue)
  ddc_pop_ests_EBS <- ddc_pop_ests$pollock_biomass_MT_ha %>% dplyr::filter(subarea != 0)
  ddc_pop_ests_NBS <- ddc_pop_ests$pollock_biomass_MT_ha %>% dplyr::filter(subarea == 0)
  
  # Length comps
  ddc_length_comps <- length_comp_f(pollock_length, hauls_survey, pollock_catch, 
                                    ddc_pop_ests$ebs_strata, 
                                    ddc_pop_ests$pollock_biomass_kg_ha, 
                                    ddc_pop_ests$pollock_biomass_MT_ha)
  ddc_cpue_length_table <- full_join(ddc_length_comps$sizecomp_cpue_stn, hauls_survey) %>%  #, by = c("vessel", "year", "haul"))
    select(cruisejoin, vessel, year, haul, start_latitude, start_longitude, stratum, stationid, sex, length, cpue_length_num_ha)
  
  # Age-length key - function name is switched!
  ddc_age_comps <- get_agecomps(pollock_specimen, ddc_length_comps$pollock_length_comp)
  
  # Age comps - function name is switched!
  ddc_al_key <- get_al_key(pollock_specimen, ddc_length_comps$pollock_length_comp)
  
  # DDC numbers-at-age by year and station ------------------------------------
  # This step takes a long time!
  ddc_alk_all <- ddc_age_comps_f(ddc_table,
                                 tables$age_comp_full_key,  # CIA: Stan uses UNCORRECTED AAL key
                                 tables$cpue_length_table)  # CIA: Stan uses UNCORRECTED CPUE number at length
  
  ddc_alk <- ddc_alk_all$ddc_age
  
  return(list(ddc_cpue = ddc_cpue, 
              ddc_pop_ests = ddc_pop_ests,
              ddc_pop_ests_EBS = ddc_pop_ests_EBS,
              ddc_pop_ests_NBS = ddc_pop_ests_NBS,
              ddc_cpue_length_table = ddc_cpue_length_table,
              ddc_length_comps = ddc_length_comps,
              ddc_cpue_length_table = ddc_cpue_length_table,
              ddc_al_key = ddc_al_key,
              ddc_alk_all = ddc_alk_all,
              ddc_alk = ddc_alk))
}

ddc <- ddc_conversion()
ddc_alk <- ddc$ddc_alk


# Apply Stan's bootstrapping method, if using design-based method -------------
bootstrapping <- function() {
  if(data_type == "mb") {
    print("Moving on! This function is only for the design-based data.")
  }
  
  if(data_type == "db") {
    bootstrap_stan(cpue_length_table = tables$cpue_length_table, cpue_info = tables$cpue_info)
  }
}

db_bootstrap <- bootstrapping()


# Check and save model-based results (for VAST) -------------------------------
# Repeated file path pieces
output <- here("output", dir_thisyr)
file_end <- paste0("_", current_year, ".csv")

if(data_type == 'mb') {
  VAST_files <- make_VAST_input(hauls = hauls_survey,
                                spec = pollock_specimen,
                                ddc_index = ddc_table,
                                ddc_age = ddc_alk,
                                slope_survey = tables$slope_survey)
  VAST_ddc_alk <- VAST_files$VAST_ddc_alk %>% 
    filter(Age >= 1) # check
  
  if(VAST_files$encounter_rate == FALSE) {print("no years with 100% encounter rate")
  }else{print("there are years with 100% encounter rate")}
  
  print(paste0("check table: sample sizes for each age should be the same"))
  print(table(VAST_ddc_alk$Age)) #check that sample size is same for each!
  
  ## Write out tables for VAST ------------------------------------------------
  write_csv(VAST_files$VAST_ddc_table,  # Biomass for EBS + NBS together - dd correction 
            here(output, paste0("VAST_ddc_all", file_end)))
  write_csv(VAST_files$VAST_ddc_table_EBS,  # Biomass for just EBS - dd correction
            here(output, paste0("VAST_ddc_EBSonly", file_end)))
  write_csv(VAST_files$VAST_ddc_table_NBS,  # Biomass for just NBS - dd correction
            here(output, paste0("VAST_ddc_NBSonly", file_end)))
  write_csv(VAST_ddc_alk,  # Age comps for EBS + NBS together - dd correction
            here(output, paste0("VAST_ddc_alk", file_end)))
}


# Save design-based results ---------------------------------------------------
# CIA: note to separate out EBS and NBS: NBS_subarea for stratum values; I set subarea of NBS = 0; all other values are EBS
if(data_type == 'db') {
  # Tables Jim Ianelli needs (all EBS + NBS):
  write_csv(tables$cpue_info$avg_cpue_yr_strat,  # CPUE (number of fish) - uncorrected 
            here(output, paste0("CPUE_uncorrected", file_end)))
  write_csv(tables$cpue_info$p_cpue,  # add CPUE by station  
            here(output, paste0("CPUE_bystn_uncorrected", file_end)))
  
  write_csv(ddc$ddc_cpue,  # CPUE (number of fish) - dd corrected  
            here(output, paste0("CPUE_densdep_corrected", file_end)))
  write_csv(ddc_table,  # add CPUE by station 
            here(output, paste0("CPUE_densdep_bystn_corrected", file_end)))
  
  write_csv(tables$cpue_length_table,  # CPUE (number of fish, by length) - uncorrected 
            here(output, paste0("length_cpue_numfish_uncorrected", file_end)))
  write_csv(tables$pop_info$pollock_biomass_MT_ha,  # Biomass (thousands of tons) - uncorrected 
            here(output, paste0("biomass_uncorrected", file_end)))
  write_csv(ddc$ddc_pop_ests$pollock_biomass_MT_ha,  # Biomass (thousands of tons) - dd corrected
            here(output, paste0("biomass_densdep_corrected", file_end)))
  write_csv(ddc$ddc_pop_ests_EBS,  # Biomass (thousands of tons) just EBS - dd corrected
            here(output, paste0("biomass_densdep_corrected_EBSonly", file_end)))
  write_csv(ddc$ddc_pop_ests_NBS,  # Biomass (thousands of tons) just NBS - dd corrected 
            here(output, paste0("biomass_densdep_corrected_NBSonly", file_end)))
  
  # Variance-covariance tables from Stan's method
  write_csv(as.data.frame(db_bootstrap$bioms), 
            here(output, paste0("bootstrap_biomass_tons_stan_method", current_year,"_", data_type, ".csv")))
  write_csv(as.data.frame(db_bootstrap$var_covar_stan), 
            here(output, paste0("var_cov_matrix_tons_stan_method", current_year,"_", data_type, ".csv")))
  
  write_csv(ddc$ddc_cpue_length_table,  # Length compositions - dd corrected 
            here(output, paste0("length_comps_densdep_corrected", file_end)))
  
  # Calculate metrics of hauls, number of lengths, etc. per year
  annual_metrics <- get_annual_metrics(hauls_ebs = tables$hauls_survey_ebs, 
                                       hauls_nbs = tables$hauls_survey_nbs, 
                                       length = tables$pollock_raw_length_ebs,
                                       length_nbs = tables$pollock_raw_length_nbs,
                                       spec = tables$pollock_specimen_ebs,
                                       spec_nbs = tables$pollock_specimen_nbs,
                                       this_yr = current_year, 
                                       nbs_subarea = NBS_subarea)
  
  write_lines(annual_metrics$num_hauls_tot,  # Number of stations (by year)
              here(output, paste0("info_hauls_peryr", file_end)))
  write_lines(annual_metrics$fish_measured,  # Number of fish measured (by year)
              here(output, paste0("info_n_fish_measured", file_end)))
  write_lines(annual_metrics$fish_aged,  # Number of fish aged (by year)
              here(output, paste0("info_n_fish_aged", file_end)))
  
  # Extra stuff
  write_csv(ddc_alk,  # Full DDC ALK
            here(output, paste0("age_length_key_full_densdep_corrected", file_end)))
  write_csv(tables$age_comp_full_key, 
            here(output, paste0("age_length_key_full_uncorrected_", file_end)))
  write_csv(ddc$ddc_al_key,  # DDC ALK summary
            here(output, paste0("age_length_key_SUMMARY_densdep_corrected", file_end)))
  
}


# Table requests by others ------------------------------------------------

# ddc age comps by sex
#Only run this when starting from scratch to create directory
dir.create(here("output",'other_requests'), showWarnings = FALSE)
write_csv(ddc$ddc_alk_all$ddc_age_sex, 
          here("output", "other_requests", 
               paste0("age_length_key_full_with_sex_densdep_corrected", current_year, ".csv")))

# length comps with haul joins
# ddc_cpue_length_table_Ivonne <- full_join(ddc_length_comps$sizecomp_cpue_stn, hauls_survey) %>%  #, by = c("vessel", "year", "haul"))
#   select(hauljoin, vessel, year, haul, start_latitude, start_longitude, stratum, stationid, sex, length, cpue_length_num_ha)
# 
# ddc_table_clean <- ddc_table %>% 
#   dplyr::select(-bs_stn_sum, -stn_mean_corr, -stn_sd_corr, -eff,
#                 -distance_fished, -net_width, -gear_temperature, -weight,
#                 -number_fish, -area_fished_ha) %>% 
#   dplyr::rename(uncorrected_cpue_kg_ha = cpue_kg_ha,
#                 uncorrected_cpue_num_ha = cpue_num_ha )
# write_csv(ddc_table_clean, here("output", "other requests", paste0("index_year_stn_FEAST", current_year, ".csv")))

# Biomass-by-length by subarea (for Ivonne) - SNW 2024-04-02
ddc_biom_by_length <- left_join(ddc$ddc_length_comps$pollock_length_comp_biomass, ddc$ddc_length_comps$stratum_use) %>%
  filter(!is.na(cpue_prop)) %>%
  mutate(biomass_kg = biomass_kg_ha * cpue_prop) %>%
  select("year", "stratum", "sex", "length", "biomass_kg")
ddc_biom_by_length <- ddc_biom_by_length[, -1]
write_csv(ddc_biom_by_length, 
          here("output", "other_requests", 
               paste0("biomass_by_length_stratum", current_year, ".csv")))

# total numbers by year, station, and 1 cm length bin
ddc_cpue_length_table <- ddc$ddc_cpue_length_table %<>%
  rename(ddc_cpue_length_num_ha = cpue_length_num_ha)
write_csv(ddc_cpue_length_table, here("output", "other_requests", paste0("length_numbers_year_stn_FEAST", current_year, ".csv")))
