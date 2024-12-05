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

# library(magrittr)
# library(lubridate)
# library(modelr)
# library(mgcv)
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
# UPDATE this section each year with the current cruise and vessels
current_year <- year(Sys.Date())
# current_year <- 2023  # choose a different year when debugging
cruise <- paste0(current_year, "01", ",", current_year, "02")
vessel_code <- paste0(162, "," , 134) #list each vessel, separated by commas
vessel_nums <- c(162, 134)

# get cruise id from here
query_command <- paste0(" select * from race_data.v_cruises where cruise in (", cruise,");")
cruise_info <- sqlQuery(channel, query_command) %>% 
  as_tibble() %>% 
  clean_names()

cruise_id_nums <- cruise_info %>% dplyr::filter(vessel_id %in% vessel_nums) 

cruise_id <- paste0(cruise_id_nums$cruise_id[1], "," , cruise_id_nums$cruise_id[2],
                    "," , cruise_id_nums$cruise_id[3], "," , cruise_id_nums$cruise_id[4])

# NBS subarea stratum-- this shouldn't change too much, but is a fixed input
NBS_subarea <- c(81, 70, 71, 99) # NBS stratum numbers; added 99 to indicate 2018 NBS emergency survey; diff survey methods

# Strata metadata year; 2022 is the latest update (use for current assessments)
strat_meta_year <- 2022

data_type <- "mb"

# Set up folder 
dir_thisyr <- paste0(current_year,"_", data_type, "_data_", strat_meta_year, "_strata")
dir.create(here("output",dir_thisyr))


# data --------------------------------------------------------------------
process_data <- function(first_run = TRUE, estimate_ages = FALSE, save_data = TRUE) {
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

  hauls_survey_ebs <- get_hauls$good_hauls_DB_EBS
  hauls_survey_nbs <- get_hauls$good_hauls_DB_NBS
  hauls_survey_bad_ebs <- get_hauls$bad_hauls_DB_ebs
  hauls_survey_bad_nbs <- get_hauls$bad_hauls_DB_nbs
  all_hauljoins <- c(hauls_survey_ebs$hauljoin, hauls_survey_bad_ebs$hauljoin)
  all_hauljoins_nbs <- c(hauls_survey_nbs$hauljoin, hauls_survey_bad_nbs$hauljoin)
  valid_hauljoins_nbs <- hauls_survey_nbs$hauljoin
  hauls_survey <- hauls_survey_ebs %>% bind_rows(hauls_survey_nbs)
  hauls_survey_bad <- hauls_survey_bad_ebs %>% bind_rows(hauls_survey_bad_nbs)
  
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
  if(data_type == "db")
  {
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
  
  ## Return tables for next step ----------------------------------------------
  return(list(hauls_survey = hauls_survey, 
              pollock_specimen = pollock_specimen,
              pollock_catch = pollock_catch,
              pollock_length = pollock_length,
              all_strata = all_strata,
              ddc_table = ddc_table))
}

tables <- process_data()






# ddc-converted tables ----------------------------------------------------



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


# population:
ddc_pop_ests <- population(ebs_strata = all_strata, avg_cpue = ddc_cpue)
ddc_pop_ests_EBS <- ddc_pop_ests$pollock_biomass_MT_ha %>% dplyr::filter(subarea != 0)
ddc_pop_ests_NBS <- ddc_pop_ests$pollock_biomass_MT_ha %>% dplyr::filter(subarea == 0)

# length comps:
ddc_length_comps <- length_comp_f(pollock_length, hauls_survey, pollock_catch, 
              ddc_pop_ests$ebs_strata, 
              ddc_pop_ests$pollock_biomass_kg_ha, 
              ddc_pop_ests$pollock_biomass_MT_ha)

ddc_cpue_length_table <- full_join(ddc_length_comps$sizecomp_cpue_stn, hauls_survey) %>%  #, by = c("vessel", "year", "haul"))
  select(cruisejoin, vessel, year, haul, start_latitude, start_longitude, stratum, stationid, sex, length, cpue_length_num_ha)

# age-length key:
ddc_age_comps <- get_agecomps(pollock_specimen, ddc_length_comps$pollock_length_comp)

# age comps:
ddc_al_key <- get_al_key(pollock_specimen, ddc_length_comps$pollock_length_comp)


# ddc numbers at age by year and station ----------------------------------
# SNW: this is where it starts taking a while

ddc_alk_all <- ddc_age_comps_f(ddc_table,
                              age_comp_full_key,  #CIA: Stan uses UNCORRECTED aal key
                              cpue_length_table) # CIA: Stan uses UNCORRECTED cpue number at length

ddc_alk <- ddc_alk_all$ddc_age

# convert new 4 -----------------------------------------------------------
# this is the conversion of stan's 'new4.R' which sets up the bootstrapping procedure to generate a v-cov matrix later
# CORRECTION: Caitlin's method did not produce the same amount of variance in the v-cov matrix, 
#   so that code is commented out and stan's code is implemented directly

##### caitlin's section ####

# # resample the data for each year and strata
# resample <- resamp_data(corr_cpue = ddc_table)
# 
# # check:
# # mean(resample$ddc_cpue_kg_ha)
# # mean(ddc_table$ddc_cpue_kg_ha)
# # mean(resample$ddc_cpue_num_ha, na.rm = T)
# # mean(ddc_table$ddc_cpue_num_ha, na.rm = T)
# 
# # cpue per station for one resample
# resample_pop <- resamp_pop(resamp = resample, strata = strata_metadata) %>% 
#   dplyr::filter(!is.na(year))
# 
# # CIA: YOU ARE HERE (last part of new4.R)
# param_sample <- sample(1000,1000,replace=T)
# 
# biomass <- matrix(nrow = length(param_sample), ncol = dim(resample_pop)[1])
# population <- matrix(nrow = length(param_sample), ncol = dim(resample_pop)[1])
# 
# 
# v_cov_units <- readline(prompt = "Which units do you want to use in the v-cov matrix (Enter: k for kg or t for th t): ")
# 
# k 
# 
# if(v_cov_units == 'k') 
# {
#   biomass <- foreach(i = 1:1000, .combine=rbind) %dopar% {
#     library(magrittr)
#     library(tidyverse)
#     
#     param_sample <- sample(1000,1,replace=T)
#     
#     # estimate ddc with random sample from parameter posterior distribution
#     backscatter_cpue_corr <- ddc_cpue_bs(equiv_bs_length, param_ests[param_sample,])
#     
#     # turn backscatter into cpue
#     ddc_table_vc <- ddc_fxn(backscatter_cpue_corr, cpue_info$p_cpue) #backscatter corr from EACH SAMPLE
#     
#     # randomly resample data from all years and all stations
#     resample <- resamp_data(ddc_table_vc) #FOR EACH SAMPLE
#     
#     # get total population est from each resample for each year (cpue*area)
#     resample_pop <- resamp_pop(resample, strata_metadata) %>%  #FOR EACH SAMPLE
#       dplyr::filter(!is.na(year))
#     
#     # save resampled biomass
#     resample_pop$total_ddc_kg #EACH ROW IS RESAMPLE, EACH COL IS YEARS
#   }
#   
#   write_csv(as.data.frame(biomass), here("output",dir_thisyr,paste0("bootstrap_biomass_kg_", data_type, ".csv")))
# } else if(v_cov_units == 't') {
#   biomass <- foreach(i = 1:1000, .combine=rbind) %dopar% {
#     library(magrittr)
#     library(tidyverse)
#     
#     param_sample <- sample(1000,1,replace=T)
#     
#     # estimate ddc with random sample from parameter posterior distribution
#     backscatter_cpue_corr <- ddc_cpue_bs(equiv_bs_length, param_ests[param_sample,])
#     
#     # turn backscatter into cpue
#     ddc_table_vc <- ddc_fxn(backscatter_cpue_corr, cpue_info$p_cpue) #backscatter corr from EACH SAMPLE
#     
#     # randomly resample data from all years and all stations
#     resample <- resamp_data(ddc_table_vc) #FOR EACH SAMPLE
#     
#     # get total population est from each resample for each year (cpue*area)
#     resample_pop <- resamp_pop(resample, strata_metadata) %>%  #FOR EACH SAMPLE
#       dplyr::filter(!is.na(year))
#     
#     # save resampled biomass
#     resample_pop$total_ddc_th_t #EACH ROW IS RESAMPLE, EACH COL IS YEARS
#   }
#   
#   write_csv(as.data.frame(biomass), here("output",dir_thisyr,paste0("bootstrap_biomass_th_t_", data_type, ".csv")))
# } else(print("Invalid selection: please select k for kg or t for th t"))

# v-cov check
# sqrt(diag(biomass))/mean(diag(biomass))

#### stan's bootstrapping method: ####
if(data_type == 'db')
{
corr_sa <-  function(aa,ab){
  
  #calculates corrected sa by tow (corrected backscatter by tow)
  
  by1=list(haul=aa$HAUL, subst=aa$SUBSTRATA, vessel=aa$VESSEL, year=aa$YEAR)
  sa_station = aggregate(aa$bt_sa, by1, sum)
  str(sa_station)
  #hist(sa_station$x, 2000, xlim=c(0,5000))
  
  ii=1
  sa_corr=sa_station$x[ii]/(ab$asp+exp(-(ab$int+ab$slope*sa_station$x[ii])))
  corr_mean=mean(sa_corr)
  corr_sd = sd(sa_corr)
  
  for(ii in 2:length(sa_station$x)){
    sa_corr=sa_station$x[ii]/(ab$asp+exp(-(ab$int+ab$slope*sa_station$x[ii])))
    corr_mean[ii]=mean(sa_corr)
    corr_sd[ii] = sd(sa_corr) #might need na.rm = T
  }
  
  #this data contains old and new estimates of sa (uncorrected and corrected backscatter estimates returned)
  sa_sum=cbind(sa_station,corr_mean,corr_sd)
  return(sa_sum)
}

corr_cpue <-  function(sa,bb){
  
  #calculates new cpue per tow (applies bs correction; eff = efficency; eff = ratio uncorr to corr, less than 1)
  
  eff=sa$x/sa$corr_mean
  sa=cbind(sa,eff)
  #str(sa)
  
  cc=merge(bb,sa,all=T)
  
  new_kgha=cc$CPUE_KGHA/cc$eff
  new_noha=cc$CPUE_NOHA/cc$eff
  
  pollock_cpue_new=cbind(cc[,c(1:11,16)],new_kgha,new_noha)
  #str(pollock_cpue_new)
  pollock_cpue_new$new_kgha[pollock_cpue_new$CPUE_KGHA==0]=0
  pollock_cpue_new$new_noha[pollock_cpue_new$CPUE_NOHA==0]=0
  pollock_cpue_new$eff[pollock_cpue_new$CPUE_KGHA==0]=1
  #write.csv(pollock_cpue_new, "pollock_cpue_new.csv")
  
  na.eff.cpue = data.frame(xx=pollock_cpue_new[which(is.na(pollock_cpue_new$eff)),]$CPUE_KGHA)
  na.eff.cpue.idx=which(is.na(pollock_cpue_new$eff))
  
  xx=pollock_cpue_new$CPUE_KGHA
  yy=pollock_cpue_new$eff
  #plot(yy~xx)
  gam1=gam(yy~s(xx))
  #summary(gam1)
  #plot(gam1)
  
  na.eff=predict(gam1,na.eff.cpue)
  pollock_cpue_new$eff[na.eff.cpue.idx]=na.eff
  pollock_cpue_new$new_kgha[na.eff.cpue.idx]=pollock_cpue_new$CPUE_KGHA[na.eff.cpue.idx]/na.eff
  pollock_cpue_new$new_noha[na.eff.cpue.idx]=pollock_cpue_new$CPUE_NOHA[na.eff.cpue.idx]/na.eff
  #write.csv(pollock_cpue_new, "pollock_cpue_new.csv")
  #str(pollock_cpue_new)
  
  
  
  pollock_cpue_new2=pollock_cpue_new[,c(1:8,13:14,11)]
  #str(pollock_cpue_new2)
  colnames(pollock_cpue_new2)[c(1:3,9,10)]=c("VESSEL","YEAR","HAUL","CPUE_KGHA","CPUE_NOHA")
  #str(pollock_cpue_new2)
  
  return(pollock_cpue_new2) #returns the DDC cpue
  
}

ncpue_sample <-  function(ncpue){
  first.y=min(ncpue$YEAR) #for each year
  last.y=max(ncpue$YEAR)
  
  ii=first.y
  ncpue.y=ncpue[ncpue$YEAR==ii,] #select a specific year
  subs=unique(ncpue.y$SUBSTRATA) #for each strata; get strata from that year
  jj=subs[1]
  ncpue.y.s=ncpue.y[ncpue.y$SUBSTRATA==jj,] #select specific year and strata
  nn=length(ncpue.y.s[,1]) #number samples wanted
  samp=sample.int(nn, size = nn, replace=T) # randomly pick order, possibly with duplicates
  ncpue.y.s.n=ncpue.y.s[samp,] #arrange in random order
  for(jj in subs[-1]){
    ncpue.y.s=ncpue.y[ncpue.y$SUBSTRATA==jj,] #
    nn=length(ncpue.y.s[,1])
    samp=sample.int(nn, size = nn, replace=T)
    ncpue.y.s.n=rbind(ncpue.y.s.n,ncpue.y.s[samp,])
    # estimate of one realization of estimating cpue; with all substrata, 1st year
  }
  ncpue.y.n=ncpue.y.s.n
  
  for(ii in (first.y+1):last.y){ #this is for all remaining years (with all substrata)
    ncpue.y=ncpue[ncpue$YEAR==ii,]
    subs=unique(ncpue.y$SUBSTRATA)
    jj=subs[1]
    ncpue.y.s=ncpue.y[ncpue.y$SUBSTRATA==jj,]
    nn=length(ncpue.y.s[,1])
    samp=sample.int(nn, size = nn, replace=T)
    ncpue.y.s.n=ncpue.y.s[samp,]
    for(jj in subs[-1]){
      ncpue.y.s=ncpue.y[ncpue.y$SUBSTRATA==jj,]
      nn=length(ncpue.y.s[,1])
      samp=sample.int(nn, size = nn, replace=T)
      ncpue.y.s.n=rbind(ncpue.y.s.n,ncpue.y.s[samp,])
    }
    ncpue.y.n=rbind(ncpue.y.n,ncpue.y.s.n)
    
  }
  return(ncpue.y.n) #one realization of cpue est for all years and all substrata
}

corr_population <-  function(aa,bb){
  
  by1=list(aa$YEAR,aa$SUBSTRATA)
  cc=aggregate(cbind(aa$CPUE_KGHA,aa$CPUE_NOHA),by1,mean,na.rm=T) #get mean by year and strata
  colnames(cc)=c("year","substrata","m_kgha","m_noha")
  colnames(bb)=c("substrata","area")
  cc1=merge(cc,bb,all=T) #add area back
  # cc1=cc1[cc1$substrata!=70,] #filter out NBS
  # cc1=cc1[cc1$substrata!=81,]
  pop_by_strata = cc1[,3:4]*cc1[,5] #get number by multiplying cpue by area fished (cpue = num/ha; area = km^2; need to *100 to get km^2 -> ha-- done 3 lines below)
  by2=list(cc1$year)
  total_pop=aggregate(pop_by_strata,by2,sum) #sum by year
  total_pop1=cbind(year=total_pop[,1],tons=total_pop[,2]/10,mlns=total_pop[,3]/10000) # years, convert to kg to tons, convernt num to millions
  return(total_pop1)
}

ab <- param_ests

bb <- cpue_info$p_cpue %>%  #cpue
  rename(CRUISE = cruise, SUBSTRATA = stratum, 
         LATITUDE = start_latitude, LONGITUDE = start_longitude,
         CPUE_KGHA = cpue_kg_ha, CPUE_NOHA = cpue_num_ha, AREA_FISHED = area_fished_ha) %>% 
  mutate(SPECIES_CODE = 21740) %>% 
  dplyr::select("vessel", "CRUISE", "year", "haul", "SUBSTRATA", "LATITUDE", "LONGITUDE", 
                "SPECIES_CODE", "CPUE_KGHA", "CPUE_NOHA", "AREA_FISHED") #%>% 
  # dplyr::filter(!SUBSTRATA %in% NBS_subarea)
  
# bb_sub <- bb %>% dplyr::select("vessel", "CRUISE", "year", "haul", "SUBSTRATA", "LATITUDE", "LONGITUDE") %>% 
#   clean_names(case = 'all_caps')

aa <- equiv_bs_length %>%  #cpue length
  clean_names(case = 'all_caps') %>% 
  rename(cpue_mile = CPUE_MILE_SQ, ts = TARGET_STRENGTH, bt_sa = BACKSCATTER, 
         CRUISE = CRUISEJOIN, SUBSTRATA = STRATUM, CPUE_LENGTH = CPUE_LENGTH_NUM_HA) %>%
  # left_join(bb_sub) %>% 
  mutate( LATITUDE = NA, LONGITUDE = NA, 
    SPECIES_CODE = 21740) %>% 
  dplyr::select("VESSEL", "YEAR", "CRUISE", "HAUL", "LATITUDE", "LONGITUDE", "SUBSTRATA", 
                "STATIONID", "SPECIES_CODE", "SEX", "LENGTH", "CPUE_LENGTH", 
                "cpue_mile",  "ts", "bt_sa") #%>% 
  # dplyr::filter(!SUBSTRATA %in% NBS_subarea)

gg <- strata_metadata %>%  #strata
  rename(area_ha = area) %>% 
  dplyr::select( "stratum", "area_ha") #%>% 
  # dplyr::filter(!stratum %in% NBS_subarea)
  
ii=sample(1000,1) #randome number btw 1 and 1000
sa_sum = corr_sa(aa,ab[ii,]) #estimates corr sa based on params from abc_exp
ncpue = corr_cpue(sa_sum,bb) #estimates cpue based on param
ncpue.s = ncpue_sample(ncpue) #gets ncpue sample
estim=corr_population(ncpue.s,gg) #get corr pop based on random param selected
pops=estim[,3] #pop vector
bioms=estim[,2] #biomass vector

for(ii in sample(1000,999,replace=T)){ #run 1000 times (999 more times), randomly sampling
  sa_sum = corr_sa(aa,ab[ii,])
  ncpue = corr_cpue(sa_sum,bb)
  ncpue.s = ncpue_sample(ncpue)
  estim=corr_population(ncpue.s,gg)
  popsii=estim[,3]
  biomsii=estim[,2]
  pops=rbind(pops,popsii) #rows = resamples, cols = yrs
  bioms=rbind(bioms,biomsii)
  print(length(pops[,1]))
}
colnames(bioms) = unique(aa$YEAR)
bioms
var_covar_stan = cov(bioms)
colnames(var_covar_stan) = unique(aa$YEAR)
rownames(var_covar_stan) = unique(aa$YEAR)

# check:
biomass_mean <- bioms %>% as_tibble() %>% summarise_if(is.numeric, mean) 
colnames(biomass_mean) <- unique(ddc_table$year)
cv <- sqrt(diag(var_covar_stan))/biomass_mean
range(cv) #range should be 5-20%
if(max(range(cv)) > 0.2){print("YOUR CV IS > 20% PLEASE CHECK")}
}

# convert new 5 -----------------------------------------------------------

# CORRECTION: this section is commented out-- v-cov matrix generated above in stan's code version

# get variance-covariance matrix from resampled data
# var_covar_pop_num = cov(population)

# var_covar_biomass = cov(biomass)
# 
# # then rename rows and columns
# colnames(var_covar_biomass) <- resample_pop$year
# rownames(var_covar_biomass) <- resample_pop$year
# 
# # check:
# biomass_mean <- biomass %>% as_tibble() %>% summarise_if(is.numeric, mean) 
# colnames(biomass_mean) <- unique(ddc_table$year)
# cv <- sqrt(diag(var_covar_biomass))/biomass_mean
# range(cv) #range should be 5-20%
# if(max(range(cv)) > 0.2){print("YOUR CV IS > 20% PLEASE RE-CHECK")}

# # ianelli check:
# c1 <- var_covar_biomass
# d1 <- sqrt(diag(c1))
# d2 <- as.vector(d1/1000)
# c2 <- cov2cor(c1)
# c3 <- sweep(sweep(c2, 1L, d2), 2L, d2)
# all.equal(c3,c1)

# measured ----------------------------------------------------------------

if(data_type == 'db')
{
  annual_metrics <- get_annual_metrics(hauls_ebs = hauls_survey_ebs, 
                                       hauls_nbs = hauls_survey_nbs, 
                                       length = pollock_raw_length_ebs,
                                       length_nbs = pollock_raw_length_nbs,
                                       spec = pollock_specimen_ebs,
                                       spec_nbs = pollock_specimen_nbs,
                                       this_yr = current_year, 
                                       nbs_subarea = NBS_subarea)
  num_hauls_tot <- annual_metrics$num_hauls_tot
  fish_measured <- annual_metrics$fish_measured
  fish_aged <- annual_metrics$fish_aged
}


# VAST prep section -------------------------------------------------------
if(data_type == 'mb')
{
  VAST_files <- make_VAST_input(hauls = hauls_survey,
                                spec = pollock_specimen,
                                ddc_index = ddc_table,
                                ddc_age = ddc_alk)
  VAST_ddc_table <- VAST_files$VAST_ddc_table
  VAST_ddc_table_EBS <- VAST_files$VAST_ddc_table_EBS
  VAST_ddc_table_NBS <- VAST_files$VAST_ddc_table_NBS
  VAST_ddc_alk <- VAST_files$VAST_ddc_alk %>% 
    dplyr::filter(Age >= 1) #check
  
  if(VAST_files$encounter_rate == FALSE) {print("no years with 100% encounter rate")
  }else{print("there are years with 100% encounter rate")}
  
  print(paste0("check table: sample sizes for each age should be the same"))
  print(table(VAST_ddc_alk$Age)) #check that sample size is same for each!
}
# final output ------------------------------------------------------------

# CIA: note to separate out EBS and NBS: NBS_subarea for stratum values; I set subarea of NBS = 0; all other values are EBS
if(data_type == 'db')
{
  #   Tables Jim Ianelli needs (all EBS + NBS):
  #     - CPUE (number of fish) - uncorrected           cpue_info$avg_cpue_yr_strat
  #     - CPUE (kg/hectare)     - uncorrected           cpue_info$avg_cpue_yr_strat
  write_csv(cpue_info$avg_cpue_yr_strat, here("output", dir_thisyr,paste0(
    "CPUE_uncorrected_", current_year, ".csv"
  )))
  #     --> add CPUE by station                         cpue_info$p_cpue
  write_csv(cpue_info$p_cpue, here(
    "output", dir_thisyr,
    paste0("CPUE_bystn_uncorrected_", current_year, ".csv")
  ))
  
  #     - CPUE (number of fish) - dd corrected          ddc_cpue
  #     - CPUE (kg/hectare)     - dd corrected          ddc_cpue
  write_csv(ddc_cpue, here(
    "output", dir_thisyr,
    paste0("CPUE_densdep_corrected_", current_year, ".csv")
  ))
  #     --> add CPUE by station                         ddc_table
  write_csv(ddc_table, here(
    "output", dir_thisyr,
    paste0("CPUE_densdep_bystn_corrected_", current_year, ".csv")
  ))
  
  
  #     - CPUE (number of fish, by length) - uncorrected  length_comps$sizecomp_cpue_stn
  write_csv(cpue_length_table, here(
    "output", dir_thisyr,
    paste0("length_cpue_numfish_uncorrected_", current_year, ".csv")
  ))
  #     - CPUE (kg/hectare, by length)     - uncorrected  length_comps$sizecomp_cpue_stn; cpue_length_table CIA: in number/ha
  #
  #     - Biomass (thousands of tons) - uncorrected     pop_info$pollock_biomass_MT_ha
  write_csv(pop_info$pollock_biomass_MT_ha, here(
    "output", dir_thisyr,
    paste0("biomass_uncorrected_", current_year, ".csv")
  ))
  #     - Biomass (thousands of tons) - dd corrected    ddc_pop_ests$pollock_biomass_MT_ha
  write_csv(ddc_pop_ests$pollock_biomass_MT_ha,
            here(
              "output", dir_thisyr,
              paste0("biomass_densdep_corrected_", current_year, ".csv")
            ))
  #     - Biomass (thousands of tons) just EBS - dd corrected   ddc_pop_ests_EBS
  write_csv(ddc_pop_ests_EBS, here(
    "output", dir_thisyr,
    paste0("biomass_densdep_corrected_EBSonly_", current_year, ".csv")
  ))
  #     - Biomass (thousands of tons) just NBS - dd corrected   ddc_pop_ests_NBS
  write_csv(ddc_pop_ests_NBS, here(
    "output", dir_thisyr,
    paste0("biomass_densdep_corrected_NBSonly_", current_year, ".csv")
  ))
  #
  #     - Co-variance matrix (biomass; thousands of tons) - dd corrected    var_covar_biomass
  #         -> additional table converted to th t
  # # Caitlin's tables are commented out; Stan's tables are the new default output
  # if (v_cov_units == 'k')
  # {
  #   write_csv(as.data.frame(var_covar_biomass),
  #             here(
  #               "output", dir_thisyr,
  #               paste0("var_cov_matrix_kg_", current_year,"_", data_type, ".csv")
  #             ))
  # } else if (v_cov_units == 't')
  # {
  #   write_csv(as.data.frame(var_covar_biomass),
  #             here(
  #               "output", dir_thisyr,
  #               paste0("var_cov_matrix_th_t_", current_year, "_", data_type, ".csv")
  #             ))
  # }
  
  write_csv(as.data.frame(bioms), here(
    "output", dir_thisyr,
    paste0("bootstrap_biomass_tons_stan_method", current_year,"_", data_type, ".csv")
  ))
  write_csv(as.data.frame(var_covar_stan), here(
    "output", dir_thisyr,
    paste0("var_cov_matrix_tons_stan_method", current_year,"_", data_type, ".csv")
  ))
  
  #
  #     - Length compositions - dd corrected            ddc_cpue_length_table
  write_csv(ddc_cpue_length_table, here(
    "output", dir_thisyr,
    paste0("length_comps_densdep_corrected_", current_year, ".csv")
  ))
  
  #     OTHER INFO:
  #       - number of stations (by year)          num_hauls_tot
  write_lines(num_hauls_tot, here("output", dir_thisyr, paste0(
    "info_hauls_peryr_", current_year, ".txt"
  )))
  #       - number of fish measured (by year)     fish_measured
  write_lines(fish_measured, here(
    "output", dir_thisyr,
    paste0("info_n_fish_measured_", current_year, ".txt")
  ))
  #       - number of fish aged (by year)         fish_aged
  write_lines(fish_aged, here("output", dir_thisyr, paste0(
    "info_n_fish_aged_", current_year, ".txt"
  )))
  
  #       - station-specific temperature (bottom temp) + add to station cpue table
  
  #       - extra: full ddc alk
  write_csv(ddc_alk, here(
    "output", dir_thisyr,
    paste0(
      "age_length_key_full_densdep_corrected",
      current_year,
      ".csv"
    )
  ))
  write_csv(age_comp_full_key, here(
    "output", dir_thisyr,
    paste0("age_length_key_full_uncorrected_", current_year, ".csv")
  ))
  # age comps summary:
  write_csv(ddc_al_key, here(
    "output", dir_thisyr,
    paste0(
      "age_length_key_SUMMARY_densdep_corrected",
      current_year,
      ".csv"
    )
  ))
  
}

if(data_type == 'mb')
{
  #     PREP for VAST:
  #      - biomass for EBS + NBS together  - dd correction   VAST_ddc_table
  write_csv(VAST_ddc_table, here("output", dir_thisyr,paste0("VAST_ddc_all_", current_year, ".csv")))
  #      - biomass for just EBS            - dd correction
  write_csv(VAST_ddc_table_EBS, here("output", dir_thisyr,paste0(
    "VAST_ddc_EBSonly_", current_year, ".csv"
  )))
  #      - biomass for just NBS            - dd correction
  write_csv(VAST_ddc_table_NBS, here("output", dir_thisyr,paste0(
    "VAST_ddc_NBSonly_", current_year, ".csv"
  )))
  #      - age comps for EBS + NBS together- dd correction
  write_csv(VAST_ddc_alk, here("output", dir_thisyr,paste0("VAST_ddc_alk_", current_year, ".csv")))
}

# Table requests by others ------------------------------------------------

# ddc age comps by sex
#Only run this when starting from scratch to create directory
dir.create(here::here("output",'other_requests'), showWarnings = FALSE)
write_csv(ddc_alk_all$ddc_age_sex, here("output", "other_requests", paste0("age_length_key_full_with_sex_densdep_corrected", current_year, ".csv")))
# drive_upload(ddc_alk_all$ddc_age_sex, 
#              path = (as_id("1FXW5eaQh28ZnmY83mJHCWnWRmSJR9lGL")), 
#              name = paste0("age_length_key_full_with_sex_densdep_corrected", current_year, ".csv"))
# 
# ddc_alk_all$ddc_age_sex %>% 
#   drive_upload(path = (as_id("1FXW5eaQh28ZnmY83mJHCWnWRmSJR9lGL")),
#                name = paste0("age_length_key_full_with_sex_densdep_corrected", current_year, ".csv"))


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
ddc_biom_by_length <- left_join(ddc_length_comps$pollock_length_comp_biomass, ddc_length_comps$stratum_use) %>%
  filter(!is.na(cpue_prop)) %>%
  mutate(biomass_kg = biomass_kg_ha * cpue_prop) %>%
  select("year", "stratum", "sex", "length", "biomass_kg")
ddc_biom_by_length <- ddc_biom_by_length[, -1]
write_csv(ddc_biom_by_length, here("output", "other_requests", paste0("biomass_by_length_stratum", current_year, ".csv")))

# total numbers by year, station, and 1 cm length bin
ddc_cpue_length_table %<>%
  rename(ddc_cpue_length_num_ha = cpue_length_num_ha)
write_csv(ddc_cpue_length_table, here("output", "other_requests", paste0("length_numbers_year_stn_FEAST", current_year, ".csv")))
