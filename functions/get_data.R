# These functions get data needed for ddc conversion
# Author: Caitlin I. Allen Akselrud
# Contact: caitlin.allen_akselrud@noaa.gov
# Date created: 2022.06.02
# Date updated: 2022.06.02

slope_survey_d <- function()
{
  query_command <- paste0(" select * from SAFE.survey;")
  all_survey_info <- sqlQuery(channel, query_command) %>% 
    as_tibble() %>% 
    clean_names()
  
  unique(all_survey_info$survey)
  
  slope_survey <- all_survey_info %>% dplyr::filter(survey == "EBS_SLOPE") %>% distinct(cruisejoin)
  return(slope_survey)
}

haul_data_d <- function(data_selection = "db", nbs_subarea = c(81, 70, 71), slope_info)
{
  # DB: design-based: ebs and nbs separate; nbs = 2010+; no 2018 nbs
  # through 2021
  query_command <- paste0("select * from SAFE.pollock_historical_design_hauls")
  stan_good_hauls <- sqlQuery(channel, query_command) %>%
    as_tibble() %>%
    clean_names()
  
  # CIA: get this table into oracle and update this section
  
  query_command <- paste0("select a.CRUISEJOIN, a.HAULJOIN, a.REGION, a.VESSEL, a.CRUISE,
                            a.HAUL, a.HAUL_TYPE, a.PERFORMANCE, a.START_TIME, a.DURATION,
                            a.DISTANCE_FISHED, a.NET_WIDTH, a.NET_MEASURED, a.NET_HEIGHT,
                            a.STRATUM, a.START_LATITUDE, a.END_LATITUDE, a.START_LONGITUDE,
                            a.END_LONGITUDE, a.STATIONID, a.GEAR_DEPTH, a.BOTTOM_DEPTH,
                            a.BOTTOM_TYPE, a.SURFACE_TEMPERATURE, a.GEAR_TEMPERATURE,
                            a.WIRE_LENGTH, a.GEAR, a.ACCESSORIES, a.SUBSAMPLE, a.AUDITJOIN,
                            floor(a.cruise/100) year
                            from racebase.haul a
                            where a.PERFORMANCE >=0 and a.haul_type = 3 and a.region = 'BS'
                            order by a.cruise, a.vessel, a.haul;")
                            # remove restriction and correct for missing stratum in 2022
                            # and a.stratum is not null and a.stationid is not null

  hauls_survey_orig <- sqlQuery(channel, query_command) %>%
    as_tibble() %>%
    clean_names()
  
  hist_hauls <- hauls_survey_orig %>% #use historical hauls based on Stan's selection of 'good hauls'
    dplyr::filter(hauljoin %in% stan_good_hauls$hauljoin)
  new_hauls <- hauls_survey_orig %>% 
    dplyr::filter(year > 2021)
  good_hauls_DB <- bind_rows(hist_hauls, new_hauls)
  #then separate EBS and NBS with strata
  good_hauls_DB_EBS <- good_hauls_DB %>% filter(!stratum %in% nbs_subarea)
  good_hauls_DB_NBS <- hauls_survey_orig %>% filter(stratum %in% nbs_subarea,
                                                year >= 2010)
  
  ##############################################################
  
  # MB: model-based: ebs and nbs together; NBS include pre 2010 and 2010+; NBS include 2018
  if(data_selection == "mb") {nbs_subarea = c(nbs_subarea, 99)} #99 is the 2018 special tow
  
  # valid_hauls <- read_rds(here("data","valid_hauls.rds"))
  # all_hauljoins <- c(valid_hauls$EBS_hauljoins, valid_hauls$NBS_hauljoins)
  
  query_command <- paste0("select * from SAFE.ebs_vast_hauls
                        where region = 'BS' 
                        and haul_type = 3 and performance>=0")
  valid_hauls_ebs <- sqlQuery(channel, query_command) %>%
    as_tibble() %>%
    clean_names()
  
  query_command <- paste0("select * from SAFE.nbs_vast_hauls
                        where region = 'BS' 
                        and haul_type = 3 and performance>=0")
  valid_hauls_nbs <- sqlQuery(channel, query_command) %>%
    as_tibble() %>%
    clean_names()
  
  all_hauljoins <- bind_rows(valid_hauls_ebs, valid_hauls_nbs) 
  #valid hauls should be updated annually in Oracle database
  
  # hist_hauls <- hauls_survey_orig %>% 
  #   dplyr::filter(hauljoin %in% all_hauljoins)
  # new_hauls <- hauls_survey_orig %>% 
  #   dplyr::filter(year > 2021)
  # good_hauls_MB <- bind_rows(hist_hauls, new_hauls) %>% 
  #   dplyr::filter(!cruisejoin %in% slope_info$cruisejoin)
  
  good_hauls_MB <- hauls_survey_orig %>% 
    dplyr::filter(hauljoin %in% all_hauljoins$hauljoin) %>% #valid vast hauls
    dplyr::filter(!cruisejoin %in% slope_info$cruisejoin) #make sure no slope survey data are included
  
  
  query_command <- paste0("select * from racebase.haul
                        where region = 'BS' and cruise = 201801
                        and haul_type = 13 and performance>=0")
  hauls_survey_18 <- sqlQuery(channel, query_command) %>%
    as_tibble() %>%
    clean_names() %>%
    mutate(stratum = 99) %>%
    mutate(year = 2018)
  # hauls_survey_18$start_time <- as.character(hauls_survey_18$start_time)
  
  good_hauls_MB_all <- bind_rows(good_hauls_MB, hauls_survey_18)
  
  ##############################################################
  
  # All standard bad hauls:
  bad_hauls <- sqlQuery(channel,
                        paste0("select a.CRUISEJOIN, a.HAULJOIN, a.REGION, a.VESSEL, a.CRUISE,
                                      a.HAUL, a.HAUL_TYPE, a.PERFORMANCE, a.START_TIME, a.DURATION,
                                      a.DISTANCE_FISHED, a.NET_WIDTH, a.NET_MEASURED, a.NET_HEIGHT,
                                      a.STRATUM, a.START_LATITUDE, a.END_LATITUDE, a.START_LONGITUDE,
                                      a.END_LONGITUDE, a.STATIONID, a.GEAR_DEPTH, a.BOTTOM_DEPTH,
                                      a.BOTTOM_TYPE, a.SURFACE_TEMPERATURE, a.GEAR_TEMPERATURE,
                                      a.WIRE_LENGTH, a.GEAR, a.ACCESSORIES, a.SUBSAMPLE, a.AUDITJOIN,
                                      floor(a.cruise/100) year
                               from racebase.haul a
                               where a.PERFORMANCE < 0 and a.haul_type = 3 and a.region = 'BS'
                                      and a.stratum is not null and a.stationid is not null
                               order by a.cruise, a.vessel, a.haul;")) %>% 
    as_tibble() %>% 
    clean_names()
  
  bad_hauls_DB_ebs <- bad_hauls %>% 
    dplyr::filter(!stratum %in% nbs_subarea)
  bad_hauls_DB_nbs <- bad_hauls %>% 
    dplyr::filter(stratum %in% nbs_subarea,
                  year >= 2010)
  
  # + 2018 NBS bad hauls
  bad18_hauls <- sqlQuery(channel,
                          paste0("select a.CRUISEJOIN, a.HAULJOIN, a.REGION, a.VESSEL, a.CRUISE,
                                      a.HAUL, a.HAUL_TYPE, a.PERFORMANCE, a.START_TIME, a.DURATION,
                                      a.DISTANCE_FISHED, a.NET_WIDTH, a.NET_MEASURED, a.NET_HEIGHT,
                                      a.STRATUM, a.START_LATITUDE, a.END_LATITUDE, a.START_LONGITUDE,
                                      a.END_LONGITUDE, a.STATIONID, a.GEAR_DEPTH, a.BOTTOM_DEPTH,
                                      a.BOTTOM_TYPE, a.SURFACE_TEMPERATURE, a.GEAR_TEMPERATURE,
                                      a.WIRE_LENGTH, a.GEAR, a.ACCESSORIES, a.SUBSAMPLE, a.AUDITJOIN,
                                      floor(a.cruise/100) year
                               from racebase.haul a
                               where a.PERFORMANCE < 0 and a.haul_type = 13 and a.region = 'BS'
                                      and a.stratum is not null and a.stationid is not null
                               order by a.cruise, a.vessel, a.haul;")) %>% 
    as_tibble() %>% 
    clean_names() %>% 
    dplyr::filter(year == 2018)
  
  bad_hauls_MB <- bind_rows(bad_hauls, bad18_hauls)
  
  ##############################################################
  if(data_selection == "db")
  {
    return(list(good_hauls_DB_EBS = good_hauls_DB_EBS,
                good_hauls_DB_NBS = good_hauls_DB_NBS,
                bad_hauls_DB_ebs = bad_hauls_DB_ebs,
                bad_hauls_DB_nbs = bad_hauls_DB_nbs))
  }else if(data_selection == "mb")
  {
    return(list(good_hauls_MB = good_hauls_MB_all,
                bad_hauls = bad_hauls_MB))
  }
  
}

specimen_data_d <- function(hauls_survey_dat)
{
  hauls_survey_dat <- hauls_survey
  hauljoins <- hauls_survey_dat$hauljoin
  # SNW: remove "and a.age is not null" to create estimated ALK for current year
  query_command <- paste0("select a.CRUISEJOIN, a.HAULJOIN, a.REGION, a.VESSEL, a.CRUISE, a.HAUL,
                              a.SPECIMENID, a.BIOSTRATUM, a.SPECIES_CODE, round(a.LENGTH/10)*10 length, a.SEX, a.WEIGHT,
                              a.AGE, a.MATURITY, a.MATURITY_TABLE, a.GONAD_WT, a.AUDITJOIN
                            from racebase.specimen a
                            where species_code = 21740 and a.age is not null
                            order by cruise, vessel, haul;")
  
  pollock_specimen_orig <- sqlQuery(channel, query_command) %>% 
    as_tibble() %>% 
    clean_names() 
  
  pollock_spec <- pollock_specimen_orig %>% 
    filter(hauljoin %in% hauljoins) %>% 
    right_join(hauls_survey_dat, by = c("cruisejoin", "hauljoin", "region", "vessel", "cruise", "haul")) %>% 
    mutate(age_strata = if_else(stratum %in% c(10, 31, 32, 50), 1, 2)) 
  
  return(pollock_spec)
}

catch_data_d <- function(hauljoins)
{
  query_command <- paste0(" SELECT * FROM racebase.catch
                           WHERE region='BS' and species_code = 21740")
  pollock_catch <- sqlQuery(channel, query_command) %>% 
    as_tibble() %>% 
    clean_names() %>% 
    filter(hauljoin %in% hauljoins)
  return(pollock_catch)
}

length_data_d <- function(hauljoins)
{
  query_command <- paste0(" SELECT * FROM racebase.length
                           WHERE region='BS' and species_code = 21740")
  pollock_length <- sqlQuery(channel, query_command) %>% 
    as_tibble() %>% 
    clean_names() %>% 
    filter(hauljoin %in% hauljoins)
  
  
  query_command <- paste0(" select c.vessel_id vessel, d.cruise, b.haul, a.LENGTH_ID, a.HAUL_ID, a.SPECIES_CODE, a.SEX, a.LENGTH,
                            a.FREQUENCY, b.cruise_id 
                            from race_data.lengths a, race_data.hauls b, race_data.vessels c, race_data.cruises d
                            where b.cruise_id in (", cruise_id, ") and a.haul_id=b.haul_id 
                            and a.species_code in (21740,21741)
                            and c.vessel_id in (", vessel_code, ") and c.vessel_id=d.vessel_id 
                            and b.cruise_id=d.cruise_id
                            order by a.length_id;")
  
  pollock_raw_length <- sqlQuery(channel, query_command) %>% 
    as_tibble() %>% 
    clean_names()
  
  # CIA NOTEs (from Rebecca H.): in racebase.length, juv pollock is extrapolated; racedata.length, juv pollock is NOT extrapolated
  #   racebase.length_frequency, nothing is extrapolated out to haul level, EXCEPT juv pollock
  #   See email from Heather
  
  return(list(pollock_length = pollock_length,
              pollock_raw_length = pollock_raw_length))
}

metadata_d <- function(cruise_id = cruise, vessel_code_id = vessel_code, pollock_specimen_data = pollock_specimen, meta_select = strat_meta_year)
{
  # get most recent year of metadata (most recent cruise and vessel)
  query_command <- paste0(" SELECT * FROM racebase.cruise WHERE cruise in (",
                          cruise_id, ") and region='BS' and vessel in(", vessel_code_id, "); ")
  ebs_shelf_cruise <- sqlQuery(channel, query_command) %>% 
    as_tibble() %>% 
    clean_names() %>% 
    mutate(year = year(start_date)) %>% 
    mutate(vessel = as.numeric(vessel)) %>%
    dplyr::select(vessel, year) %>% 
    distinct()
  
  # ensure most recent cruise info is included in metadata:
  pollock_no_fpc <- pollock_specimen_data %>% 
    dplyr::select(vessel, year) %>% 
    distinct() %>% 
    full_join(ebs_shelf_cruise) %>% 
    mutate(species_code = 21740, fpc = 1)
  
  # general metadata
  query_command <- paste0(" SELECT * FROM racebase.cruise
                        WHERE region='BS' and agency_name='USA' and vessel in(", vessel_code_id, ")")
  metadata <- sqlQuery(channel, query_command)
  
  # strata metadata
  query_command <- paste0(" SELECT * FROM racebase.stratum WHERE region='BS'") #BS incl BS and NBS 81, 70, 71 = NBS
  strata_metadata_raw <- sqlQuery(channel, query_command) %>% 
    as_tibble() %>% 
    clean_names()
  
  #  NOTE: any 999 = stratum is a cumulative number across strata
  
  strata_metadata_2010 <- strata_metadata_raw %>%  # compare to stan, uses 2010
    dplyr::filter(year == 2010) %>% 
    dplyr::filter(stratum < 100) %>% #100-160 not part of normal survey grid area
    mutate(subarea = if_else(stratum %in% NBS_subarea, 0, as.numeric(sub('^(.).*(.)$', '\\1', stratum))))
  #   NBS = subarea 0; other subareas are the first digit of the stratum
  
  strata_metadata_2019 <- strata_metadata_raw %>% # 2019 update
    dplyr::filter(year == 2019)  %>% 
    dplyr::filter(stratum < 100) %>% #100-160 not part of normal survey grid area
    mutate(subarea = if_else(stratum %in% NBS_subarea, 0, as.numeric(sub('^(.).*(.)$', '\\1', stratum))))
  
  strata_metadata_2022 <- strata_metadata_raw %>% # 2022 update
    dplyr::filter(year == 2022)  %>% 
    dplyr::filter(stratum < 100) %>% #100-160 not part of normal survey grid area
    mutate(subarea = if_else(stratum %in% NBS_subarea, 0, as.numeric(sub('^(.).*(.)$', '\\1', stratum))))
  
  if(meta_select == 2010) {strata_metadata <- strata_metadata_2010
  } else if(meta_select == 2019) {strata_metadata <- strata_metadata_2019
  } else if(meta_select == 2022) {strata_metadata <- strata_metadata_2022
  } else(print("Invalid selection: please select strata from 2010, 2019, or 2022"))
  
  return(strata_metadata)
  
}

