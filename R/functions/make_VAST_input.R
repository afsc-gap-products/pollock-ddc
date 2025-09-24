# Make data formatted for VAST input
# Author: Caitlin I. Allen Akselrud
# Maintained by: Sophia N. Wassermann
# Contact: sophia.wassermann@noaa.gov
# Date created: 2022.06.02
# Date updated: 2024.03.13

make_VAST_input <- function(hauls = hauls_survey,
                            spec = pollock_specimen,
                            ddc_index = ddc_table,
                            ddc_age = ddc_alk,
                            slope_survey = slope_survey)
{
  # add in 2018 NBS emergency survey (read_csv)
  # ddc_table_2018 <- read_csv(here("output", "VAST_ddc_2018_NBS_addon.csv")) %>% 
  #   dplyr::select(year, stratum, start_latitude, start_longitude, ddc_cpue_kg_ha, hauljoin, cruisejoin) %>% 
  #   drop_na(year, start_latitude, start_longitude) 
  # 
  # agecomp_2018 <- read_csv(here("output", "VAST_ddc_alk_2018_NBS_addon.csv")) %>% 
  #   dplyr::select(year, cruisejoin, haul, age, cpue_age_num, prop_check, 
  #                 mean_length, count_n, ddc_cpue_num_ha, cpue_num_ha, correction, age_cpue_corr, )
  # 
  # all stations- add in 0s
  all_stns_surveyed <- hauls %>% 
    # bind_rows(hauls_survey_18) %>% # add 2018 hauls
    dplyr::select(year, cruisejoin, hauljoin, start_latitude, start_longitude, stratum) %>% 
    distinct()
  
  # check this is all unique stations (tibble should be 0)
  # all_stns_surveyed %>% group_by(year, cruisejoin, hauljoin) %>% 
  #   summarise(n = n()) %>% 
  #   dplyr::filter(n>1)
  
  # catch (cpue, kg/ha), year, lat, long
  VAST_ddc_table <- ddc_index %>% 
    ungroup %>% 
    dplyr::select(year, stratum, start_latitude, start_longitude, ddc_cpue_kg_ha, hauljoin, cruisejoin) %>% 
    drop_na(year, start_latitude, start_longitude) %>% 
    # bind_rows(ddc_table_2018) %>% # bind 2018 rows
    arrange(year) %>% 
    full_join(all_stns_surveyed) %>% 
    mutate(ddc_cpue_kg_ha = if_else(is.na(ddc_cpue_kg_ha), 0, ddc_cpue_kg_ha)) %>% 
    arrange(year) #%>% 
  # dplyr::select(-hauljoin)
  
  VAST_ddc_table_EBS <- VAST_ddc_table %>% 
    filter(!stratum %in% NBS_subarea)  #Not NBS
  
  VAST_ddc_table_NBS<- VAST_ddc_table %>% 
    filter(stratum %in% NBS_subarea) #just NBS
  
  # check for years with 100% encounter rate
  all_years_surveyed <- hauls %>% distinct(year)
  enc_rate <- full_join(hauls, spec, by = "hauljoin") %>% 
    dplyr::filter(is.na(cruisejoin.y)) %>% 
    distinct(year.x) %>% 
    rename(year = year.x) %>% 
    anti_join(all_years_surveyed)
  
  if(dim(enc_rate)[1]==0)
  {
    # print("no years with 100% encounter rate")
    encounter_100 <- FALSE
  }else
  {
    # print("there are years with 100% encounter rate")
    encounter_100 <- TRUE
  }
  
  # do age comps
  all_stns_surveyed <- hauls %>% 
    # bind_rows(hauls_survey_18) %>% # add 2018 hauls
    dplyr::select(year, cruisejoin, hauljoin, start_latitude, start_longitude) %>% 
    distinct()
  # check:
  # filter(all_stns_surveyed, !hauljoin %in% all_hauljoins)
  # ddc_alk %>% anti_join(all_stns_surveyed) #%>% ungroup %>% distinct(ddc_cpue_num_ha)
  # ddc_alk %>% filter(is.na(ddc_cpue_num_ha), is.na(cruisejoin))
  
  ddc_alk <- ddc_age %>% filter(!is.na(ddc_cpue_num_ha), !is.na(cruisejoin)) #remove bad stations from age comps
  
  all_ages <- unique(ddc_alk$age) %>% sort()
  age <- as.integer(seq(from = min(all_ages, na.rm = T), to = 15, by = 1))
  fill_zeroes <- all_stns_surveyed %>% 
    expand_grid(age) %>% 
    filter(!cruisejoin %in% slope_survey$cruisejoin)
  
  VAST_ddc_alk_prelim <- ddc_alk %>%
    # bind_rows(agecomp_2018) %>% 
    arrange(year) %>% 
    ungroup() %>% 
    dplyr::select(age_cpue_corr, year, age, start_latitude, start_longitude, cruisejoin) %>% 
    drop_na(year, start_latitude, start_longitude, age) 
  
  # Add plus group 15+
  plus_alk <- VAST_ddc_alk_prelim %>% 
    dplyr::filter(age >= 15) %>% 
    group_by( year, start_latitude, start_longitude, cruisejoin) %>% 
    summarise(age_cpue_corr = sum(age_cpue_corr)) %>% 
    mutate(age = 15)
  less_alk <-  VAST_ddc_alk_prelim %>% 
    dplyr::filter(age < 15)
  
  VAST_ddc_alk_ages <- bind_rows(less_alk, plus_alk) %>%
    arrange(year, start_latitude, start_longitude, age) %>%
    drop_na(age_cpue_corr)
  
  # missing_ages <- anti_join(fill_zeroes, VAST_ddc_alk_ages, by = c("year", "hauljoin", "age")) #"start_latitude", "start_longitude", "age"))
  # 
  # VAST_ddc_alk <- bind_rows(VAST_ddc_alk_ages, missing_ages)
  # odd_stns <- anti_join(VAST_ddc_alk_ages, fill_zeroes) %>%  #%>% distinct(start_latitude, start_longitude)
  #   dplyr::select(-age_cpue_corr, -age) %>% 
  #   distinct()
  # fill_zeroes <- all_stns_surveyed %>% 
  #   bind_rows(odd_stns) %>% 
  #   expand_grid(age) %>% 
  #   filter(!cruisejoin %in% slope_survey$cruisejoin)
  
  VAST_ddc_alk <- VAST_ddc_alk_ages %>% 
    filter(!cruisejoin %in% slope_survey$cruisejoin) %>% 
    full_join(fill_zeroes) %>% 
    mutate(age_cpue_corr = if_else(is.na(age_cpue_corr), 0, age_cpue_corr)) %>% 
    rename(CPUE_num = age_cpue_corr, 
           Year= year, 
           Age = age, 
           Lat = start_latitude,
           Lon = start_longitude) %>% 
    dplyr::select(-hauljoin) %>% 
    arrange(Year, Lat, Lon, Age)
  
  # VAST_ddc_alk_prelim %>% dplyr::filter(is.na(hauljoin))
  return(list(VAST_ddc_table = VAST_ddc_table,
              VAST_ddc_table_EBS = VAST_ddc_table_EBS,
              VAST_ddc_table_NBS = VAST_ddc_table_NBS,
              VAST_ddc_alk = VAST_ddc_alk,
              encounter_rate = encounter_100))
}

# New function for when in-season ages are not available
make_VAST_input_noage <- function(hauls = hauls_survey,
                                  spec = pollock_specimen,
                                  ddc_index = ddc_table,
                                  slope_survey = slope_survey)
{
  # add in 2018 NBS emergency survey (read_csv)
  # ddc_table_2018 <- read_csv(here("output", "VAST_ddc_2018_NBS_addon.csv")) %>% 
  #   dplyr::select(year, stratum, start_latitude, start_longitude, ddc_cpue_kg_ha, hauljoin, cruisejoin) %>% 
  #   drop_na(year, start_latitude, start_longitude) 
  # 
  # agecomp_2018 <- read_csv(here("output", "VAST_ddc_alk_2018_NBS_addon.csv")) %>% 
  #   dplyr::select(year, cruisejoin, haul, age, cpue_age_num, prop_check, 
  #                 mean_length, count_n, ddc_cpue_num_ha, cpue_num_ha, correction, age_cpue_corr, )
  # 
  # all stations- add in 0s
  all_stns_surveyed <- hauls %>% 
    # bind_rows(hauls_survey_18) %>% # add 2018 hauls
    dplyr::select(year, cruisejoin, hauljoin, start_latitude, start_longitude, stratum) %>% 
    distinct()
  
  # check this is all unique stations (tibble should be 0)
  # all_stns_surveyed %>% group_by(year, cruisejoin, hauljoin) %>% 
  #   summarise(n = n()) %>% 
  #   dplyr::filter(n>1)
  
  # catch (cpue, kg/ha), year, lat, long
  VAST_ddc_table <- ddc_index %>% 
    ungroup %>% 
    dplyr::select(year, stratum, start_latitude, start_longitude, ddc_cpue_kg_ha, hauljoin, cruisejoin) %>% 
    drop_na(year, start_latitude, start_longitude) %>% 
    # bind_rows(ddc_table_2018) %>% # bind 2018 rows
    arrange(year) %>% 
    full_join(all_stns_surveyed) %>% 
    mutate(ddc_cpue_kg_ha = if_else(is.na(ddc_cpue_kg_ha), 0, ddc_cpue_kg_ha)) %>% 
    arrange(year) #%>% 
  # dplyr::select(-hauljoin)
  
  VAST_ddc_table_EBS <- VAST_ddc_table %>% 
    filter(!stratum %in% NBS_subarea)  #Not NBS
  
  VAST_ddc_table_NBS<- VAST_ddc_table %>% 
    filter(stratum %in% NBS_subarea) #just NBS
  
  # check for years with 100% encounter rate
  all_years_surveyed <- hauls %>% distinct(year)
  enc_rate <- full_join(hauls, spec, by = "hauljoin") %>% 
    dplyr::filter(is.na(cruisejoin.y)) %>% 
    distinct(year.x) %>% 
    rename(year = year.x) %>% 
    anti_join(all_years_surveyed)
  
  if(dim(enc_rate)[1]==0)
  {
    # print("no years with 100% encounter rate")
    encounter_100 <- FALSE
  }else
  {
    # print("there are years with 100% encounter rate")
    encounter_100 <- TRUE
  }
  
  # do age comps
  all_stns_surveyed <- hauls %>% 
    # bind_rows(hauls_survey_18) %>% # add 2018 hauls
    dplyr::select(year, cruisejoin, hauljoin, start_latitude, start_longitude) %>% 
    distinct()
  # check:
  # filter(all_stns_surveyed, !hauljoin %in% all_hauljoins)
  
  # VAST_ddc_alk_prelim %>% dplyr::filter(is.na(hauljoin))
  return(list(VAST_ddc_table = VAST_ddc_table,
              VAST_ddc_table_EBS = VAST_ddc_table_EBS,
              VAST_ddc_table_NBS = VAST_ddc_table_NBS,
              encounter_rate = encounter_100))
}