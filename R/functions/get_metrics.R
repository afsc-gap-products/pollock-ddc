# These functions get data needed for ddc conversion
# Author: Caitlin I. Allen Akselrud
# Contact: caitlin.allen_akselrud@noaa.gov
# Date created: 2022.06.14
# Date updated: 2022.06.14

get_annual_metrics <- function(hauls_ebs = hauls_survey, 
                   hauls_nbs = hauls_survey_nbs, 
                   length = pollock_raw_length,
                   length_nbs = pollock_raw_length_nbs,
                   spec = pollock_specimen,
                   spec_nbs = pollock_specimen_nbs,
                   this_yr = current_year, 
                   nbs_subarea = NBS_subarea)
{
  # num stations (stan does num hauls)
  sort(unique(hauls_ebs$stratum)) # BS and NBS; 81, 70, 71 = NBS; >100 = Russian waters
  num_hauls_tot <- bind_rows(hauls_ebs, hauls_nbs) %>% 
    dplyr::filter(!stratum > 100) %>% #russian waters
    group_by(year, stratum) %>% 
    summarise(num_hauls = n())
  
  EBS_hauls_thisyr <- num_hauls_tot %>% dplyr::filter(year == this_yr) %>% 
    dplyr::filter(!stratum %in% nbs_subarea) %>% 
    summarise(tot_hauls = sum(num_hauls))
  NBS_hauls_thisyr <- num_hauls_tot %>% dplyr::filter(year == this_yr) %>% 
    dplyr::filter(stratum %in% nbs_subarea) %>% 
    summarise(tot_hauls = sum(num_hauls))
  
  # num fish measured (lengths)
  fish_measured_ebs <- length %>% 
    left_join(hauls_ebs) %>% 
    # left_join(hauls_nbs) %>% 
    group_by(year, stratum) %>% 
    summarise(n_measured = n()) %>% 
    na.omit()
  
  fish_measured_nbs <- length_nbs %>% 
    # left_join(hauls_ebs) %>% 
    left_join(hauls_nbs) %>% 
    group_by(year, stratum) %>% 
    summarise(n_measured = n()) %>% 
    na.omit()
  
  EBS_fish_measured_thisyr <- fish_measured_ebs %>% 
    dplyr::filter(!stratum %in% nbs_subarea,
                  year %in% this_yr) %>% 
    summarise(tot_fish = sum(n_measured))
  NBS_fish_measured_thisyr <- fish_measured_nbs %>% 
    dplyr::filter(stratum %in% nbs_subarea,
                  year %in% this_yr) %>% 
    summarise(tot_fish = sum(n_measured))
  
  # num fish aged
  fish_aged <- spec %>% 
    group_by(year, stratum) %>% 
    dplyr::filter(year == this_yr,
                  !is.na(age)) %>% 
    summarise(n_aged = n()) %>% 
    na.omit()
  
  fish_aged_nbs <- spec_nbs %>% 
    group_by(year, stratum) %>% 
    dplyr::filter(year == this_yr,
                  !is.na(age)) %>% 
    summarise(n_aged = n()) %>% 
    na.omit()
  
  EBS_fish_aged_thisyr <- fish_aged %>% 
    dplyr::filter(!stratum %in% nbs_subarea) %>% 
    summarise(tot_fish = sum(n_aged))
  NBS_fish_aged_thisyr <- fish_aged_nbs %>% 
    dplyr::filter(stratum %in% nbs_subarea) %>% 
    summarise(tot_fish = sum(n_aged))
  
  print(paste0(this_yr, " EBS hauls: ", EBS_hauls_thisyr$tot_hauls))
  print(paste0(this_yr, " NBS hauls: ", NBS_hauls_thisyr$tot_hauls))
  
  print(paste0(this_yr, " EBS fish measured: ", EBS_fish_measured_thisyr$tot_fish))
  print(paste0(this_yr, " NBS fish measured: ", NBS_fish_measured_thisyr$tot_fish))
  
  print(paste0(this_yr, " EBS fish aged: ", EBS_fish_aged_thisyr$tot_fish))
  if(dim(NBS_fish_aged_thisyr)[1] > 0)
  {
    print(paste0(this_yr, " NBS fish aged: ", NBS_fish_aged_thisyr$tot_fish))
  }else {print(paste0("There were no fish aged in the NBS"))}
    
  
  return(list(num_hauls_tot = list(EBS_hauls = EBS_hauls_thisyr,
                                   NBS_hauls = NBS_hauls_thisyr),
              fish_measured = list(EBS_fish_measured = EBS_fish_measured_thisyr,
                                   NBS_fish_measured = NBS_fish_measured_thisyr),
              fish_aged = list(EBS_fish_aged = EBS_fish_aged_thisyr,
                               NBS_fish_aged = NBS_fish_aged_thisyr)))
}

