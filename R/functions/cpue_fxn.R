# This function produces cpue estimates for the Bering Sea Pollock assessment
# This code is converted from cpue_no_fpc.SQL files by Stan Kotwicki
# Author: Caitlin I. Allen Akselrud
# Contact: caitlin.allen_akselrud@noaa.gov
# Date created: 2021.05.05
# Date updated: 2021.05.05

# cpue_no_fpc.sql conversion
#    this script calculates mean substrata cpue 'avg_cpue' in kg/ha, and numbers/ha,
#    it also calculates variance of these estimates, count of hauls, and area fished.
#    It is important to remember that strata that substrata 82 was never sampled entirely
#      
#    table 'hauls_survey' gets all EBS shelf survey tows starting in year 1982,
#    add field 'year" to identify year of the survey

get_cpue <- function(hauls, catch, strata_filter)
{
  
  addl_vars <- hauls %>% 
    dplyr::select(year, stratum, any_of(names(catch)), - auditjoin) #keep any column names that match catch data, except auditjoin which causes join explosions
  
  # pollock catch: cruisejoin, hauljoin, catchjoin, region, substrata (stratum), vessel, cruise, year, haul, species_code, weight, number_fish, subsample_code, voucher, auditjoin
  poll_catch <- catch %>% 
    dplyr::filter(hauljoin %in% hauls$hauljoin) %>% 
    full_join(addl_vars)  # catch with stratum and year
    
    # fill in NAs in categories; fill in 0s for NAs; if there is weight but no number, should be NA
    # CIA: you are here
    
    # don't need to change weight and number of fish-- not using FPC
  
  # pollock cpue: vessel, cruise, year, haul, substrata (stratum), latitude, longitude, species_code, 
  #     cpue_kgha (weight/(dist fished * net width/10)), cpue_noha (num fish/(dist fished * net width/10)), area_fished (dist fished * net width/10)
  poll_cpue <- hauls %>% 
    dplyr::select(vessel, cruise, year, haul, stratum, start_latitude, start_longitude, distance_fished, net_width, cruisejoin, hauljoin, gear_temperature) %>% 
    left_join(dplyr::select(poll_catch, cruisejoin, hauljoin, vessel, cruise, year, haul, stratum, weight, number_fish, species_code)) %>% 
    #dplyr::filter(!is.na(weight)) %>% 
    mutate(area_fished_ha = distance_fished*net_width/10) %>% 
    mutate(cpue_kg_ha = weight/area_fished_ha,
           cpue_num_ha = number_fish/area_fished_ha) %>% 
    mutate(cpue_kg_ha = if_else(is.na(cpue_kg_ha), 0, cpue_kg_ha)) %>% 
    mutate(cpue_num_ha = if_else(is.na(cpue_num_ha), 0, cpue_num_ha)) %>% 
    mutate(cpue_num_ha = if_else(cpue_num_ha == 0 & cpue_kg_ha > 0, NA_real_, cpue_num_ha)) %>%
    dplyr::filter(stratum %in% unique(strata_filter$stratum)) %>% 
    arrange(cruise, vessel, haul)
  
  # average cpue: avg and variance cpue by strata
  #   year, stratum, avg_cpue_no, var_cpue_no, avg_cpue_kg, var_cpue_kg, count_hauls, area_fished
  avg_cpue <- poll_cpue %>% 
    ungroup() %>% 
    group_by(year, stratum) %>% 
    summarize(number_hauls_no = length(cpue_num_ha[!is.na(cpue_num_ha)]), # this is to account for NAs in rows with weight but not lengths
              number_hauls_kg = length(cpue_kg_ha[!is.na(cpue_kg_ha )]),
              avg_cpue_no_ha = mean(cpue_num_ha, na.rm = TRUE), 
              var_cpue_no_ha = var(cpue_num_ha,  na.rm = TRUE)/number_hauls_no, 
              avg_cpue_kg_ha = mean(cpue_kg_ha,  na.rm = TRUE), 
              var_cpue_kg_ha = var(cpue_kg_ha,   na.rm = TRUE)/number_hauls_kg, 
              tot_area_fished = sum(area_fished_ha, na.rm = TRUE))
  
  return(list(p_catch = poll_catch, p_cpue = poll_cpue, avg_cpue_yr_strat = avg_cpue))
  
}
