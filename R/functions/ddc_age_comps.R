# This function produces density dependent corrected ages at length by year and station for the Bering Sea Pollock assessment
# This code is converted from s1.R files by Stan Kotwicki
# Author: Caitlin I. Allen Akselrud
# Contact: caitlin.allen_akselrud@noaa.gov
# Date created: 2021.09.28
# Date updated: 2021.09.28

# Notes: this is prep for VAST input

ddc_age_comps_f <- function(ddc_cpue, age_key, length_cpue)
{
  # ddc_cpue <- ddc_table
  # age_key <- age_comp_full_key
  # length_cpue <- cpue_length_table #length_comps$sizecomp_cpue_stn
  
  years <- unique(ddc_cpue$year)
  
  # fill aa key
  # # for each year, and each row in each year of the length cpue (aa1)
  # #   get key rows with matching year, sex, and length  (aa2)
  # #   IF there is only one source for the matching key:
  # #     expand that aa1 row for each matching key row (aa2) with: age, proportion, and cpue age = aa2 proportion * aa1 cpue length
  # #   Otherwise: keep aa1 row and add: age = NA, prop = 1, and cpue age = aa1 cpue length
  # #     maybe print an alert when using this option
  # # then summarise: 
  # #    group_by(YEAR, CRUISE, VESSEL, HAUL, AGE)
  # #    cpue_age = sum(cpue_age); prop_check = sum (prop); mean_length = mean(length), count n = n()
  alk <- length_cpue %>% 
    nest_join(age_key) %>% 
    unnest(cols = age_key) #%>% 
    # group_by(year, vessel, haul, stationid, age) #%>% 
    # mutate(n_source = n_distinct(data_source)) #%>%
  
  alk_a <- alk %>% dplyr::filter(data_source == 'a')
  alk_g <- alk %>% dplyr::filter(data_source == 'g')
  alk_na <- alk %>% dplyr::filter(is.na(data_source )) 
  
  # CHECK:
  alk %>% dplyr::filter(is.na(age )) #CIA: NOTE all age == NA have prop = 1, and source = a; no special case treatment needed for NAs
  alk %>% dplyr::filter(is.na(age )) %>% dplyr::select(prop_aal) %>% unique()
  alk %>% dplyr::filter(is.na(age )) %>% dplyr::select(data_source) %>% unique()
  
# join together annual and global sources, so you get vars separated out for same location
  new_alk <-full_join(alk_a, alk_g, by = c("cruisejoin", "vessel", "year", "haul", "start_latitude", 
                                           "start_longitude", "stratum", "stationid", 
                                           "sex", "length", "age")) %>% 
    mutate(cpue_age = if_else(is.na(data_source.x), # cpue_age from annual data unless only global avail for that stn, then cpue_age from global data
                              prop_aal.y*cpue_length_num_ha.y, 
                              prop_aal.x*cpue_length_num_ha.x),
           proportion = if_else(is.na(data_source.x), # cpue_age from annual data unless only global avail for that stn, then cpue_age from global data
                              prop_aal.y, 
                              prop_aal.x))
  alk_summary <- new_alk %>% 
    # group_by(year, haul, stationid, age) %>%
    group_by(year, cruisejoin, haul, start_latitude, start_longitude, age) %>%
    summarise(cpue_age_num = sum(cpue_age), 
              prop_check = sum(proportion), 
              mean_length = mean(length), 
              count_n = n())
  
  # cpue check
  # CIA: sum cpue by station by length to check you get the same cpue sum by year/station
  # CIA: PROP CHECK
  cpue_check <- alk_summary %>% 
    group_by(year, haul) %>% 
    summarise(
      cpue = sum(cpue_age_num),
      prop_check = sum(prop_check),
      n = n()
    )
  
  # correct for density dependence
  # check:
  # cpue %>% filter(is.na(cpue_num_ha ))
  ddc_cpue %>% filter(is.na(ddc_cpue_num_ha ))
  # ddc_corr <- ddc_cpue %>% 
  #   mutate(correction = ddc_cpue_num_ha/cpue_num_ha,
  #          correction = if_else(is.na(correction), 1, correction))
  # hist(ddc_corr$correction)
  # max(ddc_corr$correction)
  
  ddc_corr <- ddc_cpue %>% 
    group_by(year, cruisejoin, haul) %>%
    summarise(ddc_cpue_num_ha = sum(ddc_cpue_num_ha),
              cpue_num_ha = sum(cpue_num_ha)) %>% 
    mutate(correction = ddc_cpue_num_ha/cpue_num_ha,
           correction = if_else(is.na(correction), 1, correction))
  
  hist(ddc_corr$correction)
  max(ddc_corr$correction)
  
  ddc_age <- alk_summary %>% full_join(ddc_corr) %>% 
    mutate(cpue_age_num = if_else(ddc_cpue_num_ha == 0, 0, cpue_age_num),
           age_cpue_corr = cpue_age_num * correction,
           age_cpue_corr = if_else(is.na(age_cpue_corr), ddc_cpue_num_ha, age_cpue_corr))
  
  # check for duplicates:
  # ddc_age %>% ungroup %>%
  #   group_by(year, cruisejoin, haul, age) %>%
  #   filter(n()>1) %>%
  #   arrange(year, cruisejoin, haul, age)
  # table(ddc_age$age)
  # 
  

# age comps by sex --------------------------------------------------------

  alk_summary_sex <- new_alk %>% 
    # group_by(year, haul, stationid, age) %>%
    group_by(year, cruisejoin, haul, start_latitude, start_longitude, age, sex) %>%
    summarise(cpue_age_num = sum(cpue_age), 
              prop_check = sum(proportion), 
              mean_length = mean(length), 
              count_n = n())  
  
  ddc_age_sex <- alk_summary_sex %>% full_join(ddc_corr) %>% 
    mutate(cpue_age_num = if_else(ddc_cpue_num_ha == 0, 0, cpue_age_num),
           age_cpue_corr = cpue_age_num * correction,
           age_cpue_corr = if_else(is.na(age_cpue_corr), ddc_cpue_num_ha, age_cpue_corr))

# return ------------------------------------------------------------------

  
  # return corrected cpue in NUMBERS at age (by year and lat/lon)
  return(list(ddc_age = ddc_age, ddc_age_sex = ddc_age_sex))
  
}
