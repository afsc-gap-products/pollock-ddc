# These functions produce a density-dependent corrected variance-covariance matrix index for the Bering Sea Pollock assessment
# This code is converted from new4.R and new5.R files by Stan Kotwicki
# Author: Caitlin I. Allen Akselrud
# Contact: caitlin.allen_akselrud@noaa.gov
# Date created: 2021.08.23
# Date updated: 2021.08.23

# get density-dependent corrected variance-covariance matrix

# one cpue estimate for all years and all strata
#   keeping all years and all stations, randomly ordered, with some stations duplicated and some not selected

resamp_data <- function(corr_cpue)
{
  # corr_cpue <- ddc_table
  
  resamp <- corr_cpue %>% 
    group_by(year, stratum) %>% #for each year and stratum
    modelr::resample(idx = 1:dim(corr_cpue)[1]) #resample the rows randomly, allowing for some to be duplicated and some not selected
  
  return(resample = resamp$data)
}

resamp_pop <- function(resamp, strata)
{
  # resamp <- resample
  # strata <- strata_metadata
  
  strata <- strata %>% dplyr::select(stratum, area) %>% 
    mutate(area_ha = area*100)
  
  numbers <- resamp %>%  
    group_by(year, stratum) %>% 
    summarise(mean_kg_ha_ddc = mean(ddc_cpue_kg_ha, na.rm = TRUE), #mean by year and stratum
              mean_num_ha_ddc = mean(ddc_cpue_num_ha, na.rm = TRUE)) %>% 
    full_join(strata) %>% #add area fished
    mutate(ddc_kg = mean_kg_ha_ddc * area_ha, #multiply kg or num by area to get total 
           ddc_num = mean_num_ha_ddc * area_ha) %>% 
    ungroup() %>% 
    group_by(year) %>% 
    summarise(total_ddc_kg = sum(ddc_kg, na.rm = TRUE), # get total pop for each year over all strata
              total_ddc_num = sum(ddc_num, na.rm = TRUE),
              total_ddc_th_t = sum(ddc_kg, na.rm = TRUE)*0.001*0.001) #convert kg to th tons
  
  return(resamp_pop = numbers)
  
}




