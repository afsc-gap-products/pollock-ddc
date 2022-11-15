# This function produces proportional areas of strata for the Bering Sea Pollock assessment
# This code is converted from part of population.SQL files by Stan Kotwicki
# Author: Caitlin I. Allen Akselrud
# Contact: caitlin.allen_akselrud@noaa.gov
# Date created: 2021.05.20
# Date updated: 2021.05.20

# strata proportions from population.sql conversion


get_stratum_proportions <- function(ebs_strata)
{
  # CIA: check all for consistency/inconsistency with Stan's tables
  
  # ebs_strata <- strata_metadata
  ebs_subarea <- ebs_strata %>% 
    group_by(region, year, subarea) %>% 
    summarise(tot_area = sum(area)) %>%  #note: area = km^2
    mutate(tot_area_ha = tot_area * 100) # area by hectares
  
  ebs_stratum_prop <- full_join(ebs_strata, ebs_subarea) %>% 
    mutate(prop_area = area/tot_area) %>% 
    arrange(stratum)
  
  return(ebs_stratum_prop = ebs_stratum_prop)
}