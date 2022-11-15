# This function produces length compositions for the Bering Sea Pollock assessment
# This code is converted from sizecomp.SQL files by Stan Kotwicki
# Author: Caitlin I. Allen Akselrud
# Contact: caitlin.allen_akselrud@noaa.gov
# Date created: 2021.05.11
# Date updated: 2021.12.16

# sizecomp.SQL conversion
# this script calculates pollock size composition for entire survey area
# table pollock_length joins racebase. Length table with hauls_survey table
  

length_comp_f <- function(length, hauls, catch, stratum, biomass_kg, biomass_mt)
{
  # length <- pollock_length
  # hauls <- hauls_survey
  # catch <- pollock_catch
  # stratum <- pop_info$ebs_strata
  # biomass_kg <- pop_info$pollock_biomass_kg_ha
  # biomass_mt <- pop_info$pollock_biomass_MT_ha
  
  options(dplyr.summarise.inform = FALSE)
  
  length_info <- left_join(dplyr::select(length, cruisejoin, hauljoin, catchjoin, region, vessel, length, frequency, sex, sample_type, length_type),
                           dplyr::select(hauls, cruisejoin, hauljoin, region, vessel, year, haul, distance_fished, net_width, start_latitude, start_longitude, stratum, stationid)) %>%
    dplyr::filter(!is.na(year)) %>% 
    arrange(year, vessel, haul)
  # length_info <- left_join(dplyr::select(length, hauljoin, catchjoin, length, frequency, sex, sample_type, length_type, auditjoin), 
  #                          dplyr::select(hauls, cruisejoin, hauljoin, region, vessel, year, haul, distance_fished, net_width, start_latitude, start_longitude, stratum, stationid)) %>% 
  #   arrange(year, vessel, haul)
  
  total_measured <- length_info %>% 
    group_by(hauljoin) %>% 
    summarize(tot_measured = sum(frequency))
  
  # table pollock_length2 appends total_measured, 
  #   and number_fish (from catch) to pollock length table
  length_info_2 <- full_join(length_info, total_measured) %>% 
    left_join(dplyr::select(catch, cruisejoin, hauljoin, catchjoin, region, vessel, haul, number_fish)) %>% 
    arrange(year, vessel, haul)
  
  
  # pollock_sizecomp calculates cpue at length for each station
  pollock_sizecomp <- length_info_2 %>% 
    mutate(cpue_length_num_ha = ((frequency/tot_measured)*number_fish)/(distance_fished * net_width/10)) #%>% #distance_fished = km; net_width = m; then convert to ha: m -> km = *1000, km -> ha: /100; => 1000/100 = 10; (amount)/10
  # note: lengths are in mm
    # CIA: the following lines are a possible temp fix for the duplicate data rows problem with diff frequencies
    # group_by(cruisejoin, hauljoin, catchjoin, length, sex, sample_type, length_type,  year,  haul, distance_fished, net_width, start_latitude, start_longitude, stratum, stationid, tot_measured, number_fish) %>%
    # summarize(frequency = sum(frequency),
    #           cpue_length_num_ha = sum(cpue_length_num_ha ))
  
  haul_count <- hauls %>% 
    group_by(year, stratum) %>% 
    summarize(haul_count = length(haul))
  
  # pollock_sizecomp2 sums at length cpue by stratum
  pollock_sizecomp_summary <- pollock_sizecomp %>% 
    group_by(year, stratum, sex, length) %>% 
    summarize(sum_cpue_length = sum(cpue_length_num_ha))
  
  # pollock_sizecomp3 calculates average cpue at length in stratum 
  pollock_sizecomp_sumbyhaul <- full_join(pollock_sizecomp_summary, haul_count) %>% 
    mutate(avg_cpue_length = sum_cpue_length/haul_count)
  
  # weighted_cpue_length weights cpue at length by stratum area ratio
  stratum_use <- stratum %>% dplyr::select(stratum, subarea, prop_area)
  weighted_cpue_length <- pollock_sizecomp_sumbyhaul %>% 
    left_join(stratum_use) %>% 
    mutate(cpue_ratio_no_ha = prop_area * avg_cpue_length)
  
  # mean_cpue_length calculates mean cpue at length in subarea
  mean_cpue_length <- weighted_cpue_length %>%
    group_by(year, subarea, sex, length) %>% 
    summarise(cpue_subarea_ratio_no_ha = sum(cpue_ratio_no_ha, na.rm=TRUE))
  
  # sum_mean_cpue sums all cpues within the strata, to be used to 
  #   calculate proportion of cpue at length to sum of all cpue prop_cpue_length, 
  #   this proportion needs to be used in the events when in some haul fish 
  #   has been caught but not measured
  sum_mean_cpue <- mean_cpue_length %>% 
    group_by(year, subarea) %>% 
    summarise(sum_cpue_subarea_no_ha = sum(cpue_subarea_ratio_no_ha, na.rm = TRUE))
  
  prop_cpue_length <- full_join(mean_cpue_length, sum_mean_cpue) %>% 
    mutate(cpue_prop = cpue_subarea_ratio_no_ha/sum_cpue_subarea_no_ha)
  
  # pollock_sizecomp4 used proportion calculated above and total popolation estimate 
  #   from pollock_biomass, to calculate sizecomp in the event of hauls 
  #   with polllock catch but no length sample
  pollock_sizecomp_pop <- prop_cpue_length %>% 
    dplyr::select(-cpue_subarea_ratio_no_ha, -sum_cpue_subarea_no_ha ) %>% 
    left_join(dplyr::select(biomass_kg, year, subarea, population_num_ha )) %>% 
    mutate(pop_prop_num_ha = cpue_prop * population_num_ha)
  # CIA: check some of the joins here- problems with some subareas as NAs
  
  # CIA: added section for cpue biomass at length (kg/ha) (on 12/16)
  # pollock_sizecomp_pop_weight <- prop_cpue_length %>% 
  #   dplyr::select(-cpue_subarea_ratio_no_ha, -sum_cpue_subarea_no_ha ) %>% 
  #   left_join(dplyr::select(biomass_kg, year, subarea, biomass_kg_ha )) %>% 
  #   mutate(pop_prop_num_ha = cpue_prop * biomass_kg_ha)
  
  
  # pollock_sizecomp_st_survey sums strata sizecomp into one sizecomp 
  #   for STANDARD survey area
  # NOT USED
  pollock_sizecomp_standard_survey <- pollock_sizecomp_pop %>% 
    dplyr::filter(subarea < 7) %>% # CIA: WHY IS THIS?? standard survey area? NOTE that NBS subarea assignment = 0; so this incl NBS FYI
    dplyr::filter(subarea != 0) %>% # this EXCLUDES NBS, so you are only seeing *standard* survey area
    group_by(year, sex, length) %>% 
    summarise(pop_prop_sum = sum(pop_prop_num_ha, na.rm = TRUE))
  
  # pollock_sizecomp_nw_survey sums strata sizecomp into one sizecomp 
  #   for survey area with nw stations, starting in 1987
  # NOT USED
  pollock_sizecomp_nw_survey_noNBS <- pollock_sizecomp_pop %>% 
    dplyr::filter(year > 1986) %>% 
    dplyr::filter(subarea != 0) %>% # this EXCLUDES NBS, so you are only seeing EBS (NW) survey area
    group_by(year, sex, length) %>% 
    summarise(pop_prop_sum = sum(pop_prop_num_ha, na.rm = TRUE))
  
  # population checks
  pop_check <- pollock_sizecomp_pop %>%
    group_by(year, subarea) %>%
    summarise(pop_prop_sum = sum(pop_prop_num_ha, na.rm = TRUE))
  
  pop_check_yr <- pollock_sizecomp_pop %>%
    group_by(year) %>%
    summarise(pop_prop_sum = sum(pop_prop_num_ha, na.rm = TRUE))
  
  orig_pop_yr <- biomass_kg %>% 
    group_by(year) %>%
    summarise(pop_summary = sum(population_num_ha, na.rm = TRUE))
  
  print(paste("********************************************"))
  print(paste("PLEASE PAY ATTENTION TO THE FOLLOWING LINES:"))
  print(paste("********************************************"))
  print(paste("Checking for differences btw calc values and expected values: ", anti_join(pop_check_yr, orig_pop_yr)))
  print(paste("Checking for differences btw calc values and expected values: ", anti_join(orig_pop_yr, pop_check_yr)))
  print(paste("If your difference results are: NA, 0, integer(0), or numeric(0), you have passed this check"))
  
  return(list(pollock_length_comp = pollock_sizecomp_pop, 
              # pollock_length_comp_biomass = pollock_sizecomp_pop_weight,
              sizecomp_cpue_stn = pollock_sizecomp)) #CIA: which tables needed to return vs. needed as output to .csv?
}