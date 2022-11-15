# This function produces population and biomass estimates by strata for the Bering Sea Pollock assessment
# This code is converted from population.SQL files by Stan Kotwicki
# Author: Caitlin I. Allen Akselrudi
# Contact: caitlin.allen_akselrud@noaa.gov
# Date created: 2021.05.04
# Date updated: 2021.05.11

# population.sql conversion
#   this script calculates biomass and population by strata in EBS 
#   and variance of these estimates. 
#   substrata 70 and 81 were excluded from these calculations

population <- function(ebs_strata, avg_cpue)
{
  
  # CIA: check all for consistency/inconsistency with Stan's tables
  
  # ebs_strata <- nw_strata
  # avg_cpue <- cpue_info$avg_cpue_yr_strat
  # 
  strat_area <- ebs_strata %>% 
    distinct(subarea, .keep_all = TRUE) %>% 
    ungroup() %>% 
    dplyr::select(subarea, tot_area, tot_area_ha) #note: tot_area = km
  
  
  # cpue weighted by stratum area (distance weighted)
  weighted_cpue <- avg_cpue %>% 
    full_join(dplyr::select(ebs_strata, stratum, subarea, prop_area)) %>% 
    mutate(cpue_no_dist_wt = avg_cpue_no_ha * prop_area,
           cpue_no_var_dist_wt = var_cpue_no_ha * (prop_area^2), # why prop area ^2? area squared
           cpue_kg_dist_wt = avg_cpue_kg_ha * prop_area,
           cpue_kg_var_dist_wt = var_cpue_kg_ha * (prop_area^2)) # why prop area ^2? area squared
  
  total_weighted_cpue_subarea <- weighted_cpue %>% 
    group_by(year, subarea) %>% 
    summarise(tot_num_hauls = sum(number_hauls_kg),
              tot_cpue_no_ha_dist = sum(cpue_no_dist_wt, na.rm = TRUE),
              tot_cpue_no_ha_var_dist = sum(cpue_no_var_dist_wt, na.rm = TRUE),
              tot_cpue_kg_ha_dist = sum(cpue_kg_dist_wt, na.rm = TRUE),
              tot_cpue_kg_ha_var_dist = sum(cpue_kg_var_dist_wt, na.rm = TRUE))
  

  pollock_biomass_MT_ha <- total_weighted_cpue_subarea %>%
    left_join(strat_area) %>% 
    mutate(biomass_MT_ha = tot_area * tot_cpue_kg_ha_dist /10, # why /10? km -> ha = *100; kg ->MT = /1000; => 100/1000 = /10
           biomass_th_t = biomass_MT_ha/1000,
           biomass_MT_ha_var = ((tot_area/0.01)^2) * tot_cpue_kg_ha_var_dist /1000000, # why /1000000? for kg -> MT variance: /(1000^2)
           population_num_ha = tot_area * tot_cpue_no_ha_dist  * 100, # why *100? km -> ha = *100
           population_num_ha_var = ((tot_area/0.01)^2) * tot_cpue_no_ha_var_dist ) # why no adjustment? km -> ha = *100 OR tot_area/0.01
  
  pollock_biomass <- total_weighted_cpue_subarea %>%
    left_join(strat_area) %>% 
    mutate(biomass_kg_ha = tot_area_ha * tot_cpue_kg_ha_dist, # biomass in kg
           biomass_kg_ha_var = (tot_area_ha^2) * tot_cpue_kg_ha_var_dist,
           population_num_ha = tot_area_ha * tot_cpue_no_ha_dist, 
           population_num_ha_var = (tot_area_ha^2) * tot_cpue_no_ha_var_dist)
    

# empirical distribution function (edf) -----------------------------------

  #step 1: fi table: strata, subarea, year, counth, area fished, calc: fi
  fi <- dplyr::select(avg_cpue, year, stratum, number_hauls_kg, tot_area_fished) %>% 
    left_join(dplyr::select(ebs_strata, stratum, subarea, tot_area)) %>% 
    mutate(fi = ((tot_area/tot_area_fished)*number_hauls_kg) * 
             (((tot_area/tot_area_fished)*number_hauls_kg)-number_hauls_kg)/number_hauls_kg)
  
  # CIA: you are here: number of hauls and total area fished is less than Stan's table (even incl NBS...??)
  
  # step 2: create edf 1
  edf_1 <- dplyr::select(avg_cpue, year, stratum, var_cpue_kg_ha, var_cpue_no_ha) %>% 
    full_join(fi) %>% 
    mutate(fi_kg = fi*var_cpue_kg_ha*number_hauls_kg,
           fi_kg_2 = (fi*var_cpue_kg_ha*number_hauls_kg)^2,
           fi_num = fi*var_cpue_no_ha*number_hauls_kg, 
           fi_num_2 =(fi*var_cpue_no_ha*number_hauls_kg)^2)
    
  
  # step 3: final edf table
  edf <- edf_1 %>% 
    dplyr::filter(fi_kg_2 != 0) %>% 
    group_by(year, subarea) %>% 
    summarize (num_hauls_subarea = sum(number_hauls_kg),
               edf_biomass = (sum(fi_kg)^2)/(sum(fi_kg_2/(number_hauls_kg - 1))),
               edf_pop =  (sum(fi_num)^2)/(sum(fi_num_2/(number_hauls_kg - 1)))) 
  # CIA: very different values from Stan's tables...
    
  
# confidence intervals ----------------------------------------------------
  
  # use t-dist function: qt(p = prob, df = edf_biomass); 
  #   range (low to high): biomass - 95% CI to biomass + 95% CI
  # CIA: NOTE- IMPORTANT: using qt() function results in more precise t-stat values than those in SLQ t-table
  conf_interval_biomass_kg <- dplyr::select(pollock_biomass, year, subarea, 
                                            biomass_kg_ha, biomass_kg_ha_var, 
                                            population_num_ha, population_num_ha_var) %>% 
    full_join(edf) %>% 
    mutate(df_bm = floor(round(edf_biomass, digits = 6)), #floor goes down to eg. 3 if number is 3.9999999, but displays as 4.
           df_num = floor(round(edf_pop, digits = 6)),
           t_stat_bm = qt(p = 0.975, df = df_bm),    # need p = 0.975 for 95% CI: 1-0.5/2
           t_stat_num = qt(p = 0.975, df = df_num),  # need p = 0.975 for 95% CI: 1-0.5/2
           low95_bm = biomass_kg_ha - t_stat_bm,
           high95_bm = biomass_kg_ha + t_stat_bm,
           low95_num = population_num_ha - t_stat_num,
           high95_num = population_num_ha + t_stat_num) 
  
  conf_interval_biomass_MT <- dplyr::select(pollock_biomass_MT_ha, year, subarea, 
                                            biomass_MT_ha, biomass_MT_ha_var, 
                                            population_num_ha, population_num_ha_var) %>% 
    full_join(edf) %>% 
    mutate(df_bm = floor(round(edf_biomass, digits = 6)), #floor goes down to eg. 3 if number is 3.9999999, but displays as 4.
           df_num = floor(round(edf_pop, digits = 6)),
           t_stat_bm = qt(p = 0.975, df = df_bm),    # need p = 0.975 for 95% CI: 1-0.5/2
           t_stat_num = qt(p = 0.975, df = df_num),  # need p = 0.975 for 95% CI: 1-0.5/2
           low95_bm = biomass_MT_ha - t_stat_bm,
           high95_bm = biomass_MT_ha + t_stat_bm,
           low95_num = population_num_ha - t_stat_num,
           high95_num = population_num_ha + t_stat_num) 
    
  all_pollock_biomass <- full_join(pollock_biomass, conf_interval_biomass_kg)
  
  all_pollock_biomass_MT_ha <- full_join(pollock_biomass_MT_ha, conf_interval_biomass_MT)
  
  return(list(pollock_biomass_kg_ha = all_pollock_biomass, 
              pollock_biomass_MT_ha = all_pollock_biomass_MT_ha,
              ebs_strata = ebs_strata))
  
}
