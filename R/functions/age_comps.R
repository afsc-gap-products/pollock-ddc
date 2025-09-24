# This function produces age compositions for the Bering Sea Pollock assessment
# This code is converted from combined agecomp_11.SQL then agecomp_with_global_key_no_loop.SQL files by Stan Kotwicki
# Author: Caitlin I. Allen Akselrud
# Contact: caitlin.allen_akselrud@noaa.gov
# Date created: 2021.05.25
# Date updated: 2021.05.25

# Produce age comps
# agecomp_11.SQL conversion
# agecomp_with_global_key_no_loop.SQL conversion (added to agecomp_11.SQL)
# NOTE: '_g' stands for "global" indicates combined over years

# CIA: generates age-length key and does laa regression

get_agecomps <- function(specimen, pop_lengths)
{
  # specimen <- pollock_specimen_all
  # pop_lengths <- length_comps$pollock_length_comp # ddc_length_comps$pollock_length_comp #aka pollock_sizecomp4
 
  # years
  survey_years <- sort(unique(specimen$year))
  survey_sex <- sort(unique(specimen$sex))
  survey_length <- sort(unique(specimen$length))
  survey_age <- sort(unique(specimen$age))
  full_set <- expand_grid(year = survey_years, sex = survey_sex, length = survey_length, age = survey_age)
  
  # age_key1: number of ages at length, by year and sex
  aal_key <- specimen %>% 
    dplyr::select(year, sex, length, age) %>% 
    group_by(year, sex, length, age) %>%
    count() %>% 
    rename(num_aal = n) %>% 
    arrange(year, sex, length, age)
  
  # weight at age and length, by year and sex
  # extra:
  waal_key <- specimen %>% 
    dplyr::select(year, sex, length, age, weight) %>% 
    group_by(year, sex, length, age) %>%
    summarise(weight_age_length = mean(weight, na.rm = T)) %>% 
    arrange(year, sex, length)
  
  # age_key_global
  aal_key_g <- specimen %>% 
    dplyr::select(sex, length, age) %>% 
    group_by(sex, length, age) %>%
    count(age) %>% 
    rename(num_aal = n) %>% 
    arrange(sex, length, age)
  
  # weight at age and length, by sex
  # extra:
  waal_key_g <- specimen %>% 
    dplyr::select(sex, length, age, weight) %>% 
    group_by(sex, length, age) %>%
    summarise(weight_age_length = mean(weight, na.rm = T)) %>% 
    arrange(sex, length)
  
  # age_sum_at_length number of ages at each length bin (for year and sex)
  yr_num_in_length_bin <- specimen %>%
    dplyr::select(year, sex, length, age) %>% 
    group_by(year, sex, length) %>%
    count() %>% 
    rename(num_in_length_bin = n) %>% 
    arrange(year, sex, length)
  
  # age_sum_at_len_gl
  # age_sum_at_length number of ages at each length bin (for sex)
  num_in_length_bin_g <- specimen %>%
    dplyr::select(sex, length, age) %>% 
    group_by(sex, length) %>%
    count() %>% 
    rename(num_in_length_bin = n) %>% 
    arrange(sex, length)
  
  # age_key2: year, sex, length, age, freq, proportion, source = 's' (for survey?)
  # proportion: what proportion of each length are x age, given year and sex
  
  aal_prop_key <- aal_key %>% 
    left_join(yr_num_in_length_bin) %>% 
    mutate(prop_aal = num_aal/num_in_length_bin) %>% 
    mutate(data_source = 'a') #a = annual 
  
  # age_key_global2
  # sex, length, age, freq, proportion, source = 'g' (for global?)
  # proportion: what proportion of each length are x age, given year and sex
  aal_prop_key_g <- aal_key_g %>% 
    left_join(num_in_length_bin_g) %>% 
    mutate(prop_aal = num_aal/num_in_length_bin) %>% 
    mutate(data_source = 'g') %>% 
    # age_key_global3: expand global for each year
    expand_grid(year = survey_years) %>%
    arrange(year, sex, length, age)
  
  # age_key2_combined: 
  #   where values are missing in aal_prop_key, fill them in with values from aal_prop_key_g
    # anti-join aal annual and global; where annual is missing data, bind_rows from global source
    fill_key <- anti_join(aal_prop_key_g, aal_prop_key, by = c("year", "sex", "length"))#, "age"))
    # fill_key <- right_join(aal_prop_key_g, aal_prop_key, by = c("year", "sex", "length"))
    
    aal_full_key <- bind_rows(aal_prop_key, fill_key) %>% 
      arrange(year, sex, length, age) 

    fill_key_2 <- anti_join(full_set, aal_full_key)
    aal_full_key_allcombos <- bind_rows(aal_full_key, fill_key_2) %>% 
      arrange(year, sex, length, age)
    # table(is.na(aal_full_key$prop_aal))
    # unique(aal_full_key$data_source)
  
  # agecomp1: combine pollock_sizecomp4 and age_key2
  #   year, strata, sex, length, age, fre, prop; pop_at_age = prop * pop; p_sum = population (by year, sex, length)
  # CIA: problematic not to have subareas in aal?
  # CIA: exploding # rows (BUT same in .sql code)
  age_comps <- pop_lengths %>% 
    # ungroup() %>% 
    # dplyr::select(year, subarea, sex, length) %>% 
    dplyr::filter(!is.na(year), !is.na(sex), !is.na(length)) %>% 
    group_by(year, sex, length) %>% 
    left_join(aal_full_key) %>% 
    mutate(pop_aal = prop_aal*pop_prop_num_ha ) %>% 
    mutate(age = if_else(is.na(age),  -9L,  age))#,
           # pop_aal = if_else(age == -9L ,population_num_ha, pop_aal)) # CIA: why does Stan set to -9?? Why not leave as NA?; this line is to set p_sum = pop_aal (not other way around)
  
  # agecomp2: group year, subarea, sex, age, and take mean length
  mean_laa <- age_comps %>% 
    ungroup() %>% 
    group_by(year, subarea, sex, age) %>% 
    summarise(pop_at_age = sum(pop_aal),
              mean_length = sum(length * pop_aal)/ sum(pop_aal))
  
  
  # agecomp3: a.pop_at_age*power(a.length - b.mean_length,2) sqr  
  #   not sure what sqr stands for?
  laa_sqr <- full_join(age_comps, mean_laa) %>% 
    mutate(sqr = pop_aal * ((length - mean_length)^2))
  
  # agecomp4: 
  laa_summary <- laa_sqr %>% 
    group_by(year, subarea, sex, age) %>% 
    summarise(pop_at_age_sum = sum(pop_at_age),
              mean_length_sum = sum(length*pop_at_age)/sum(pop_at_age),
              sd_length = sqrt(sum(sqr)/sum(pop_at_age)))
  
  # agecomp5:
  age_comps_summary <- age_comps %>% 
    group_by(year, sex, age) %>% 
    summarise(pop_sum = sum(population_num_ha),
              mean_length_sum = sum(length*population_num_ha)/sum(population_num_ha))
  
  # agecomp6:
  age_comp_sqr <- full_join(age_comps, age_comps_summary) %>% 
    mutate(sqr_age_comps = pop_aal * ((length-mean_length_sum)^2) )
  
  # agecomp7:
  age_comp_sqr_summary <- age_comp_sqr %>% 
    group_by(year, sex, age) %>% 
    mutate(total_pop_sum = sum(pop_aal),
           total_mean_length_sum = sum(length*pop_aal)/sum(pop_aal),
           total_sd_length = (sum(sqr_age_comps)/sum(pop_aal))^2)
  
 
  # agecomp 8-10 same thing again (as agecomp 5-7), but filtered for EBS subareas (!=0)
  # CIA: check these (stratum included or summarized over- compare to stan's code)
  # agecomp8: 
  age_comps_summary_EBS <- age_comps %>% 
    dplyr::filter(!is.na(subarea)) %>% #filter out NBS and NAs
    group_by(year, sex, age) %>% 
    summarise(pop_sum = sum(population_num_ha),
              mean_length_sum = sum(length*population_num_ha)/sum(population_num_ha))
  
  # agecomp9:
  age_comp_sqr_EBS <- full_join(age_comps, age_comps_summary_EBS) %>% 
    dplyr::filter(!is.na(subarea)) %>% #filter out NBS and NAs
    mutate(sqr_age_comps = pop_aal * ((length-mean_length_sum)^2) )
  
  # agecomp10:
  age_comp_sqr_summary_EBS <- age_comp_sqr_EBS %>% 
    group_by(year, sex, age) %>% 
    summarize(total_pop_sum = sum(pop_aal),
           total_mean_length_sum = sum(length*pop_aal)/sum(pop_aal),
           total_sd_length = (sum(sqr_age_comps)/sum(pop_aal))^2)
  
  # CIA: need to return relevant age comp tables
  # shared ebs_pollock_agecomp_g: laa_summary
  # shared ebs_pollock_agecomp_nw_g: age_comp_sqr_summary_EBS
  # ebs_pollock_cpue_length does not need data from this function; created in main code
  
  return(list(laa_summary_g = laa_summary, 
              laa_summary_g_EBS_only = age_comp_sqr_summary_EBS, 
              aal_filled_key = aal_full_key, 
              aal_key_incl_NAs = aal_full_key_allcombos) )
}
