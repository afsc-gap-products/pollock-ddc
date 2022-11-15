# This function produces age-length keys for the Bering Sea Pollock assessment
# This code is converted from agecomp2key.SQL files by Stan Kotwicki
# Author: Caitlin I. Allen Akselrud
# Contact: caitlin.allen_akselrud@noaa.gov
# Date created: 2021.06.01
# Date updated: 2021.06.01

# agecomp2key.SQL conversion
# 
# need: pollock_specimen, pollock_sizecomp4 

get_al_key <- function(specimen, length_comp)
{
  # specimen <- pollock_specimen
  # length_comp <- ddc_length_comps$pollock_length_comp
  
  # age_2key1
  al_freq <- specimen %>% 
    group_by(year, age_strata, sex, length, age) %>% 
    summarise(age_freq = n(), 
              avg_weight = mean(weight, na.rm = T))
  
  # age_2key_sum_length
  length_freq <- specimen %>% 
    group_by(year, age_strata, sex, length) %>% 
    summarise(length_freq = n(), 
              avg_weight = mean(weight, na.rm = T))
  
  # age_2key2
  al_freq_prop <- full_join(al_freq, length_freq, by = c("year", "age_strata", "sex", "length")) %>% 
    mutate(al_prop = age_freq/length_freq) #%>% 
    # select(-avg_weight.x, -avg_weight.y, -length_freq)
  
  # test to make sure proportions are right: summed values = 1, all OK
  # test <- al_freq_prop %>% group_by(year, sex, length, age_strata) %>% 
  #   summarise(sum_al_prop = sum(al_prop)) 
  # 
  # sprintf("%.10f", unique(test$sum_al_prop))
  
  # agecomp1_2key
  # # STAN note: decode returns 1 if b.strata are 1 or 3 or 5, and if not it returns 2
  alk <- length_comp %>% 
    mutate(age_strata = case_when(subarea %in% c(1,3,5) ~ 1L,
                                  !(subarea %in% c(1,3,5)) ~ 2L)) %>% 
    drop_na(year, sex, length, age_strata) %>%
    left_join(al_freq_prop, by = c("year", "sex", "length", "age_strata")) %>%
    mutate(pop_at_age = al_prop* pop_prop_num_ha)
  # CIA: tibble doubles in size, meaning there are two matching rows from al_freq_prop for every row from length_comp
  #     # There can be more than one age for each year, sex, length, and age strata
  
  # agecomp2_2key
  alk_summary_bysubarea <- alk %>% 
    group_by(year, subarea, sex, age) %>% 
    summarise(pop_at_age_sum = sum (pop_at_age), 
           mean_laa = sum(length*pop_at_age)/pop_at_age_sum)
  
  # agecomp3_2key & agecomp4_2key
  alk_sd_bysubarea <- full_join(alk_summary_bysubarea, alk) %>% 
    mutate(sqr = pop_at_age * ((length-mean_laa)^2)) %>% 
    group_by(year, subarea, sex, age) %>% 
    summarize(pop_at_age_total = sum(pop_at_age),
           mean_pop_length = sum(length*pop_at_age_total)/pop_at_age_total,
           sd_pop_length = sqrt(sum(sqr)/pop_at_age_total))
  # 
  
  # agecomp5_2key
  alk_overall <- alk %>% 
    group_by(year, sex, age) %>% 
    summarise(pop_at_age_sum = sum (pop_at_age), 
              mean_laa = sum(length*pop_at_age)/pop_at_age_sum)
  
  
  # agecomp6_2key & agecomp7_2key
  
  alk_sd_overall <- full_join(alk_overall, alk) %>% 
    mutate(sqr = pop_at_age * ((length-mean_laa)^2)) %>% 
    group_by(year, sex, age) %>% 
    summarize(pop_at_age_total = sum(pop_at_age),
              mean_pop_length = sum(length*pop_at_age_total)/pop_at_age_total,
              sd_pop_length = sqrt(sum(sqr)/pop_at_age_total))
  
  
  # agecomp8_2key; no nbs or strata 90, 92- northwest extension
  # # CIA: skipping this
  # alk_no_n_extension
  
  # agecomp9_2key
  
  
  # agecomp10_2key
  
  
  return(age_length_key = alk_sd_overall)
}