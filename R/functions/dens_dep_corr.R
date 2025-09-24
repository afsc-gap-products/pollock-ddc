# These functions produces a density-dependent corrected index for the Bering Sea Pollock assessment
# This code is converted from new1.R, new2.R, and new3.R files by Stan Kotwicki
# Author: Caitlin I. Allen Akselrud
# Contact: caitlin.allen_akselrud@noaa.gov
# Date created: 2021.08.13
# Date updated: 2021.08.13

# get density-dependent correction
# new1.R: estimate equivalent backscatter by length for fish we caught
# new2.R: estimate the denisty-dependent corrected CPUE in the form of backscatter
# new3.R: transforms corrected backscatter back into numbers

backscatter_est <- function(length_cpue)
{
  # length_cpue <- cpue_length_table
  
  bs_length <- length_cpue %>% 
    mutate(# convert to miles bc backscatter untis are miles^2
      cpue_mile_sq = 1852^2/10000*cpue_length_num_ha, 
      # formula to estimate how much backscatter comes from 1 fish of given length; target strength; 
      target_strength = 20*log10(length/10)-66,
      #assign the backscatter we expect for each CPUE at length
      backscatter = 4*pi*cpue_mile_sq*10^(target_strength/10) )
  return(bs_length)
}

ddc_cpue_bs <- function(bs_length_equiv, params)
{
  # dd-corrected CPUE in the form of backscatter
  # bs_length_equiv <- equiv_bs_length
  # params <- param_ests
  
  # get sum of backscatter by station using year, vessel, stratum, and haul
  bs_stations <- bs_length_equiv %>% 
    group_by(haul, stratum, vessel, year) %>% 
    summarise(bs_stn_sum = sum(backscatter)) %>% 
    arrange(year, vessel, stratum, haul)# %>% 
    # drop_na()
  # hist(bs_stations$bs_stn_sum, 2000, xlim=c(0,5000))
  # ggplot(bs_stations, aes(x = bs_stn_sum)) +
  #   geom_histogram(binwidth = 100) +
  #   scale_x_continuous(name = "sum of backscatter per station") 
  
  # corrected backscatter
  bs_correction <- bs_stations %>%
    mutate(stn_mean_corr = NA,
           stn_sd_corr = NA) # set up columns to fill in below
  
  for(i in 1:dim(bs_correction)[1])
  {
    sa_corr <- bs_correction$bs_stn_sum[i] / (params$asp + exp(-(params$int + params$slope * bs_correction$bs_stn_sum[i])))
    bs_correction$stn_mean_corr[i] = mean(sa_corr)
    bs_correction$stn_sd_corr[i] = sd(sa_corr)
  } #CIA note- tried a cleaner version w/o loop, but numbers didn't match stan's code
  
  return(backscatter_corr_ests = bs_correction)
}

ddc_fxn <- function(bs_corr, pollock_cpue)
{
  # bs_corr <- backscatter_cpue_corr
  # pollock_cpue <- cpue_info$p_cpue
  
  correction <- bs_corr %>% 
    mutate(eff = bs_stn_sum/stn_mean_corr) %>% # get eff (effect?) of backscatter on numbers
    full_join(pollock_cpue) %>% 
    mutate(ddc_cpue_kg_ha = cpue_kg_ha/eff, #apply ddc
           ddc_cpue_num_ha = cpue_num_ha/eff) %>% 
    mutate(ddc_cpue_kg_ha = if_else(cpue_kg_ha == 0, 0, ddc_cpue_kg_ha), # where cpue is 0, ddc cpue should also be 0 (not NA)
           ddc_cpue_num_ha = if_else(cpue_num_ha == 0, 0, ddc_cpue_num_ha),
           eff = if_else(cpue_kg_ha == 0, 1, eff, missing = eff)) # where cpue is 0, effect = 1 (not NA)
  
  # CIA YOU ARE HERE
  # next step: you will take all places where eff = NA and 
  #   predict that value with a GAM, then feed pred values back into correction
  #   following stan's method (cleaner version didn't work with mgcv::gam package)
  cpue_na <-  data.frame(xx = correction[which(is.na(correction$eff)),]$cpue_kg_ha)
  cpue_na_loc <- which(is.na(correction$eff))
  
  xx <- correction$cpue_kg_ha
  yy <- correction$eff
  plot(yy~xx, xlab = "CPUE (kg/ha)", ylab = "Effect of backscatter")
  eff_gam <- mgcv::gam(yy~s(xx)) #doing gam using smoother s()
  summary(eff_gam)
  plot(eff_gam)
  
  na_eff <- predict(eff_gam, cpue_na) #predict cpue: kg/ha
  
  correction$eff[cpue_na_loc] = na_eff
  correction$ddc_cpue_kg_ha[cpue_na_loc] =  correction$cpue_kg_ha[cpue_na_loc]/na_eff
  correction$ddc_cpue_num_ha[cpue_na_loc] = correction$cpue_num_ha[cpue_na_loc]/na_eff
  
  return(correction)
}
