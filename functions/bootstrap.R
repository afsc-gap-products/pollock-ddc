#' This function is an updated (mostly snytax) version of the boostrapping 
#' methods for creating the variance-covariance matrix for the design-based
#' application of the density-dependent correction. Also included is the code
#' that Caitlin Allen-Akselrud commented out, which includes her attempt at 
#' creating the var-covar matrix. 
#' 
#' The function returns two dataframes that are saved in the design-based data
#' produced in the main code.
#' 
# Author: Caitlin I. Allen Akselrud
# Updated by: Sophia N. Wassermann
# Contact: sophia.wassermann@noaa.gov
# Date created: 2024.12.06

# convert new 4 -----------------------------------------------------------
# this is the conversion of stan's 'new4.R' which sets up the bootstrapping procedure to generate a v-cov matrix later
# CORRECTION: Caitlin's method did not produce the same amount of variance in the v-cov matrix, 
#   so that code is commented out and stan's code is implemented directly

##### caitlin's section ####

# # resample the data for each year and strata
# resample <- resamp_data(corr_cpue = ddc_table)
# 
# # check:
# # mean(resample$ddc_cpue_kg_ha)
# # mean(ddc_table$ddc_cpue_kg_ha)
# # mean(resample$ddc_cpue_num_ha, na.rm = T)
# # mean(ddc_table$ddc_cpue_num_ha, na.rm = T)
# 
# # cpue per station for one resample
# resample_pop <- resamp_pop(resamp = resample, strata = strata_metadata) %>% 
#   dplyr::filter(!is.na(year))
# 
# # CIA: YOU ARE HERE (last part of new4.R)
# param_sample <- sample(1000,1000,replace=T)
# 
# biomass <- matrix(nrow = length(param_sample), ncol = dim(resample_pop)[1])
# population <- matrix(nrow = length(param_sample), ncol = dim(resample_pop)[1])
# 
# 
# v_cov_units <- readline(prompt = "Which units do you want to use in the v-cov matrix (Enter: k for kg or t for th t): ")
# 
# k 
# 
# if(v_cov_units == 'k') 
# {
#   biomass <- foreach(i = 1:1000, .combine=rbind) %dopar% {
#     library(magrittr)
#     library(tidyverse)
#     
#     param_sample <- sample(1000,1,replace=T)
#     
#     # estimate ddc with random sample from parameter posterior distribution
#     backscatter_cpue_corr <- ddc_cpue_bs(equiv_bs_length, param_ests[param_sample,])
#     
#     # turn backscatter into cpue
#     ddc_table_vc <- ddc_fxn(backscatter_cpue_corr, cpue_info$p_cpue) #backscatter corr from EACH SAMPLE
#     
#     # randomly resample data from all years and all stations
#     resample <- resamp_data(ddc_table_vc) #FOR EACH SAMPLE
#     
#     # get total population est from each resample for each year (cpue*area)
#     resample_pop <- resamp_pop(resample, strata_metadata) %>%  #FOR EACH SAMPLE
#       dplyr::filter(!is.na(year))
#     
#     # save resampled biomass
#     resample_pop$total_ddc_kg #EACH ROW IS RESAMPLE, EACH COL IS YEARS
#   }
#   
#   write_csv(as.data.frame(biomass), here("output",dir_thisyr,paste0("bootstrap_biomass_kg_", data_type, ".csv")))
# } else if(v_cov_units == 't') {
#   biomass <- foreach(i = 1:1000, .combine=rbind) %dopar% {
#     library(magrittr)
#     library(tidyverse)
#     
#     param_sample <- sample(1000,1,replace=T)
#     
#     # estimate ddc with random sample from parameter posterior distribution
#     backscatter_cpue_corr <- ddc_cpue_bs(equiv_bs_length, param_ests[param_sample,])
#     
#     # turn backscatter into cpue
#     ddc_table_vc <- ddc_fxn(backscatter_cpue_corr, cpue_info$p_cpue) #backscatter corr from EACH SAMPLE
#     
#     # randomly resample data from all years and all stations
#     resample <- resamp_data(ddc_table_vc) #FOR EACH SAMPLE
#     
#     # get total population est from each resample for each year (cpue*area)
#     resample_pop <- resamp_pop(resample, strata_metadata) %>%  #FOR EACH SAMPLE
#       dplyr::filter(!is.na(year))
#     
#     # save resampled biomass
#     resample_pop$total_ddc_th_t #EACH ROW IS RESAMPLE, EACH COL IS YEARS
#   }
#   
#   write_csv(as.data.frame(biomass), here("output",dir_thisyr,paste0("bootstrap_biomass_th_t_", data_type, ".csv")))
# } else(print("Invalid selection: please select k for kg or t for th t"))

# v-cov check
# sqrt(diag(biomass))/mean(diag(biomass))

#### stan's bootstrapping method: ####
bootstrap_stan <- function() {
  if(data_type == "db") {
    corr_sa <-  function(aa,ab){
      # calculates corrected sa by tow (corrected backscatter by tow)
      by1 <- list(haul = aa$HAUL, subst = aa$SUBSTRATA, vessel = aa$VESSEL, year = aa$YEAR)
      sa_station <- aggregate(aa$bt_sa, by1, sum)
      str(sa_station)
      # hist(sa_station$x, 2000, xlim=c(0,5000))
      
      ii <- 1 
      sa_corr <- sa_station$x[ii] / (ab$asp + exp(-(ab$int + ab$slope * sa_station$x[ii])))
      corr_mean <- mean(sa_corr)
      corr_sd <- sd(sa_corr)
      
      for(ii in 2:length(sa_station$x)) {
        sa_corr <- sa_station$x[ii] / (ab$asp + exp(-(ab$int + ab$slope * sa_station$x[ii])))
        corr_mean[ii]=mean(sa_corr)
        corr_sd[ii] = sd(sa_corr) #might need na.rm = T
      }
      
      # this data contains old and new estimates of sa (uncorrected and corrected backscatter estimates returned)
      sa_sum <- cbind(sa_station,corr_mean,corr_sd)
      return(sa_sum)
    }
    
    corr_cpue <- function(sa,bb) {
      # calculates new cpue per tow (applies bs correction; eff = efficency; eff = ratio uncorr to corr, less than 1)
      eff <- sa$x / sa$corr_mean
      sa <- cbind(sa,eff)
      # str(sa)
      
      cc <- merge(bb,sa,all=T)
      
      new_kgha <- cc$CPUE_KGHA/cc$eff
      new_noha <- cc$CPUE_NOHA/cc$eff
      
      pollock_cpue_new <- cbind(cc[, c(1:11, 16)], new_kgha, new_noha)
      # str(pollock_cpue_new)
      pollock_cpue_new$new_kgha[pollock_cpue_new$CPUE_KGHA == 0] <- 0
      pollock_cpue_new$new_noha[pollock_cpue_new$CPUE_NOHA == 0] <- 0
      pollock_cpue_new$eff[pollock_cpue_new$CPUE_KGHA == 0] <- 1
      # write.csv(pollock_cpue_new, "pollock_cpue_new.csv")
      
      na.eff.cpue <- data.frame(xx = pollock_cpue_new[which(is.na(pollock_cpue_new$eff)), ]$CPUE_KGHA)
      na.eff.cpue.idx <- which(is.na(pollock_cpue_new$eff))
      
      xx <- pollock_cpue_new$CPUE_KGHA
      yy <- pollock_cpue_new$eff
      # plot(yy~xx)
      gam1 <- gam(yy ~ s(xx))
      # summary(gam1)
      # plot(gam1)
      
      na.eff <- predict(gam1, na.eff.cpue)
      pollock_cpue_new$eff[na.eff.cpue.idx] <- na.eff
      pollock_cpue_new$new_kgha[na.eff.cpue.idx] <- pollock_cpue_new$CPUE_KGHA[na.eff.cpue.idx] / na.eff
      pollock_cpue_new$new_noha[na.eff.cpue.idx] <- pollock_cpue_new$CPUE_NOHA[na.eff.cpue.idx] / na.eff
      # write.csv(pollock_cpue_new, "pollock_cpue_new.csv")
      # str(pollock_cpue_new)
      
      pollock_cpue_new2 <- pollock_cpue_new[, c(1:8, 13:14, 11)]
      # str(pollock_cpue_new2)
      colnames(pollock_cpue_new2)[c(1:3, 9, 10)] <- c("VESSEL", "YEAR", "HAUL", 
                                                      "CPUE_KGHA","CPUE_NOHA")
      # str(pollock_cpue_new2)
      
      return(pollock_cpue_new2) #returns the DDC cpue
    }
    
    ncpue_sample <- function(ncpue) {
      first.y <- min(ncpue$YEAR) #for each year
      last.y <- max(ncpue$YEAR)
      
      ii <- first.y
      ncpue.y <- ncpue[ncpue$YEAR == ii, ]  # select a specific year
      subs <- unique(ncpue.y$SUBSTRATA)  # for each strata; get strata from that year
      jj <- subs[1]
      ncpue.y.s <- ncpue.y[ncpue.y$SUBSTRATA == jj, ]  # select specific year and strata
      nn <- length(ncpue.y.s[, 1])  # number samples wanted
      samp <- sample.int(nn, size = nn, replace = T)  # randomly pick order, possibly with duplicates
      ncpue.y.s.n <- ncpue.y.s[samp,]  # arrange in random order
      for(jj in subs[-1]) {
        ncpue.y.s <- ncpue.y[ncpue.y$SUBSTRATA == jj, ] 
        nn <- length(ncpue.y.s[, 1])
        samp <- sample.int(nn, size = nn, replace = T)
        ncpue.y.s.n <- rbind(ncpue.y.s.n,ncpue.y.s[samp, ])
        # estimate of one realization of estimating cpue; with all substrata, 1st year
      }
      ncpue.y.n <- ncpue.y.s.n
      
      for(ii in (first.y + 1):last.y) {  # this is for all remaining years (with all substrata)
        ncpue.y <- ncpue[ncpue$YEAR == ii, ]
        subs <- unique(ncpue.y$SUBSTRATA)
        jj <- subs[1]
        ncpue.y.s <- ncpue.y[ncpue.y$SUBSTRATA == jj, ]
        nn <- length(ncpue.y.s[, 1])
        samp <- sample.int(nn, size = nn, replace = T)
        ncpue.y.s.n <- ncpue.y.s[samp, ]
        for(jj in subs[-1]){
          ncpue.y.s <- ncpue.y[ncpue.y$SUBSTRATA == jj, ]
          nn <- length(ncpue.y.s[, 1])
          samp <- sample.int(nn, size = nn, replace = T)
          ncpue.y.s.n <- rbind(ncpue.y.s.n, ncpue.y.s[samp, ])
        }
        ncpue.y.n <- rbind(ncpue.y.n,ncpue.y.s.n)
        
      }
      return(ncpue.y.n) #one realization of cpue est for all years and all substrata
    }
    
    corr_population <- function(aa,bb) {
      by1 <- list(aa$YEAR, aa$SUBSTRATA)
      cc <- aggregate(cbind(aa$CPUE_KGHA, aa$CPUE_NOHA), by1, mean, na.rm = T)  # get mean by year and strata
      colnames(cc) <- c("year", "substrata", "m_kgha", "m_noha")
      colnames(bb) <- c("substrata", "area")
      cc1 <- merge(cc, bb, all = T)  # add area back
      # cc1=cc1[cc1$substrata!=70,]  # filter out NBS
      # cc1=cc1[cc1$substrata!=81,]
      pop_by_strata <- cc1[, 3:4] * cc1[, 5]  # get number by multiplying cpue by area fished (cpue = num/ha; area = km^2; need to *100 to get km^2 -> ha-- done 3 lines below)
      by2 <- list(cc1$year)
      total_pop <- aggregate(pop_by_strata, by2, sum)  # sum by year
      total_pop1 <- cbind(year = total_pop[, 1], tons = total_pop[, 2] / 10, mlns = total_pop[, 3] / 10000) # years, convert to kg to tons, convernt num to millions
      return(total_pop1)
    }
    
    ab <- param_ests
    
    bb <- cpue_info$p_cpue %>%  #cpue
      rename(CRUISE = cruise, SUBSTRATA = stratum, 
             LATITUDE = start_latitude, LONGITUDE = start_longitude,
             CPUE_KGHA = cpue_kg_ha, CPUE_NOHA = cpue_num_ha, AREA_FISHED = area_fished_ha) %>% 
      mutate(SPECIES_CODE = 21740) %>% 
      select("vessel", "CRUISE", "year", "haul", "SUBSTRATA", "LATITUDE", "LONGITUDE", 
             "SPECIES_CODE", "CPUE_KGHA", "CPUE_NOHA", "AREA_FISHED") #%>% 
    # dplyr::filter(!SUBSTRATA %in% NBS_subarea)
    
    # bb_sub <- bb %>% dplyr::select("vessel", "CRUISE", "year", "haul", "SUBSTRATA", "LATITUDE", "LONGITUDE") %>% 
    #   clean_names(case = 'all_caps')
    
    aa <- equiv_bs_length %>%  #cpue length
      clean_names(case = 'all_caps') %>% 
      rename(cpue_mile = CPUE_MILE_SQ, ts = TARGET_STRENGTH, bt_sa = BACKSCATTER, 
             CRUISE = CRUISEJOIN, SUBSTRATA = STRATUM, CPUE_LENGTH = CPUE_LENGTH_NUM_HA) %>%
      # left_join(bb_sub) %>% 
      mutate(LATITUDE = NA, LONGITUDE = NA, 
             SPECIES_CODE = 21740) %>% 
      select("VESSEL", "YEAR", "CRUISE", "HAUL", "LATITUDE", "LONGITUDE", "SUBSTRATA", 
             "STATIONID", "SPECIES_CODE", "SEX", "LENGTH", "CPUE_LENGTH", 
             "cpue_mile",  "ts", "bt_sa") #%>% 
    # dplyr::filter(!SUBSTRATA %in% NBS_subarea)
    
    gg <- strata_metadata %>%  #strata
      rename(area_ha = area) %>% 
      select( "stratum", "area_ha") #%>% 
    # dplyr::filter(!stratum %in% NBS_subarea)
    
    ii <- sample(1000, 1)  # random number btw 1 and 1000
    sa_sum <- corr_sa(aa, ab[ii, ])  # estimates corr sa based on params from abc_exp
    ncpue <- corr_cpue(sa_sum, bb)  # estimates cpue based on param
    ncpue.s <- ncpue_sample(ncpue)  # gets ncpue sample
    estim <- corr_population(ncpue.s, gg)  # get corr pop based on random param selected
    pops <- estim[, 3]  # pop vector
    bioms <- estim[, 2]  # biomass vector
    
    for(ii in sample(1000, 999, replace = T)) {  # run 1000 times (999 more times), randomly sampling
      sa_sum <- corr_sa(aa, ab[ii, ])
      ncpue <- corr_cpue(sa_sum, bb)
      ncpue.s <- ncpue_sample(ncpue)
      estim <- corr_population(ncpue.s, gg)
      popsii <- estim[, 3]
      biomsii <- estim[, 2]
      pops <- rbind(pops, popsii)  # rows = resamples, cols = yrs
      bioms <- rbind(bioms, biomsii)
      print(length(pops[, 1]))
    }
    
    colnames(bioms) <- unique(aa$YEAR)
    bioms
    var_covar_stan <- cov(bioms)
    colnames(var_covar_stan) <- unique(aa$YEAR)
    rownames(var_covar_stan) <- unique(aa$YEAR)
    
    # check:
    biomass_mean <- bioms %>% 
      as_tibble() %>% 
      summarise_if(is.numeric, mean) 
    colnames(biomass_mean) <- unique(ddc_table$year)
    cv <- sqrt(diag(var_covar_stan)) / biomass_mean
    range(cv)  # range should be 5-20%
    if(max(range(cv)) > 0.2) {print("YOUR CV IS > 20% PLEASE CHECK")}
  }
  
  return(list(bioms = bioms, var_covar_stan = var_covar_stan))
}