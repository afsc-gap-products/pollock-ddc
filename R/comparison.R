#' Script to demonstrate the effect of the density-dependent correction on
#' pollock CPUE. The uncorrected uncorrected CPUE is downloaded using the 
#' standard gapindex workflow used for other groundfish stocks in Alaska. For 
#' the moment, this comparison is only for the EBS.

library(here)
library(dplyr)
library(RODBC)
library(ggplot2)
library(viridis)
library(sf)
library(rnaturalearth)

if (!requireNamespace("gapindex", quietly = TRUE)) {
  install_github("afsc-gap-products/gapindex", build_vignettes = TRUE)
}
library(gapindex)

# Set ggplot theme
if (!requireNamespace("ggsidekick", quietly = TRUE)) {
  devtools::install_github("seananderson/ggsidekick")
}
library(ggsidekick)
theme_set(theme_sleek())

# Pull uncorrected pollock CPUE (biomass) -------------------------------------
# Connect to Oracle
if (file.exists("Z:/Projects/ConnectToOracle.R")) {
  source("Z:/Projects/ConnectToOracle.R")
} else {
  # For those without a ConnectToOracle file
  channel <- odbcConnect(dsn = "AFSC", 
                         uid = rstudioapi::showPrompt(title = "Username", 
                                                      message = "Oracle Username", 
                                                      default = ""), 
                         pwd = rstudioapi::askForPassword("Enter Password"),
                         believeNRows = FALSE)
}

# Plot each -------------------------------------------------------------------

odbcGetInfo(channel)  # check connection

# First, pull data from the standard EBS stations
species_code <- c(21740, 21741)
ebs_standard_data <- get_data(year_set = 1982:as.integer(format(Sys.Date(), "%Y")),
                              survey_set = "EBS",
                              spp_codes = species_code,
                              pull_lengths = FALSE, 
                              haul_type = 3, 
                              abundance_haul = "Y",
                              channel = channel,
                              remove_na_strata = TRUE)

#' Next, pull data from hauls that are not included in the design-based index
#' production (abundance_haul == "N") but are included in VAST. By default, the 
#' gapindex::get_data() function will filter out hauls with negative performance 
#' codes (i.e., poor-performing hauls).
ebs_other_data <- get_data(year_set = c(1994, 2001, 2005, 2006),
                           survey_set = "EBS",
                           spp_codes = species_code,
                           pull_lengths = FALSE, 
                           haul_type = 3, 
                           abundance_haul = "N",
                           channel = channel, 
                           remove_na_strata = TRUE)

# Combine the EBS standard and EBS other data into one list. 
ebs_data <- list(survey = ebs_standard_data$survey,
                 survey_design = ebs_standard_data$survey_design,
                 #' Some cruises are shared between the standard and other EBS cruises, so the 
                 #' unique() wrapper is there to remove duplicate cruise records. 
                 cruise = unique(rbind(ebs_standard_data$cruise, ebs_other_data$cruise)),
                 haul = rbind(ebs_standard_data$haul, ebs_other_data$haul),
                 catch = rbind(ebs_standard_data$catch, ebs_other_data$catch),
                 species = ebs_standard_data$species,
                 strata = ebs_standard_data$strata)

# Calculate CPUE and export 
ebs_cpue <- calc_cpue(gapdata = ebs_data) %>%
  select("YEAR", "LATITUDE_DD_START",
         "LONGITUDE_DD_START", "CPUE_KGKM2") %>%
  transmute(Lat = LATITUDE_DD_START,
            Lon = LONGITUDE_DD_START,
            Year = as.integer(YEAR),
            CPUE =  CPUE_KGKM2,
            data = "base") 

# write.csv(ebs_cpue, file = here("data", "uncorrected_cpue.csv"), row.names = FALSE)

# Read in density-corrected cpue and combine for comparisons ------------------
ddc_cpue <- read.csv(here("output", "2025-11-25_mb", "VAST_ddc_EBSonly_2025.csv")) %>%
  transmute(Lat = start_latitude,
            Lon = start_longitude,
            Year = year,
            CPUE = ddc_cpue_kg_ha * 100,  # convert to kg/km2
            data = "DDC") 

# Plot both, separately -------------------------------------------------------
cpue_sf <- st_as_sf(bind_rows(ebs_cpue, ddc_cpue), coords = c("Lon", "Lat"), crs = 4326) %>%
  mutate(logCPUE = log(CPUE)) %>%
  filter(Year >= 2021)
world <- ne_countries(scale = "medium", returnclass = "sf")
sf_use_s2(FALSE)  # turn off spherical geometry
ggplot(cpue_sf) +
  geom_sf(data = world) +
  geom_sf(data = cpue_sf, aes(color = logCPUE), shape = "square") +
  scale_color_viridis(na.value = NA, 
                      guide = guide_colourbar(title.position = "bottom",
                                              title.hjust = 0.5)) +
  coord_sf(xlim = c(-180, -157), ylim = c(53.8, 63.2), expand = FALSE) +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank()) +
  theme(legend.position = "bottom") +
  labs(color = "log(CPUE (kg/km2))") +
  facet_grid(data ~ Year)

# Plot difference between them ------------------------------------------------
diff_cpue <- inner_join(ddc_cpue, ebs_cpue, join_by(Lat == Lat, Lon == Lon, Year == Year)) %>%
  group_by(Lat, Lon, Year) %>%
  summarize(CPUE = (CPUE.x - CPUE.y) / CPUE.y) %>%
  filter(Year >= 2021) 
diff_cpue <- st_as_sf(diff_cpue, coords = c("Lon", "Lat"), crs = 4326)

ggplot(diff_cpue) +
  geom_sf(data = world) +
  geom_sf(data = diff_cpue, aes(color = CPUE), shape = "square") +
  scale_color_viridis(na.value = NA, 
                      limits = c(0, NA),
                      guide = guide_colourbar(title.position = "bottom", 
                                              title.hjust = 0.5)) +
  coord_sf(xlim = c(-180, -157), ylim = c(53.8, 63.2), expand = FALSE) +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank()) +
  theme(legend.position = "bottom") +
  labs(color = "Relative difference in CPUE") +
  facet_wrap(~ Year)
