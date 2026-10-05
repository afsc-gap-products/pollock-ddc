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
  pak::pkg_install("afsc-gap-products/gapindex", build_vignettes = TRUE)
}
library(gapindex)

# Set ggplot theme
if (!requireNamespace("ggsidekick", quietly = TRUE)) {
  pak::pkg_install("seananderson/ggsidekick")
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
this_year <- as.integer(format(Sys.Date(), "%Y"))
odbcGetInfo(channel)  # check connection

# First, pull data from the standard EBS stations
species_code <- c(21740, 21741)
ebs_standard_data <- get_data(
  year_set = 1982:this_year,
  survey_set = "EBS",
  spp_codes = species_code,
  pull_lengths = FALSE, 
  haul_type = 3, 
  abundance_haul = "Y",
  channel = channel,
  remove_na_strata = TRUE
)

#' Next, pull data from hauls that are not included in the design-based index
#' production (abundance_haul == "N") but are included in VAST. By default, the 
#' gapindex::get_data() function will filter out hauls with negative performance 
#' codes (i.e., poor-performing hauls).
ebs_other_data <- get_data(
  year_set = c(1994, 2001, 2005, 2006),
  survey_set = "EBS",
  spp_codes = species_code,
  pull_lengths = FALSE, 
  haul_type = 3, 
  abundance_haul = "N",
  channel = channel, 
  remove_na_strata = TRUE
)

# Combine the EBS standard and EBS other data into one list. 
ebs_data <- list(
  survey = ebs_standard_data$survey,
  survey_design = ebs_standard_data$survey_design,
  #' Some cruises are shared between the standard and other EBS cruises, so the 
  #' unique() wrapper is there to remove duplicate cruise records. 
  cruise = unique(rbind(ebs_standard_data$cruise, ebs_other_data$cruise)),
  haul = rbind(ebs_standard_data$haul, ebs_other_data$haul),
  catch = rbind(ebs_standard_data$catch, ebs_other_data$catch),
  species = ebs_standard_data$species,
  strata = ebs_standard_data$strata
)

# Calculate CPUE and export 
cpue_out <- calc_cpue(gapdata = ebs_data)

ebs_cpue <- cpue_out %>%
  select("YEAR", "LATITUDE_DD_START",
         "LONGITUDE_DD_START", "CPUE_KGKM2") %>%
  transmute(
    Lat = LATITUDE_DD_START,
    Lon = LONGITUDE_DD_START,
    Year = as.integer(YEAR),
    CPUE =  CPUE_KGKM2,
    data = "base"
  ) 

# write.csv(ebs_cpue, file = here("data", "uncorrected_cpue.csv"), row.names = FALSE)

# Read in density-corrected cpue and combine for comparisons ------------------
ddc_cpue <- read.csv(here("output", "2026-09-25_db", "VAST_ddc_EBSonly_2026.csv")) %>%
  transmute(
    Lat = start_latitude,
    Lon = start_longitude,
    Year = year,
    CPUE = ddc_cpue_kg_ha * 100,  # convert to kg/km2
    data = "DDC"
  ) 

# Plot both, separately -------------------------------------------------------
cpue_sf <- st_as_sf(bind_rows(ebs_cpue, ddc_cpue), coords = c("Lon", "Lat"), crs = 4326) %>%
  mutate(logCPUE = log(CPUE)) %>%
  filter(Year >= (this_year - 4))
world <- ne_countries(scale = "medium", returnclass = "sf")
sf_use_s2(FALSE)  # turn off spherical geometry
ggplot(cpue_sf) +
  geom_sf(data = world) +
  geom_sf(data = cpue_sf, aes(color = logCPUE), size = 1) +
  scale_color_viridis(na.value = NA) +
  coord_sf(xlim = c(-180, -157), ylim = c(53.8, 63.2), expand = FALSE) +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank()) +
  labs(color = "log(CPUE (kg/km2))") +
  facet_grid(data ~ Year)

ggsave(filename = here("output", "comparison", "comparison_maps.png"), width = 10, height = 3, units = "in", dpi = 300)

# Plot difference between them ------------------------------------------------
diff_cpue <- inner_join(ddc_cpue, ebs_cpue, join_by(Lat == Lat, Lon == Lon, Year == Year)) %>%
  group_by(Lat, Lon, Year) %>%
  summarize(CPUE = ((CPUE.x - CPUE.y) / CPUE.y) * 100) %>%
  filter(Year >= (this_year - 4)) 
diff_cpue <- st_as_sf(diff_cpue, coords = c("Lon", "Lat"), crs = 4326)

ggplot(diff_cpue) +
  geom_sf(data = world) +
  geom_sf(data = diff_cpue, aes(color = CPUE), size = 1) +
  scale_color_viridis(na.value = NA) +
  coord_sf(xlim = c(-180, -157), ylim = c(53.8, 63.2), expand = FALSE) +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank()) +
  # theme(legend.position = "bottom") +
  labs(color = "Percent Change") +
  facet_wrap(~ Year, ncol = 6)

ggsave(filename = here("output", "comparison", "difference_maps.png"), width = 10, height = 2, units = "in", dpi = 300)

# Compare indices -------------------------------------------------------------
# Read in DDC biomass estimate
ddc_db <- read.csv(here("output", "2026-09-25_db", "biomass_densdep_corrected_2026.csv")) %>%
  group_by(year) %>%
  summarize(
    biomass = sum((biomass_MT_ha)),
    # The variance is HUGE, but these are the same as the biomass estimate....
    low95_bm = sum(low95_bm),
    high95_bm = sum(high95_bm)
  ) %>% 
  mutate(method = "DDC")

# Read in uncorrected biomass estimate
uncorrected <- read.csv(here("output", "2026-09-25_db", "biomass_uncorrected_2026.csv")) %>%
  group_by(year) %>%
  summarize(
    biomass = sum((biomass_MT_ha)),
    low95_bm = sum(low95_bm),
    high95_bm = sum(high95_bm)
  ) %>% 
  mutate(method = "uncorrected")

ggplot(bind_rows(ddc_db, uncorrected)) +
  geom_line(aes(x = year, y = biomass, color = method), alpha = 0.4) +
  geom_point(aes(x = year, y = biomass, color = method, shape = method)) +
  scale_color_manual(values = c("darkblue", "lightslateblue")) 

ggsave(filename = here("output", "comparison", "biomass_ddc.png"), width = 8, height = 5, units = "in", dpi = 300)
