# Plotting the design-based index from the density-dependent correction code

library(ggplot2)
library(here)
library(dplyr)

if (!requireNamespace("ggsidekick", quietly = TRUE)) {
  pak::pkg_install("seananderson/ggsidekick")
}
library(ggsidekick)
theme_set(theme_sleek())

current_year <- as.numeric(format(Sys.Date(), "%Y"))

# New index 
index <- read.csv(here(
  "output",
  "2026-09-01_db",
  paste0("biomass_densdep_corrected_", current_year, ".csv")
)) %>%
  select(year, biomass_MT_ha, low95_bm, high95_bm) %>%
  group_by(year) %>%
  summarize(
    biomass_MT_ha = sum(biomass_MT_ha) / 1e6,
    low95_bm = sum(low95_bm) / 1e6,
    high95_bm = sum(high95_bm) / 1e6
  ) 

ggplot() +
  geom_pointrange(
    data = index, 
    aes(x = year, y = biomass_MT_ha, ymin = low95_bm, ymax = high95_bm)
  ) +
  geom_line(
    data = index, 
    aes(x = year, y = biomass_MT_ha)
  ) +
  xlab("") + ylab("Biomass (Mt)")

# Read in old index and compare
old_index <- read.csv(here(
  "output",
  "2026-03-19_db",
  paste0("biomass_densdep_corrected_", current_year - 1, ".csv")
)) %>%
  select(year, biomass_MT_ha, low95_bm, high95_bm) %>%
  group_by(year) %>%
  summarize(
    biomass_MT_ha = sum(biomass_MT_ha) / 1e6,
    low95_bm = sum(low95_bm) / 1e6,
    high95_bm = sum(high95_bm) / 1e6
  ) %>%
  mutate(index_year = as.character(current_year - 1))

both_indices <- bind_rows(
  index %>% mutate(index_year = as.character(current_year)),
  old_index
) %>%
  ggplot(.) +
  geom_pointrange(
    aes(x = year, 
      y = biomass_MT_ha, 
      ymin = low95_bm, ymax = high95_bm,
      color = index_year
    )) +
  geom_line(aes(x = year, y = biomass_MT_ha, color = index_year)) +
  xlab("") + ylab("Biomass (Mt)")
both_indices
