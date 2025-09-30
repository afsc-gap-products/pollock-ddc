# Plotting output from density-dependent correction

library(ggplot2)
library(here)
library(dplyr)

library(ggsidekick)
theme_set(theme_sleek())

current_year <- as.numeric(format(Sys.Date(), "%Y"))
# Strata metadata year; 2022 is the latest update (use for current assessments)
strat_meta_year <- 2022
# Set output - model- or design-based
data_type <- "db"

alk_summary <- read.csv(here("output", 
                             paste0(current_year,"_", data_type, "_data_", strat_meta_year, "_strata"),
                             paste0("age_length_key_SUMMARY_densdep_corrected_", current_year, ".csv")))

alk <- na.omit(alk_summary) %>%
  group_by(year, age) %>%
  summarize(pop_at_age = sum(pop_at_age_total)) %>%
  ungroup() %>%
  group_by(year) %>%
  mutate(total = sum(pop_at_age)) %>%
  ungroup() %>%
  group_by(year, age) %>%
  mutate(proportion = pop_at_age / total) %>%
  ungroup()

ggplot(alk, aes(x = age, y = proportion)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~ year)
