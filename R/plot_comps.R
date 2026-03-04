# Plotting output from density-dependent correction

library(ggplot2)
library(here)
library(dplyr)

if (!requireNamespace("ggsidekick", quietly = TRUE)) {
  devtools::install_github("seananderson/ggsidekick")
}
library(ggsidekick)
theme_set(theme_sleek())

# current_year <- as.numeric(format(Sys.Date(), "%Y"))
current_year <- 2025
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

# Compare ALK files for model-based comps -------------------------------------
old_alk <- read.csv(here("output", "archive", "2024_mb_data_2022_strata", "VAST_ddc_alk_2024.csv"))
new_alk <- read.csv(here("output", "2026-03-02_mb", "VAST_ddc_alk_2026.csv")) 

# Check if new ages have been added
old_counts <- old_alk %>%
  group_by(Year, Age) %>%
  summarize(count = n())
new_counts <- new_alk %>%
  group_by(Year, Age) %>%
  summarize(count = n()) %>%
  filter(Year %in% old_counts$Year)

counts_diff <- bind_cols(Year = new_counts$Year, 
                         Age = new_counts$Age, 
                         count = (new_counts$count - old_counts$count)) %>%
  filter(count > 0)

# Compare total year/age across the datasets
old_comps <- old_alk %>%
  group_by(Year, Age) %>%
  summarize(comp = sum(CPUE_num))

new_comps <- new_alk %>%
  group_by(Year, Age) %>%
  summarize(comp = sum(CPUE_num)) %>%
  filter(Year %in% old_comps$Year)

compare_alk <- bind_cols(Year = new_comps$Year, 
                         Age = new_comps$Age, 
                         comp = (new_comps$comp - old_comps$comp)) %>%
  mutate(sign = case_when(comp >= 0 ~ "positive",
                          comp < 0 ~ "negative"))

# Plot both regions together and without 2020
ggplot(compare_alk, aes(x = Age, y = comp, fill = sign)) +
  geom_bar(stat = "identity", show.legend = FALSE) +
  scale_fill_manual(values = c("cornflowerblue", "darkred")) +
  scale_x_discrete(breaks = c(1, 5, 10, 15)) +
  ylab(paste0("difference between ", 2026, " and ", 2024, " DDC ALK")) +
  theme(legend.title = element_blank()) +
  facet_wrap(~ Year, ncol = 6, dir = "v") 

ggsave(filename = here("figures", "alk_diff.png"),
       width=200, height=200, units="mm", dpi=300)
