# This function saves raw data from Oracle for the Bering Sea Pollock assessment
# Author: Caitlin I. Allen Akselrud
# Contact: caitlin.allen_akselrud@noaa.gov
# Date created: 2021.08.03
# Date updated: 2021.08.03

#

save_raw_data <- function(...)
  # hauls_survey, pollock_specimen, pollock_catch,
                          # pollock_length, ebs_shelf_cruise, metadata, 
                          # strata_metadata_raw)
{
  dots <- list(...)
  dots_names <- names(dots)
  length(dots)
  for(i in 1:length(dots))
  {
    write_csv(dots[[i]], here("data", paste0("raw_data_", dots_names[i], "_", today(), ".csv")))
  }
}