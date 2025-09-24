# This function maps data
# Author: Caitlin I. Allen Akselrud
# Contact: caitlin.allen_akselrud@noaa.gov
# Date created: 2022.05.11
# Date updated: 2022.05.11

map_data <- function(input, stn_color = "red", stn_title = "Station map", stn_filename = "stn_map", stn_path = here("figures"))
{
  library(ggmap)
  lat <- c(min(input$start_latitude), max(input$start_latitude))
  long <- c(min(input$start_longitude), max(input$start_longitude))
  bbox <- make_bbox(long,lat,f=0.05)
  b <- get_map(bbox,maptype="toner-lite",source="stamen")
  
  p6 <- ggmap(b) + geom_point(data = input, 
                              aes(start_longitude,start_latitude),size=.2,alpha=0.7, color = stn_color) +
    labs(x = "Longitude", y = "Latitude",
         title=stn_title) +
    facet_wrap(~ year)
  ggsave(p6, filename = paste0(stn_filename, Sys.Date(), ".png"), path = stn_path, height = 10, width = 10)
  
}
