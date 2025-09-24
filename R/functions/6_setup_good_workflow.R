# Good workflow: an easy way to set up consistent file folders
# Created by: Caitlin Allen Akselrud
# Contact: caitlin.allen_akselrud@noaa.gov
# Created: 2020-12-18
# Modified: 2021-01-18


# here are the names of the file folders you want to create:
dirs = c("code", "data", "documentation", "figures", "functions", "output")

# for each file folder, do the following:
for(i in 1:length(dirs)){
  # if the file folder already exists, do nothing
  if(dir.exists(dirs[i])==FALSE){
    # if the file folder does not exist, create it
    dir.create(dirs[i])
  }
}

# IF you have set up an R Project, you can run these lines of code
#  and it will create these file folders where your project already is

# WOW! That's nifty.