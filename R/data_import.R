## Script to import VAWT 10 Hz data from text file ############################
# 
# Raw data from an ultrasonic anemometer (81000, Young, USA) at 10 Hz
# Author: Yusri Yusup, PhD
# Date created: 2015-07-27
# Date modified: 2015-10-08
#
################################################################################

#### 1. Preliminaries #### To call up the function
source('R/vawt_importdata.R')

#### 2. Import data ####

# Listing all the files in the data folder
filename <- list.files('data/rawdata')
path <- paste('data/rawdata/',filename, sep='')

# Import the first file first
df_final <- vawt_importdata() # Default is the first file

# Import 10 Hz data
for (i in 2:length(filename)){
  df <- vawt_importdata(path[i],filename[i])
  df_final <- rbind(df_final,df)
}
