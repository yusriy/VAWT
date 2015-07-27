## Script to import VAWT 10 Hz data from text file ############################
# 
# Raw data from an ultrasonic anemometer (81000, Young, USA) at 10 Hz
# Author: Yusri Yusup, PhD
# Date created: 2015-07-27
#
################################################################################

#### 1. Preliminaries ####



#### 2. Import data ####


# Creates a column for presence of VAWT, TRUE or FALSE
vawt <- rep(TRUE,nrow(rawdata))
# Creates a column for measurement position
pos <- rep('A1',nrow(rawdata))
# Creates a dataframe of the above two parameters
df <- data.frame(vawt,pos)
# Import 10 Hz data
rawdata <- read.delim('rawdata/DW_3_O/D3O_A1.TXT',sep='',header=FALSE)
# Combines into one dataframe
df <- cbind(df,rawdata)
# Name the columns
names(df) <- c('vawt','pos','u','v','w','T')
