## Script to import VAWT 10 Hz data from text file ############################
# 
# Raw data from an ultrasonic anemometer (81000, Young, USA) at 10 Hz
# Author: Yusri Yusup, PhD
# Date created: 2015-07-27
#
################################################################################

#### 1. Preliminaries ####




#### 2. Import data ####

# Listing all the files in the data folder
path <- list.files('data/rawdata')
# To remove the '.TXT' at the end
path1<-as.data.frame(strsplit(path,split='.TXT'))

# Import 10 Hz data
rawdata <- read.delim('data/rawdata/D3O_A1.TXT',sep='',header=FALSE)

# Creates a column for presence of VAWT, TRUE or FALSE
vawt <- rep(FALSE,nrow(rawdata)) # Need to get from text file name

# Creates a column for speed setting, 1 = low, 3 = high
speed <- rep(1,nrow(rawdata)) # Need to get from the text file name

# Creates a column for measurement position
pos <- rep('A1',nrow(rawdata)) # Need to get from the text file name

# Creates a column for measurement distance, upwind or downwind
dist <- rep('downwind',nrow(rawdata)) # Need to get from the text file name

# Creates a dataframe of the above two parameters
df <- data.frame(vawt,speed,dist,pos)

# Combines into one dataframe
df <- cbind(df,rawdata)

# Name the columns
names(df) <- c('vawt','speed','dist','pos','u','v','w','T')


# Temporary variable cleanup
rm(pos,speed,vawt,rawdata,dist)
