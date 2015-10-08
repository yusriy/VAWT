## Script to create function to import VAWT 10 Hz data from text file ############################
# 
# Raw data from an ultrasonic anemometer (81000, Young, USA) at 10 Hz
# Author: YinMun H'ng
# Date created: 2015-10-28
#
################################################################################

# Create the function to import the data from VAWT file
vawt_importdata <- function(path='data/rawdata/D1O_A1.TXT',filename='D1O_A1.TXT'){  
  
  # Import 10 Hz data
  rawdata <- read.delim(path,sep='',header=FALSE)
  
  # Creates a column for presence of VAWT, TRUE or FALSE
  # Need to get from text file name
  vawt <- if (substring(filename,3, 3) <= "O")
    rep(FALSE,nrow(rawdata)) else rep(TRUE,nrow(rawdata))
  
  # Creates a column for speed setting, 1 = low, 3 = high
  # Need to get from the text file name
  speed <- if (substring(filename,2, 2) <= "1")
    rep(1,nrow(rawdata)) else rep(3,nrow(rawdata))
  
  # Creates a column for measurement position
  # Need to get from the text file name
  pos <- rep(substring(filename,5, 6),nrow(rawdata)) 
  
  # Creates a column for measurement distance, upwind or downwind
  # Need to get from the text file name
  dist <- if (substring(filename,1, 1) <= "D")
    rep('downwind',nrow(rawdata)) else rep('upwind',nrow(rawdata))
  
  # Creates a dataframe of the above two parameters
  df <- data.frame(vawt,speed,dist,pos)
  
  # Combines into one dataframe
  df <- cbind(df,rawdata)
  
  # Name the columns
  names(df) <- c('vawt','speed','dist','pos','u','v','w','T')
  
  # Return the final data
  return(df)
  
}


  