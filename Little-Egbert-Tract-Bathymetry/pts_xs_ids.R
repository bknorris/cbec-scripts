# load libraries
library(raster)
library(rgdal)
library(sp)
library(maptools)
library(rgeos)

######################################################################################################################
#### Load data #######################################################################################################
######################################################################################################################
#
# 
# csv_Path<-"E:\\19-1019_Little_Egbert_Bathy\\Survey_Data\\Hypack"
# dir(csv_Path)
# csv_fn<-"Hypack_test.txt"
# 
# pts <- read.table(file.path(csv_Path, csv_fn), header=FALSE, as.is=TRUE)
# 
# # time difference which triggers a new xs
# time_threshold<-20 #seconds
#
######################################################################################################
# # Use of function...
# 
# pts<-pts_xs_ids(pts, ext_dist)
# 
#####################################################################################################
#function

pts_xs_ids<-function(pts, time_threshold){
  
  # remove points where z_Transducer or Depth is NA
  pts<-pts[which(!is.na(pts$z_transducer)),]
  pts<-pts[which(!is.na(pts$Depth)),]
  
  # add in time - uses current computer time
  pts$Time<-as.POSIXct(paste(Sys.Date(), pts$Time), tz="GMT")
  
  # find time difference in each survey point
  xs_diff<-diff(pts$Time)
  
  # find the index of points which are greater than the time threshold
  idx<-which(xs_diff>time_threshold)
  
  # add the new ids!
  pts$id<-1
  for(i in 1:length(idx)){
    if(i!=length(idx)){
      pts$id[(idx[i]+1):idx[(i+1)]]<-i+1
    } else {
      pts$id[(idx[i]+1):nrow(pts)]<-i+1
    }# close if statement
  }# close i loop
  
  # convert to spatial feature
  coordinates(pts)<- ~ x + y
  
  # plot xs data to make sure time threshold works (after 7 or so xs the colors may repeat)
  plot(pts, pch=20, cex=2, axes=TRUE, col=pts$id)
  
  return(pts)
  
} # close function
