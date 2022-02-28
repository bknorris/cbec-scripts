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
# #read in bathy point data
# #must have id column
# pts<-readOGR(dsn=getwd(),layer="Feather_Survey_pts")
# 
# #line extention
# ext_dist<-200 # map units
#
######################################################################################################
# # Use of function...
# 
# xs_lines<-survey_pts_To_xs_lines(pts, ext_dist)
# 
#####################################################################################################
#function

survey_pts_To_xs_lines<-function(pts, ext_dist){
  
  # check if data has id column
  if(!"id" %in% colnames(pts@data)){
    stop("Points file needs \"id\" column")
  } # close if
  
  for(i in 1:length(unique(pts$id))){
    
    print(paste0("Extracting Cross section ", i))
    
    temp.pts<-pts[pts$id==i,]
    #plot(temp.pts, pch=20, axes=TRUE, col="blue", panel.first=grid(), main=paste0("Cross section ", i))
    
    xy<-data.frame(coordinates(temp.pts))

    new.x<-data.frame(x=xy$x[c(1,nrow(xy))])
    
    # define linear model
    model<-lm(y ~ x, data=xy)
    
    # slope
    m<-model$coefficients[2]
    
    # find change in x which matches line extension distance (ext_dist)
    dx<-ext_dist / sqrt(m^2 + 1)
    
    # find the new x end points
    # x coord
    if(new.x[1,1]>new.x[2,1]){
      new.x[1,1]<-new.x[1,1] + dx
      new.x[2,1]<-new.x[2,1] - dx
    } else {
      new.x[1,1]<-new.x[1,1] - dx
      new.x[2,1]<-new.x[2,1] + dx
    } # close if 
    
    # predict model
    pred.y<-predict(lm(y ~ x, data=xy), new.x)

    temp.line<-data.frame(x=new.x, y=pred.y)
    
    #plot(temp.line, lwd=2, type="l", panel.first=grid(), main=paste0("Cross section ", i), axes=TRUE, )
    #points(temp.pts, pch=20, col="blue")

    temp.line<-SpatialLines(list(Lines(list(Line(temp.line)), 1)))
    temp.line<-SpatialLinesDataFrame(temp.line, data.frame(id=i))
    
    if(i==1){
      xs_lines <- temp.line
    } else {
      xs_lines<-rbind(xs_lines, temp.line)
    }
    
  }# close i loop
    
  plot(xs_lines, lwd=2, panel.first=grid(), axes=TRUE)
  points(pts, pch=20, col=pts$id)  

  return(xs_lines)

}# close function
  
########################################################################################################
    
    