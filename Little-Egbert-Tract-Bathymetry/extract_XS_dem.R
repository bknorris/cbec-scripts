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
# # xs lines to snap data to and extract xs from
# xs<-readOGR(dsn=getwd(),layer="Feather_River_CVFED_XS")
# 
# # read in DEM data for over banks
# dem<-raster("Feather_River_2017_DEM.tif")
#
# # set sample distance for raster extraction
# samp.d<-3
#
# #read in bathy point data
# pts<-readOGR(dsn=getwd(),layer="Feather_Survey_pts")
#
#
######################################################################################################
# # Use of function...
# 
# xs_merge<-merge_XS_dem_bathy(xs_lines, dem, samp.d, pts)
# 
#####################################################################################################
#function

extract_XS_dem<-function(xs_lines, dem, samp.d, pts){
  
  # ensure data is already all in the same projection (this fixes minor differences in proj4 strings)
  if(projection(dem)!=projection(pts) | projection(xs_lines)!=projection(pts)){
    warning("Check input projection to ensure files match!")
  }
  
  # extract dem data at points current locations (before snapping to line)
  pts$dem_z<-extract(dem, pts)
  
  ######################################################################################################################
  #### Extract dem points along xs lines ###############################################################################
  ######################################################################################################################
  
  xs_table<-data.frame()
  
  for(index in 1:nrow(xs_lines)){
    
    print(paste0("Extracting Cross section ", xs_lines@data[index,"id"]))
    
    temp.line<-xs_lines[index,]
    lines(temp.line, col="red", lwd=2)
    
    # find points along the line that fit the sample distance
    d<-seq(1,floor(gLength(temp.line)),samp.d)
    samp.pts<-gInterpolate(temp.line,d)
    #plot(pts)
    
    ###########################################################################################################
    # Process XS data
    
    xs<-data.frame(x=samp.pts@coords[,1], y=samp.pts@coords[,2], station=d, elev=extract(dem, samp.pts), source="DEM")
    
    if(is.null(xs)) {
      print("DEM does not overlap with cross-section!")
      next
    }
    
    ###########################################################################################################
    # Process bathy data
    
    # subset point data
    temp.pts<-pts[pts$id==index,]
    
    if(nrow(temp.pts)==0) {
      print("Points do not overlap with cross-section!")
      next
    }
    
    #plot(temp.line)
    #points(temp.pts)
    
    bathy<-data.frame(x=temp.pts@coords[,1], y=temp.pts@coords[,2], station=gProject(temp.line, temp.pts), elev=temp.pts$dem_z, source="DEM_Bathy")
    
    ###############################################################################################
    
    # round to the nearest tens place
    dem.round.sta<-round(xs$station/10,0)*10
    pts.round.sta<-round(bathy$station/10,0)*10
    
    # remove points from dem xs that the bathy cover
    xs<-xs[! dem.round.sta %in% pts.round.sta,]
  
    # merge data and plot
    xs_pts<-rbind(xs, bathy)
    xs_pts<-xs_pts[order(xs_pts$station),]
    row.names(xs_pts)<-seq(1,nrow(xs_pts),1)
    
    plot(xs_pts[,c("station","elev")], type="l", main=paste0("Cross section ", index))
    points(xs_pts[,c("station","elev")], pch=20, col=ifelse(xs_pts$source=="DEM","dark green","blue"))
    
    xs_pts$id<-index
    
    #add xs name
    xs_pts$xsNames<-temp.line@data[,"id"]

    xs_table<-rbind(xs_table, xs_pts)
    
  }#close i loop
    
  return(xs_table)

}#close function