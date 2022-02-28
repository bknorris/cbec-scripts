# install packages - this only needs to be done once
#install.packages(c("raster", "rgdal", "sp", "maptools", "rgeos"))

# load libraries (do this even after installation)
library(raster)
library(rgdal)
library(sp)
library(maptools)
library(rgeos)

####################################################################################################
#### Bathymetery Toolbox Example Workflow - XS data extraction and comparison  #####################
####################################################################################################
# Set your working directory (set to root folder "RAS_Tools")
# Add a second slash "\\" to the path 
setwd("C:\\Work\\19-1019_Little_Egbert")

#look at what is in the wd
dir()

####################################################################################################
# import functions

# folder with R scripts
codePath<-"C:\\Work\\19-1019_Little_Egbert\\Program_R"

# import custom funcions stored in other scripts
source(file.path(codePath, "pts_xs_ids.R"))
source(file.path(codePath, "survey_pts_To_xs_lines.R"))
source(file.path(codePath, "extract_XS_dem.R"))
source(file.path(codePath, "merge_XS_dem_bathy.R"))

####################################################################################################
# import your GIS data

# Read in base raster data
# ensure survey and DEM are already in same projection!!!
DEM_Path<-".\\GIS\\Raster"
dir(DEM_Path, pattern=".tif$") # print tif names to the console
#DEM_Name<-"liberty_1m_bathy_SP2_ft_V2_Zft.tif" # with file extension USGS version
DEM_Name<-"Prospect_Island_v5_SP2_FT_Z_FT.tif" #cbec ver

dem<-raster(file.path(DEM_Path, DEM_Name))
#plot(dem, col=rev(rainbow(100)))

# check projection of DEM (must match point data!)
crs(dem)

#####################################################################################################
# Read in survey point data
#####################################################################################################

# path to csv files 
csv_Path<-".\\Survey_Data\\" #may want a new folder w/ date
dir(csv_Path) # check what files are in the directory
# survey data file name (Columns = "Time", "x", "y", "Z_transducer", "Depth", "Z_bed")
csv_fn<-"061919-_000_1424 from hypack.txt" # set to the file name

# read in data
survey_pts <- read.table(file.path(csv_Path, csv_fn), header=FALSE, as.is=TRUE, fill=NA)

# Add file names
names(survey_pts)<-c("Time", "x", "y", "z_transducer", "Depth", "z_bed") # file must have z_bed for functions to work!

# Time difference between when points are collected which triggers a new xs id
time_threshold<-20 #seconds - may need to change if survey split up too much

# break up points to line ids based on time
survey_pts<-pts_xs_ids(survey_pts, time_threshold)

# set projection of pts layer to be the same as the DEM
crs(survey_pts)<-crs(dem) 

#####################################################################################################
# Create lines from the survey points
#####################################################################################################

#line extention distance
ext_dist<-100 # map units, may need to extend to LIDAR on bank

# create xs_lines from points
xs_lines<-survey_pts_To_xs_lines(survey_pts, ext_dist)

# set projection of xs_lines layer to be the same as the DEM
crs(xs_lines)<-crs(dem) 

# Write xs lines to a shapefile.  
# gis.outDir<-"E:\\19-1019_Little_Egbert_Bathy\\GIS\\Shapefiles"
# fn.out<-"Survey_XS_Lines" # no file extension
# writeOGR(xs_lines, dsn=gis.outDir, layer=fn.out, driver="ESRI Shapefile", overwrite_layer=TRUE)

#####################################################################################################
# extract and merge bathy data and xs data from DEM along cross sections 

# set sample distance for raster extraction
samp.d<-3 #uses linear referencing eg, feet

# extract survey data
survey<-merge_XS_dem_bathy(xs_lines, dem, samp.d, survey_pts)

#####################################################################################################

# extract dem data
dem_xs<-extract_XS_dem(xs_lines, dem, samp.d, survey_pts)

###################################################################################################
###################################################################################################
###################################################################################################
# PLOT DATA!

# create plot folder
plotDir<-".\\Plots"
dir.create(plotDir)

# create a subfolder in the plot dir. 
plotDir<-paste0(".\\Plots\\", substr(csv_fn, 1, nchar(csv_fn)-4))
dir.create(plotDir)

xsNames<-as.character(unique(survey$xsNames))

for(index in 1:length(xsNames)){
  
  name<-xsNames[index]

  print(paste0("Now plotting cross section ", name))
  
  fn<-paste0(plotDir,"\\", substr(csv_fn, 1, nchar(csv_fn)-4), "_XS_Plot_", gsub(" ","_",name), ".jpeg")
  jpeg(filename = fn, width = 3000, height = 1800, quality = 100)
  
    par(mar=c(10, 12, 4, 4), mgp = c(8,3,0))
    
    yl<-c(min(dem_xs[dem_xs$id==index, "elev"], survey[survey$id==index & survey$source=="Bathy", "elev"], na.rm=TRUE),
          max(dem_xs[dem_xs$id==index, "elev"], survey[survey$id==index & survey$source=="Bathy", "elev"], na.rm=TRUE))
    
    plot(dem_xs[dem_xs$id==index,c("station","elev")] , type="l", col="black", lty=1, lwd=5,
         panel.first=grid(), main=name, ylim=yl,
         xlab="Station (ft)", ylab="Elevation (ft)",
         cex.lab=4, cex.axis=4, cex.main=4)
    
    lines(survey[survey$id==index & survey$source=="Bathy", c("station","elev")] , col="red", lty=3, lwd=5)
    points(survey[survey$id==index & survey$source=="Bathy", c("station","elev")] , col="blue", cex=3, pch=20)
    
    legend("bottomleft", c("DEM", "Survey"), lty=c(1,3), 
           lwd=rep(5,2), col=c("black", "red"), cex=4)
  
  dev.off()

}# close i loop

################################################################################################################

# plot map data

fn<-paste0(plotDir,"\\", substr(csv_fn, 1, nchar(csv_fn)-4), "_Location_Map.jpeg")
jpeg(filename = fn, width = 3000, height = 1800, quality = 100)

par(mar=c(10, 12, 4, 4), mgp = c(8,3,0))

  plot(xs_lines, lwd=3, panel.first=grid(), axes=TRUE, col="darkgreen",
       main=substr(csv_fn, 1, nchar(csv_fn)-4),
       cex.lab=4, cex.axis=4, cex.main=4)
  
  points(survey_pts, pch=20, col=survey_pts$id, cex=3)  
  
  # find the start point of each line
  startPts<-data.frame()
  for(i in 1:length(xs_lines)){
    temp<-data.frame(xs_lines@lines[[i]]@Lines[[1]]@coords)
    #identify the name column
    temp$id<-i
    startPts<-rbind(startPts,temp[1,])
  }
  names(startPts)<-c("X","Y","ID")
  coordinates(startPts)<- ~ X + Y
  
  text(startPts, startPts$ID, halo=TRUE, hc='blue', col='white', hw=1, cex=4)

dev.off()

###################################################################################################

