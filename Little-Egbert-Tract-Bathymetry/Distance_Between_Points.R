# load libraries (do this even after installation)
library(raster)
library(rgdal)
library(sp)
library(maptools)
library(rgeos)

data<-data.frame()
for(i in unique(survey$id)){
  
  tmp.pts<-survey[survey$id==i & survey$source=="Bathy",]
  
  for(j in 2:nrow(tmp.pts)){
    tmp.dist<-tmp.pts[i,"station"] - tmp.pts[i-1,"station"]
    
    data<-rbind(data, tmp.dist)
    
  }
  
}


hist(data[,1])

summary(data[,1])
