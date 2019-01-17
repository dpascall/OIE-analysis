###code for reading in OIE data on a per sequence level and generating 
library(sp)
library(raster)
library(rgdal)
library(rgeos)
library(maptools)
library(spatstat)
library(sf)
library(lwgeom)
library(geosphere)
library(lubridate)

OIE<-read.csv("~/Desktop/OIEdatasetforcovariate.csv")[,c(4,5,7,8,10,11,12,14)]
OIE<-OIE[!OIE$Cases%in%0,]
OIE<-OIE[!is.na(OIE$Cases),]

##convert dates to readable format

##convert dates to readable format

OIE$Latitude<-as.numeric(OIE$Latitude)
OIE$Longitude<-as.numeric(OIE$Longitude)
OIE$OBStartDate<-as.character(OIE$OBStartDate)
OIE$OBEndDate<-as.character(OIE$OBEndDate)

OIEreduced<-OIE[OIE$OBEndDate=="",]
OIE<-OIE[OIE$OBEndDate!="",]

splitstart<-do.call("rbind",strsplit(OIE$OBStartDate,"/"))
splitstart[,3]<-paste("20",splitstart[,3],sep="")
OIE$OBStartDate<-paste(splitstart[,3],splitstart[,2],splitstart[,1],sep="-")

splitend<-do.call("rbind",strsplit(OIE$OBEndDate,"/"))
splitend[,3]<-paste("20",splitend[,3],sep="")
OIE$OBEndDate<-paste(splitend[,3],splitend[,2],splitend[,1],sep="-")

splitstart<-do.call("rbind",strsplit(OIEreduced$OBStartDate,"/"))
splitstart[,3]<-paste("20",splitstart[,3],sep="")
OIEreduced$OBStartDate<-paste(splitstart[,3],splitstart[,2],splitstart[,1],sep="-")

OIE<-rbind(OIE,OIEreduced)
rm(OIEreduced)

OIE$OBStartDate<-as.POSIXct(OIE$OBStartDate,format="%Y-%m-%d")
OIE$OBEndDate<-as.POSIXct(OIE$OBEndDate,format="%Y-%m-%d")

##convert to spatial object

coordinates(OIE) <- ~ Longitude + Latitude
proj4string(OIE) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"

lengthlats<-rep(NA,nrow(OIE))
lengthlongs<-rep(NA,nrow(OIE))

for (i in 1:nrow(OIE)) {
  lengthlats[i]<-length(strsplit(as.character(coordinates(OIE)[i,1]),"[.]",perl=T)[[1]])
}

for (i in 1:nrow(OIE)) {
  lengthlongs[i]<-length(strsplit(as.character(coordinates(OIE)[i,2]),"[.]",perl=T)[[1]])
}

NULLlats<-which(lengthlats==1)
NULLlongs<-which(lengthlongs==1)

OIE$LatDecimals<-nchar(do.call("rbind",strsplit(as.character(coordinates(OIE)[,1]),"[.]",perl=T))[,2])
OIE$LatDecimals[NULLlats]<-NA
OIE$LongDecimals<-nchar(do.call("rbind",strsplit(as.character(coordinates(OIE)[,2]),"[.]",perl=T))[,2])
OIE$LongDecimals[NULLlongs]<-NA

rm(lengthlats, lengthlongs, NULLlats, NULLlongs)

targets<-data.frame(as.character(rep(c(2000:2018),each=4)),rep(c("03","06","09","12"),times=19),as.character(c(20,21,22,21,20,21,22,21,20,21,23,22,21,21,23,22,20,21,22,21,20,21,22,21,20,21,23,22,21,21,23,22,20,20,22,21,20,21,22,21,20,21,23,21,20,21,23,22,20,20,22,21,20,21,22,21,20,21,23,21,20,21,23,22,20,20,22,21,20,21,22,21,20,21,23,21)))
targets<-targets[-nrow(targets),]
dates<-paste(targets[,1],"-",targets[,2],"-",targets[,3],sep="")
dates<-c("1999-12-22",dates)
dates<-as.Date(dates)

datesforsearch<-data.frame(rep(NA,75),start=dates[c(1:75)],end=dates[c(2:76)])
colnames(datesforsearch)<-c("mid","start","end")

datesforsearch$mid<-datesforsearch$start+floor((datesforsearch$end-datesforsearch$start)/2)
datesforsearch$mid<-as.POSIXct(datesforsearch$mid,format="%Y-%m-%d")

lat<-rep(NA,nrow(seqinfo))
long<-rep(NA,nrow(seqinfo))
count<-rep(0,nrow(seqinfo))

for (i in 1:(nrow(datesforsearch)-1)) {
  ##test if points in known administrative region
  working<-OIE
  
  ##filter on dates
  start<-datesforsearch$mid[i]
  end<-datesforsearch$mid[(i+1)]
  workingcom<-working[!is.na(working$OBEndDate),]
  workingred<-working[is.na(working$OBEndDate),]
  workingcom<-workingcom[(workingcom$OBStartDate<=start & workingcom$OBEndDate>=end)|(workingcom$OBStartDate>=start & workingcom$OBStartDate<=end)|(workingcom$OBEndDate>=start & workingcom$OBEndDate<=end),]
  workingred<-workingred[workingred$OBStartDate<=start,]
  rm(working)
  if (nrow(workingcom)!=0|nrow(workingred)!=0) {
    working<-rbind(workingcom,workingred)
  }
  rm(workingred,workingcom)
  sppolygons<-list()
  
  ##generate polygons
  if (exists("working")) {
    if (nrow(working)!=0) {
      for (j in 1:nrow(working)) {
        coordsmat<-matrix(NA,4,2)
        if (!is.na(working$LatDecimals[j]) & !is.na(working$LongDecimals[j])) {
          coordsmat[1,1]<-coordinates(working)[j,1]+as.numeric(paste("0.", strrep("0",working$LatDecimals[j]),"5",sep=""))
          coordsmat[1,2]<-coordinates(working)[j,2]+as.numeric(paste("0.", strrep("0",working$LongDecimals[j]),"5",sep=""))
          coordsmat[2,1]<-coordinates(working)[j,1]+as.numeric(paste("0.", strrep("0",working$LatDecimals[j]),"5",sep=""))
          coordsmat[2,2]<-coordinates(working)[j,2]-as.numeric(paste("0.", strrep("0",working$LongDecimals[j]),"5",sep=""))
          coordsmat[3,1]<-coordinates(working)[j,1]-as.numeric(paste("0.", strrep("0",working$LatDecimals[j]),"5",sep=""))
          coordsmat[3,2]<-coordinates(working)[j,2]-as.numeric(paste("0.", strrep("0",working$LongDecimals[j]),"5",sep=""))
          coordsmat[4,1]<-coordinates(working)[j,1]-as.numeric(paste("0.", strrep("0",working$LatDecimals[j]),"5",sep=""))
          coordsmat[4,2]<-coordinates(working)[j,2]+as.numeric(paste("0.", strrep("0",working$LongDecimals[j]),"5",sep=""))
        } else if (!is.na(working$LatDecimals[j]) & is.na(working$LongDecimals[j])) {
          coordsmat[1,1]<-coordinates(working)[j,1]+as.numeric(paste("0.", strrep("0",working$LatDecimals[j]),"5",sep=""))
          coordsmat[1,2]<-coordinates(working)[j,2]+0.5
          coordsmat[2,1]<-coordinates(working)[j,1]+as.numeric(paste("0.", strrep("0",working$LatDecimals[j]),"5",sep=""))
          coordsmat[2,2]<-coordinates(working)[j,2]-0.5
          coordsmat[3,1]<-coordinates(working)[j,1]-as.numeric(paste("0.", strrep("0",working$LatDecimals[j]),"5",sep=""))
          coordsmat[3,2]<-coordinates(working)[j,2]-0.5
          coordsmat[4,1]<-coordinates(working)[j,1]-as.numeric(paste("0.", strrep("0",working$LatDecimals[j]),"5",sep=""))
          coordsmat[4,2]<-coordinates(working)[j,2]+0.5
        } else if (is.na(working$LatDecimals[j]) & !is.na(working$LongDecimals[j])) {
          coordsmat[1,1]<-coordinates(working)[j,1]+0.5
          coordsmat[1,2]<-coordinates(working)[j,2]+as.numeric(paste("0.", strrep("0",working$LongDecimals[j]),"5",sep=""))
          coordsmat[2,1]<-coordinates(working)[j,1]+0.5
          coordsmat[2,2]<-coordinates(working)[j,2]-as.numeric(paste("0.", strrep("0",working$LongDecimals[j]),"5",sep=""))
          coordsmat[3,1]<-coordinates(working)[j,1]-0.5
          coordsmat[3,2]<-coordinates(working)[j,2]-as.numeric(paste("0.", strrep("0",working$LongDecimals[j]),"5",sep=""))
          coordsmat[4,1]<-coordinates(working)[j,1]-0.5
          coordsmat[4,2]<-coordinates(working)[j,2]+as.numeric(paste("0.", strrep("0",working$LongDecimals[j]),"5",sep=""))
        } else {
          coordsmat[1,1]<-coordinates(working)[j,1]+0.5
          coordsmat[1,2]<-coordinates(working)[j,2]+0.5
          coordsmat[2,1]<-coordinates(working)[j,1]+0.5
          coordsmat[2,2]<-coordinates(working)[j,2]-0.5
          coordsmat[3,1]<-coordinates(working)[j,1]-0.5
          coordsmat[3,2]<-coordinates(working)[j,2]-0.5
          coordsmat[4,1]<-coordinates(working)[j,1]-0.5
          coordsmat[4,2]<-coordinates(working)[j,2]+0.5
        }
        sppolygons[[j]]<-Polygons(list(Polygon(coordsmat)),j)
      }
      sppolygons<-SpatialPolygons(sppolygons)
      proj4string(sppolygons)<-CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
      
      #EPSG:3035
      sppolygons<-spTransform(sppolygons,CRSobj = CRS("+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))
      
      latlon<-data.frame(Lat=rep(NA,length(sppolygons)),Lon=rep(NA,length(sppolygons)))
      average<-data.frame(x=rep(NA,50),y=rep(NA,50))
      count[i]<-length(sppolygons)
      for (q in 1:50) {
        for (j in 1:length(sppolygons)) {
          latlon[j,]<-spsample(sppolygons[j],1,"random", iter=100)@coords
          print(paste(i,q,round(j/length(sppolygons)*100,2),sep=" "))
        }
        average[q,]<-colMeans(latlon)
      }
      temppoints<-SpatialPoints(t(as.matrix(colMeans(average))),proj4string=CRS("+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))
      temppoints<-spTransform(temppoints,CRSobj = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
      lat[i]<-temppoints@coords[1,2]
      long[i]<-temppoints@coords[1,1]
    }
  }
  rm(sppolygons,start,end,working,latlon,average,temppoints,coordsmat)
}

final<-data.frame(start=datesforsearch$mid[1:(length(datesforsearch$mid)-1)],end=datesforsearch$mid[2:length(datesforsearch$mid)],mid=datesforsearch$end[1:(length(datesforsearch$mid)-1)],long=long[c(1:(length(datesforsearch$mid)-1))],lat=lat[c(1:(length(datesforsearch$mid)-1))],count=count[c(1:(length(datesforsearch$mid)-1))])
final$mid<-as.Date(final$mid)

abovebelow<-rep(NA,nrow(final))
for(i in 1:nrow(final)) {
  above<-0
  below<-0
  if (i==1) {
    if (is.na(final$long[i])) {
      abovebelow[i]<-"below"
    }
  } else {
    if (is.na(final$long[i])) {
      j<-i
      while (is.na(final$long[j])&j>=1) {
        above<-above+1
        if (j==1) {
          above<-Inf
          break()
        } else {
          j<-j-1
        }
      }
      j<-i
      while (is.na(final$long[j])&j<(nrow(final)+1)) {
        below<-below+1
        j<-j+1
      }
      if (above<=below) {
        abovebelow[i]<-"above"
      } else {
        abovebelow[i]<-"below"
      }
    }
  }
}

for (i in 1:nrow(final)) {
  if (!is.na(abovebelow[i])) {
    if (abovebelow[i]=="above") {
      j<-i
      while (is.na(final$long[j])) {
        j<-j-1
      }
      final$lat[i]<-final$lat[j]
      final$long[i]<-final$long[j]
    } else if (abovebelow[i]=="below") {
      j<-i
      while (is.na(final$long[j])) {
        j<-j+1
      }
      final$lat[i]<-final$lat[j]
      final$long[i]<-final$long[j]
    }
  }
}

final$daylength<-daylength(final$lat,final$mid)
final$start<-date_decimal(final$start)
final$end<-decimal_date(final$end)

write.csv(final,"~/Desktop/covariatelocations.csv",row.names = F)
