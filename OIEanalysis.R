###code for reading in OIE data on a per sequence level and generating 
library(sp)
library(raster)
library(rgdal)
library(rgeos)
library(maptools)
library(spatstat)
library(sf)
library(lwgeom)

OIE<-read.csv("~/Desktop/OIEdataset.csv")[,c(4,5,7,8,10,11,12,14)]
OIE<-OIE[!OIE$Cases%in%0,]
OIE<-OIE[!is.na(OIE$Cases),]

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

#coordinates(OIE)[1,1]-as.numeric(paste("0.", strrep("0",6),"5",sep=""))

##read in sequence IDs, KMLs and dates
#seqinfo<-read.csv()
seqinfo<-data.frame(2006,NA,NA)
colnames(seqinfo)<-c("Year","Month","Day")
seqinfo$Code<-"BEL"
seqinfo$Species<-"Cattle"

for (i in 1:nrow(seqinfo)) {
  KML<-readOGR("~/Desktop/Glasgow Work/KMLs/Belgium.kml","Belgium")
  
  ##test if points in known administrative region
  working<-OIE
  
  ##filter on dates
  if (is.na(seqinfo$Month[i])) {
    start<-paste(seqinfo$Year[i],"-01-01",sep="")
    end<-paste(seqinfo$Year[i]+1,"-01-01",sep="")
    workingcom<-working[!is.na(working$OBEndDate),]
    workingred<-working[is.na(working$OBEndDate),]
    workingcom<-workingcom[(workingcom$OBStartDate<=start & workingcom$OBEndDate>=end)|(workingcom$OBStartDate>=start & workingcom$OBStartDate<=end)|(workingcom$OBEndDate>=start & workingcom$OBEndDate<=end),]
    workingred<-workingred[workingred$OBStartDate<=start,]
    working<-rbind(workingcom,workingred)
    rm(workingred,workingcom)
  } else if (is.na(seqinfo$Day[i])) {
    start<-paste(seqinfo$Year[i],"-",seqinfo$Month[i],"-01",sep="")
    if (seqinfo$Month[i]!=12) {
      end<-paste(seqinfo$Year[i],"-",seqinfo$Month[i]+1,"-01",sep="")
    } else {
      end<-paste(seqinfo$Year[i]+1,"-01-01",sep="")
    }
    workingcom<-working[!is.na(working$OBEndDate),]
    workingred<-working[is.na(working$OBEndDate),]
    workingcom<-workingcom[(workingcom$OBStartDate<=start & workingcom$OBEndDate>=end)|(workingcom$OBStartDate>=start & workingcom$OBStartDate<=end)|(workingcom$OBEndDate>=start & workingcom$OBEndDate<=end),]
    workingred<-workingred[workingred$OBStartDate<=start,]
    working<-rbind(workingcom,workingred)
    rm(workingred,workingcom)
  } else {
    date<-paste(seqinfo$Year[i],"-",seqinfo$Month[i],"-",seqinfo$Day[i],sep="")
    workingcom<-working[!is.na(working$OBEndDate),]
    workingred<-working[is.na(working$OBEndDate),]
    workingcom<-workingcom[(workingcom$OBStartDate<=date & workingcom$OBEndDate>=date),]
    workingred<-workingred[workingred$OBStartDate<=start,]
    working<-rbind(workingcom,workingred)
    rm(workingred,workingcom)
  }
  
  ##filter on country
  #working<-working[as.character(working$Code)%in%as.character(seqinfo$Code),]
  working<-working[as.character(working$Code)%in%"BEL",]
  
  ##filter on species where possible
  if (!is.na(seqinfo$Species[i])) {
    working<-working[grep(as.character(seqinfo$Species[i]),as.character(working$Species),ignore.case = T),]
  }
  
  sppolygons<-list()
  
  ##generate polygons
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
    sppolygons<-gIntersection(sppolygons,KML,byid = T)
  }
}

# sppolygons<-spTransform(sppolygons,CRS("+proj=lcc +lat_1=51.16666723333333 +lat_2=49.8333339 +lat_0=90 +lon_0=4.367486666666666 +x_0=150000.013 +y_0=5400088.438 +ellps=intl +towgs84=-106.869,52.2978,-103.724,0.3366,-0.457,1.8422,-1.2747 +units=m +no_defs"))
# 
# prob<-working@data$Cases/sum(working@data$Cases)
# temptest<-st_as_sf(sppolygons)
# las<-st_intersection(temptest)
# 
# st_area(temptest)
# st_area(las)[do.call("c", lapply(las$origins,function(x){if(1%in%x) {return(1)} else {return(0)}}))]
# 
# do.call("c", lapply(las$origins,function(x){if(1%in%x) {return(1)} else {return(0)}}))
# 


##remove duplicate polygons and add probabilities

#fortesting<-sppolygons
#sppolygons<-fortesting
#sppolygonsdf<-SpatialPolygonsDataFrame(sppolygons,data.frame(Prob=(working@data$Cases/sum(working@data$Cases))),match.ID = F)
prob<-working@data$Cases/sum(working@data$Cases)
i<-1
end<-length(sppolygons)
while (i<end) {
  a1<-sppolygons[i]@polygons[[1]]@area
  for (j in (i+1):length(sppolygons)) {
    a2<-sppolygons[j]@polygons[[1]]@area
    print(paste(i,j))
    if (is.null(intersect(sppolygons[i],sppolygons[j]))) {
      
    } else {
      if (length(intersect(sppolygons[i],sppolygons[j]))==1&length(unique(c(a1,a2,intersect(sppolygons[i],sppolygons[j])@polygons[[1]]@area)))==1) {
        prob[i]<-prob[i]+prob[j]
        sppolygons<-sppolygons[-j]
        prob<-prob[-j]
        i<-i-1
        break()
      }
    }
  }
  end<-length(sppolygons)
  i<-i+1
}

sppolygons<-spTransform(sppolygons,CRS("+proj=lcc +lat_1=49.83333333333334 +lat_2=51.16666666666666 +lat_0=50.797815 +lon_0=4.359215833333333 +x_0=649328 +y_0=665262 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))

##shard polygons
#fortesting<-sppolygons
sppolygons<-fortesting
#temptest<-st_as_sf(sppolygons)
#las<-st_intersection(temptest)

marker1<-0
marker2<-0
i<-1
while(i<length(sppolygons)) {
  for (j in (i+1):length(sppolygons)) {
    if (!is.null(gIntersection(sppolygons[i],sppolygons[j]))) {
      if (class(gIntersection(sppolygons[i],sppolygons[j]))=="SpatialCollections") { 
        if (!is.null(gIntersection(sppolygons[i],sppolygons[j])@polyobj)&is.null(gIntersection(sppolygons[i],sppolygons[j])@ringobj)) {
          if (area(gIntersection(sppolygons[i],sppolygons[j])@polyobj)!=0) {
            if (!is.null(gDifference(sppolygons[i],sppolygons[j]))) {
              if (area(gDifference(sppolygons[i],sppolygons[j]))!=0) {
                marker1<-1
              }
            }  
            if (!is.null(gDifference(sppolygons[j],sppolygons[i]))) {
              if (area(gDifference(sppolygons[j],sppolygons[i]))!=0) {
                marker2<-1
              }
            }
            ##if both have different regions
            if (marker1==1&marker2==1) {
              temp1<-gDifference(sppolygons[i],sppolygons[j])
              temp2<-gDifference(sppolygons[j],sppolygons[i])
              temp3<-gIntersection(sppolygons[i],sppolygons[j])@polyobj
              area1<-area(sppolygons[i])
              area2<-area(sppolygons[j])
              workingprob1<-rep(0,length(temp1))
              workingprob2<-rep(0,length(temp2))
              workingprob3<-rep(0,length(temp3))
              for (q in 1:seq_along(temp1)) {
                workingprob1[q]<-prob[i]*(area(temp1[q])/area1)
              }
              for (q in 1:seq_along(temp2)) {
                workingprob2[q]<-prob[j]*(area(temp2[q])/area2)
              }
              for (q in 1:seq_along(temp3)) {
                workingprob3[q]<-prob[i]*(area(temp3[q])/area1)+prob[j]*(area(temp3[q])/area2)
              }
              prob<-c(prob,workingprob1,workingprob2,workingprob3)
              prob<-prob[-c(i,j)]
              sppolygons<-c(sppolygons,temp1,temp2,temp3)
              sppolygons<-do.call(bind, sppolygons)
              sppolygons<-SpatialPolygonsDataFrame(sppolygons,data.frame(q=c(1:length(sppolygons)),match.ID = F))
              sppolygons<-sppolygons[-c(i,j),]
              sppolygons<-SpatialPolygons(sppolygons@polygons,proj4string=sppolygons@proj4string)
              ##if j inside i
            } else if (marker1==1&marker2==0) {
              temp1<-gDifference(sppolygons[i],sppolygons[j])
              temp2<-gIntersection(sppolygons[i],sppolygons[j])@polyobj
              area1<-area(sppolygons[i])
              area2<-area(sppolygons[j])
              workingprob1<-rep(0,length(temp1))
              workingprob2<-rep(0,length(temp2))
              for (q in 1:seq_along(temp1)) {
                workingprob1[q]<-prob[i]*(area(temp1[q])/area1)
              }
              for (q in 1:seq_along(temp2)) {
                workingprob2[q]<-prob[i]*(area(temp2[q])/area1)+prob[j]*(area(temp2[q])/area2)
              }
              prob<-c(prob,workingprob1,workingprob2)
              sppolygons<-c(sppolygons,temp1,temp2)
              sppolygons<-do.call(bind, sppolygons)
              sppolygons<-SpatialPolygonsDataFrame(sppolygons,data.frame(q=c(1:length(sppolygons)),match.ID = F))
              sppolygons<-sppolygons[-c(i,j),]
              sppolygons<-SpatialPolygons(sppolygons@polygons,proj4string=sppolygons@proj4string)
              ##if i inside j
            } else if (marker1==0&marker2==1) {
              temp1<-gDifference(sppolygons[j],sppolygons[i])
              temp2<-gIntersection(sppolygons[i],sppolygons[j])@polyobj
              area2<-area(sppolygons[i])
              area1<-area(sppolygons[j])
              workingprob1<-rep(0,length(temp1))
              workingprob2<-rep(0,length(temp2))
              for (q in 1:seq_along(temp1)) {
                workingprob1[q]<-prob[j]*(area(temp1[q])/area2)
              }
              for (q in 1:seq_along(temp2)) {
                workingprob2[q]<-prob[i]*(area(temp2[q])/area1)+prob[j]*(area(temp2[q])/area2)
              }
              prob<-c(prob,workingprob1,workingprob2)
              sppolygons<-c(sppolygons,temp1,temp2)
              sppolygons<-do.call(bind, sppolygons)
              sppolygons<-SpatialPolygonsDataFrame(sppolygons,data.frame(q=c(1:length(sppolygons)),match.ID = F))
              sppolygons<-sppolygons[-c(i,j),]
              sppolygons<-SpatialPolygons(sppolygons@polygons,proj4string=sppolygons@proj4string)
            } else {
              temp1<-gIntersection(sppolygons[i],sppolygons[j])@polyobj
              area1<-area(sppolygons[i])
              area2<-area(sppolygons[j])
              workingprob1<-rep(0,length(temp1))
              for (q in 1:seq_along(temp1)) {
                workingprob1[q]<-prob[i]*(area(temp1[q])/area1)+prob[j]*(area(temp1[q])/area2)
              }
              prob<-c(prob,workingprob1)
              sppolygons<-c(sppolygons,temp1)
              sppolygons<-do.call(bind, sppolygons)
              sppolygons<-SpatialPolygonsDataFrame(sppolygons,data.frame(q=c(1:length(sppolygons)),match.ID = F))
              sppolygons<-sppolygons[-c(i,j),]
              sppolygons<-SpatialPolygons(sppolygons@polygons,proj4string=sppolygons@proj4string)
            }
            i<-i-1
            break()
          }
        }
      }
      if (class(gIntersection(sppolygons[i],sppolygons[j]))=="SpatialPolygons") {
        if (area(gIntersection(sppolygons[i],sppolygons[j]))!=0) {
          if (!is.null(gDifference(sppolygons[i],sppolygons[j]))) {
            if (area(gDifference(sppolygons[i],sppolygons[j]))!=0) {
              marker1<-1
            }
          }  
          if (!is.null(gDifference(sppolygons[j],sppolygons[i]))) {
            if (area(gDifference(sppolygons[j],sppolygons[i]))!=0) {
              marker2<-1
            }
          }
          ##if both have different regions
          if (marker1==1&marker2==1) {
            temp1<-gDifference(sppolygons[i],sppolygons[j])
            temp2<-gDifference(sppolygons[j],sppolygons[i])
            temp3<-gIntersection(sppolygons[i],sppolygons[j])
            area1<-area(sppolygons[i])
            area2<-area(sppolygons[j])
            workingprob1<-rep(0,length(temp1))
            workingprob2<-rep(0,length(temp2))
            workingprob3<-rep(0,length(temp3))
            for (q in 1:seq_along(temp1)) {
              workingprob1[q]<-prob[i]*(area(temp1[q])/area1)
            }
            for (q in 1:seq_along(temp2)) {
              workingprob2[q]<-prob[j]*(area(temp2[q])/area2)
            }
            for (q in 1:seq_along(temp3)) {
              workingprob3[q]<-prob[i]*(area(temp3[q])/area1)+prob[j]*(area(temp3[q])/area2)
            }
            prob<-c(prob,workingprob1,workingprob2,workingprob3)
            prob<-prob[-c(i,j)]
            sppolygons<-c(sppolygons,temp1,temp2,temp3)
            sppolygons<-do.call(bind, sppolygons)
            sppolygons<-SpatialPolygonsDataFrame(sppolygons,data.frame(q=c(1:length(sppolygons)),match.ID = F))
            sppolygons<-sppolygons[-c(i,j),]
            sppolygons<-SpatialPolygons(sppolygons@polygons,proj4string=sppolygons@proj4string)
            ##if j inside i
          } else if (marker1==1&marker2==0) {
            temp1<-gDifference(sppolygons[i],sppolygons[j])
            temp2<-gIntersection(sppolygons[i],sppolygons[j])
            area1<-area(sppolygons[i])
            area2<-area(sppolygons[j])
            workingprob1<-rep(0,length(temp1))
            workingprob2<-rep(0,length(temp2))
            for (q in 1:seq_along(temp1)) {
              workingprob1[q]<-prob[i]*(area(temp1[q])/area1)
            }
            for (q in 1:seq_along(temp2)) {
              workingprob2[q]<-prob[i]*(area(temp2[q])/area1)+prob[j]*(area(temp2[q])/area2)
            }
            prob<-c(prob,workingprob1,workingprob2)
            sppolygons<-c(sppolygons,temp1,temp2)
            sppolygons<-do.call(bind, sppolygons)
            sppolygons<-SpatialPolygonsDataFrame(sppolygons,data.frame(q=c(1:length(sppolygons)),match.ID = F))
            sppolygons<-sppolygons[-c(i,j),]
            sppolygons<-SpatialPolygons(sppolygons@polygons,proj4string=sppolygons@proj4string)
            ##if i inside j
          } else if (marker1==0&marker2==1) {
            temp1<-gDifference(sppolygons[j],sppolygons[i])
            temp2<-gIntersection(sppolygons[i],sppolygons[j])
            area2<-area(sppolygons[i])
            area1<-area(sppolygons[j])
            workingprob1<-rep(0,length(temp1))
            workingprob2<-rep(0,length(temp2))
            for (q in 1:seq_along(temp1)) {
              workingprob1[q]<-prob[j]*(area(temp1[q])/area2)
            }
            for (q in 1:seq_along(temp2)) {
              workingprob2[q]<-prob[i]*(area(temp2[q])/area1)+prob[j]*(area(temp2[q])/area2)
            }
            prob<-c(prob,workingprob1,workingprob2)
            sppolygons<-c(sppolygons,temp1,temp2)
            sppolygons<-do.call(bind, sppolygons)
            sppolygons<-SpatialPolygonsDataFrame(sppolygons,data.frame(q=c(1:length(sppolygons)),match.ID = F))
            sppolygons<-sppolygons[-c(i,j),]
            sppolygons<-SpatialPolygons(sppolygons@polygons,proj4string=sppolygons@proj4string)
          } else {
            temp1<-gIntersection(sppolygons[i],sppolygons[j])
            area1<-area(sppolygons[i])
            area2<-area(sppolygons[j])
            workingprob1<-rep(0,length(temp1))
            for (q in 1:seq_along(temp1)) {
              workingprob1[q]<-prob[i]*(area(temp1[q])/area1)+prob[j]*(area(temp1[q])/area2)
            }
            prob<-c(prob,workingprob1)
            sppolygons<-c(sppolygons,temp1)
            sppolygons<-do.call(bind, sppolygons)
            sppolygons<-SpatialPolygonsDataFrame(sppolygons,data.frame(q=c(1:length(sppolygons)),match.ID = F))
            sppolygons<-sppolygons[-c(i,j),]
            sppolygons<-SpatialPolygons(sppolygons@polygons,proj4string=sppolygons@proj4string)
          }
          i<-i-1
          break()
        }
      }
    }
  }
  print(i)
  marker1<-0
  marker2<-0
  i<-i+1
}

m<-0
k<-0
i<-1
while(i<length(sppolygons)) {
  a1<-sppolygons[i]@polygons[[1]]@area
  for (j in (i+1):length(sppolygons)) {
    a2<-sppolygons[j]@polygons[[1]]@area
    if (is.null(intersect(sppolygons[i],sppolygons[j]))) {
      
    } else {
      if (length(intersect(sppolygons[i],sppolygons[j]))==1&length(unique(c(a1,a2,intersect(sppolygons[i],sppolygons[j])@polygons[[1]]@area)))==1) {
      }
      print(paste(i,j))
    }
  }
  i<-i+1
}


----

for(i in 1:(length(sppolygons)-1)) {
  a1<-sppolygons[i,]@polygons[[1]]@area
  for (j in (i+1):length(sppolygons)) {
    a2<-sppolygons[j,]@polygons[[1]]@area
    if (is.null(intersect(sppolygons[i,],sppolygons[j,]))) {
      
    } else {
      m<-m+1
      if (length(intersect(sppolygons[i,],sppolygons[j,]))==1&length(unique(c(a1,a2,intersect(sppolygons[i,],sppolygons[j,])[1,]@polygons[[1]]@area)))==1) {
        k<-k+1
      }
      print(paste(m,k))
    }
  }
}

