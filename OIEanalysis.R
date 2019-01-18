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

##read in sequence IDs, KMLs and dates
#seqinfo<-read.csv()
seqinfo<-data.frame(2006,NA,NA)
colnames(seqinfo)<-c("Year","Month","Day")
seqinfo$Code<-"BEL"
seqinfo$Species<-"Cattle"

##iterate through sequences
for (i in 1:nrow(seqinfo)) {
  lostprobability<-0
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
    if (nchar(as.character(seqinfo$Month[i]))==1) {
      start<-paste(seqinfo$Year[i],"-0",seqinfo$Month[i],"-01",sep="")
      if (seqinfo$Month[i]!=9) {
        end<-paste(seqinfo$Year[i],"-0",seqinfo$Month[i]+1,"-01",sep="")
      } else {
        end<-paste(seqinfo$Year[i]+1,"-10-01",sep="")
      }
    } else {
      start<-paste(seqinfo$Year[i],"-",seqinfo$Month[i],"-01",sep="")
      if (seqinfo$Month[i]!=12) {
        end<-paste(seqinfo$Year[i],"-",seqinfo$Month[i]+1,"-01",sep="")
      } else {
        end<-paste(seqinfo$Year[i]+1,"-01-01",sep="")
      }
    }
    workingcom<-working[!is.na(working$OBEndDate),]
    workingred<-working[is.na(working$OBEndDate),]
    workingcom<-workingcom[(workingcom$OBStartDate<=start & workingcom$OBEndDate>=end)|(workingcom$OBStartDate>=start & workingcom$OBStartDate<=end)|(workingcom$OBEndDate>=start & workingcom$OBEndDate<=end),]
    workingred<-workingred[workingred$OBStartDate<=start,]
    working<-rbind(workingcom,workingred)
    rm(workingred,workingcom)
  } else {
    if (nchar(as.character(seqinfo$Month[i]))!=1&nchar(as.character(seqinfo$Day[i]))!=1) {
      date<-paste(seqinfo$Year[i],"-",seqinfo$Month[i],"-",seqinfo$Day[i],sep="")
    } else if (nchar(as.character(seqinfo$Month[i]))==1&nchar(as.character(seqinfo$Day[i]))!=1) {
      date<-paste(seqinfo$Year[i],"-0",seqinfo$Month[i],"-",seqinfo$Day[i],sep="")
    } else if (nchar(as.character(seqinfo$Month[i]))!=1&nchar(as.character(seqinfo$Day[i]))==1) {
      date<-paste(seqinfo$Year[i],"-",seqinfo$Month[i],"-0",seqinfo$Day[i],sep="")
    } else {
      date<-paste(seqinfo$Year[i],"-0",seqinfo$Month[i],"-0",seqinfo$Day[i],sep="")
    }
    workingcom<-working[!is.na(working$OBEndDate),]
    workingred<-working[is.na(working$OBEndDate),]
    workingcom<-workingcom[(workingcom$OBStartDate<=date & workingcom$OBEndDate>=date),]
    workingred<-workingred[workingred$OBStartDate<=date,]
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
  
  ##generate polygons adding the implicit uncertainty in the OIE records
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
    ##give WGS84 CRS
    proj4string(sppolygons)<-CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
    ##remove areas excluded from known administrative region
    sppolygons<-gIntersection(sppolygons,KML,byid = T)
  }
  
  
  ##generate probabilities - requires finessing when filtering by admin region removes polygons
  prob<-working@data$Cases/sum(working@data$Cases)
  ##change to equal area CRS - depends on input region
  proj<-CRS("+proj=lcc +lat_1=49.83333333333334 +lat_2=51.16666666666666 +lat_0=50.797815 +lon_0=4.359215833333333 +x_0=649328 +y_0=665262 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")
  
  sppolygons<-spTransform(sppolygons,proj)
  
  ##for testing purposes keep a record to minimise reruns
  fortesting<-sppolygons
  sppolygons<-fortesting
  
  ##initialise pointer and markers
  marker1<-0
  marker2<-0
  i<-1
  
  ##while less than the length of the polygons object
  while(i<length(sppolygons)) {
    #search from the target for intersecting polygons
    for (j in (i+1):length(sppolygons)) {
      if (!is.null(gIntersection(sppolygons[i],sppolygons[j]))) {
        ##if the intersection includes points and lines, as well as a polygon on the basis that two polygons that differ only 
        ##by lines and points are equal from the probability perspective
        if (class(gIntersection(sppolygons[i],sppolygons[j]))=="SpatialCollections") {
          ##check for ring objects if present error out
          if (!is.null(gIntersection(sppolygons[i],sppolygons[j])@polyobj)&!is.null(gIntersection(sppolygons[i],sppolygons[j])@ringobj)) {
            stop("Ring object")
          }
          ##if no ring objects
          if (!is.null(gIntersection(sppolygons[i],sppolygons[j])@polyobj)&is.null(gIntersection(sppolygons[i],sppolygons[j])@ringobj)) {
            ##check area of intersection is not zero for safety
            if (area(gIntersection(sppolygons[i],sppolygons[j])@polyobj)!=0) {
              print(paste(i,j))
              ##see if polygons are internal to one another
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
              ##if implied equal
              if (marker1==0&marker2==0) {
                ##take the intersection
                temp1<-gIntersection(sppolygons[i],sppolygons[j])@polyobj
                ##record the areas
                area1<-area(sppolygons[i])
                area2<-area(sppolygons[j])
                ##initialise object for new probabilities
                workingprob1<-rep(NA,length(temp1@polygons[[1]]@Polygons))

                ##initialise object for new polygons
                bindingpolygons1<-list()
                ##save new polygons and correct probability
                for (q in 1:length(temp1@polygons[[1]]@Polygons)) {
                  workingprob1[q]<-prob[i]*(temp1@polygons[[1]]@Polygons[[q]]@area/area1)+prob[j]*(temp1@polygons[[1]]@Polygons[[q]]@area/area2)
                  bindingpolygons1[[q]]<-Polygons(list(temp1@polygons[[1]]@Polygons[[q]]),q)
                }
                bindingpolygons1<-SpatialPolygons(bindingpolygons1)
                proj4string(bindingpolygons1)<-proj

                ##add probability to the list and remove original probabilities
                prob<-c(prob,workingprob1)
                prob<-prob[-c(i,j)]

                ##add polygon to the list and remove original polygons
                sppolygons<-c(sppolygons,bindingpolygons1)
                sppolygons<-do.call(bind, sppolygons)
                sppolygons<-SpatialPolygonsDataFrame(sppolygons,data.frame(q=c(1:length(sppolygons)),match.ID = F))
                print(nrow(sppolygons))
                sppolygons<-sppolygons[-c(i,j),]
                print(nrow(sppolygons))
                sppolygons<-SpatialPolygons(sppolygons@polygons,proj4string=sppolygons@proj4string)

                ##reset iterator as polygon at i now different
                i<-i-1
              }

              ##polygon validity checking
              
              ##record areas of all polygons
              areas<-rep(NA,length(sppolygons))
              for (m in 1:length(sppolygons)) {
                areas[m]<-area(sppolygons[m])
              }

              ##subset out polygons with areas of less than 1m^2
              sppolygons<-SpatialPolygonsDataFrame(sppolygons,data.frame(q=c(1:length(sppolygons))),match.ID = F)
              sppolygons<-sppolygons[(areas>1),]
              sppolygons<-SpatialPolygons(sppolygons@polygons,proj4string=sppolygons@proj4string)

              ##subset out probabilities of polygons with areas of less than 1m^2 and record lost probability
              lostprobability<-lostprobability+sum(prob[(areas<1)])
              prob<-prob[(areas>1)]

              ##find any polygons with invalid geometries
              ######very very hacky - change
              valid<-rep(T,length(sppolygons))
              for (m in 1:length(sppolygons)) {
                if (!gIsValid(sppolygons[m])) {
                  valid[m]<-F
                  print("Invalid geometry detected")
                }
              }

              ##remove polygons with invalid geometries
              sppolygons<-SpatialPolygonsDataFrame(sppolygons,data.frame(q=c(1:length(sppolygons))),match.ID = F)
              sppolygons<-sppolygons[valid,]
              sppolygons<-SpatialPolygons(sppolygons@polygons,proj4string=sppolygons@proj4string)

              ##remove probabilites of polygons with invalid geometries and record lost probability
              lostprobability<-lostprobability+sum(prob[!valid])
              prob<-prob[valid]

              ##break to next level
              break()
            }
          }
        }
        ##if the intersection is just a polygon
        if (class(gIntersection(sppolygons[i],sppolygons[j]))=="SpatialPolygons") {
          ##check area of intersection is not zero for safety
          if (area(gIntersection(sppolygons[i],sppolygons[j]))!=0) {
            print(paste(i,j))
            ##see if polygons are internal to one another
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
            ##if implied equal
            if (marker1==0&marker2==0) {
              ##take the intersection
              temp1<-gIntersection(sppolygons[i],sppolygons[j])
              ##record the areas
              area1<-area(sppolygons[i])
              area2<-area(sppolygons[j])
              ##initialise object for new probabilities
              workingprob1<-rep(NA,length(temp1@polygons[[1]]@Polygons))
              
              ##initialise object for new polygons
              bindingpolygons1<-list()
              ##save new polygons and correct probability
              for (q in 1:length(temp1@polygons[[1]]@Polygons)) {
                workingprob1[q]<-prob[i]*(temp1@polygons[[1]]@Polygons[[q]]@area/area1)+prob[j]*(temp1@polygons[[1]]@Polygons[[q]]@area/area2)
                bindingpolygons1[[q]]<-Polygons(list(temp1@polygons[[1]]@Polygons[[q]]),q)
              }
              bindingpolygons1<-SpatialPolygons(bindingpolygons1)
              proj4string(bindingpolygons1)<-proj
              
              ##add probability to the list and remove original probabilities
              prob<-c(prob,workingprob1)
              prob<-prob[-c(i,j)]
              
              ##add polygon to the list and remove original polygons
              sppolygons<-c(sppolygons,bindingpolygons1)
              sppolygons<-do.call(bind, sppolygons)
              sppolygons<-SpatialPolygonsDataFrame(sppolygons,data.frame(q=c(1:length(sppolygons)),match.ID = F))
              print(nrow(sppolygons))
              sppolygons<-sppolygons[-c(i,j),]
              print(nrow(sppolygons))
              sppolygons<-SpatialPolygons(sppolygons@polygons,proj4string=sppolygons@proj4string)
              
              ##reset iterator as polygon at i now different
              i<-i-1
            }
            
            ##polygon validity checking
            
            ##record areas of all polygons
            areas<-rep(NA,length(sppolygons))
            for (m in 1:length(sppolygons)) {
              areas[m]<-area(sppolygons[m])
            }
            
            ##subset out polygons with areas of less than 1m^2
            sppolygons<-SpatialPolygonsDataFrame(sppolygons,data.frame(q=c(1:length(sppolygons))),match.ID = F)
            sppolygons<-sppolygons[(areas>1),]
            sppolygons<-SpatialPolygons(sppolygons@polygons,proj4string=sppolygons@proj4string)
            
            ##subset out probabilities of polygons with areas of less than 1m^2 and record lost probability
            lostprobability<-lostprobability+sum(prob[(areas<1)])
            prob<-prob[(areas>1)]
            
            ##find any polygons with invalid geometries
            ######very very hacky - change
            valid<-rep(T,length(sppolygons))
            for (m in 1:length(sppolygons)) {
              if (!gIsValid(sppolygons[m])) {
                valid[m]<-F
                print("Invalid geometry detected")
              }
            }
            
            ##remove polygons with invalid geometries
            sppolygons<-SpatialPolygonsDataFrame(sppolygons,data.frame(q=c(1:length(sppolygons))),match.ID = F)
            sppolygons<-sppolygons[valid,]
            sppolygons<-SpatialPolygons(sppolygons@polygons,proj4string=sppolygons@proj4string)
            
            ##remove probabilites of polygons with invalid geometries and record lost probability
            lostprobability<-lostprobability+sum(prob[!valid])
            prob<-prob[valid]
            
            ##break to next level
            break()
          }
        }
      }
    }
    ##reset markers update pointer
    marker1<-0
    marker2<-0
    i<-i+1
  }
  
  
  
  fortesting<-sppolygons
  sppolygons<-fortesting
  marker1<-0
  marker2<-0
  i<-1
  while(i<length(sppolygons)) {
    for (j in (i+1):length(sppolygons)) {
      if (!is.null(gIntersection(sppolygons[i],sppolygons[j]))) {
        if (class(gIntersection(sppolygons[i],sppolygons[j]))=="SpatialCollections") {
          if (!is.null(gIntersection(sppolygons[i],sppolygons[j])@polyobj)&!is.null(gIntersection(sppolygons[i],sppolygons[j])@ringobj)) { 
            stop("Ring object")
          }
          if (!is.null(gIntersection(sppolygons[i],sppolygons[j])@polyobj)&is.null(gIntersection(sppolygons[i],sppolygons[j])@ringobj)) {
            if (area(gIntersection(sppolygons[i],sppolygons[j])@polyobj)!=0) {
              print(paste(i,j))
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
                print("Double difference")
                temp1<-gDifference(sppolygons[i],sppolygons[j])
                temp2<-gDifference(sppolygons[j],sppolygons[i])
                temp3<-gIntersection(sppolygons[i],sppolygons[j])@polyobj
                area1<-area(sppolygons[i])
                area2<-area(sppolygons[j])
                workingprob1<-rep(NA,length(temp1@polygons[[1]]@Polygons))
                workingprob2<-rep(NA,length(temp2@polygons[[1]]@Polygons))
                workingprob3<-rep(NA,length(temp3@polygons[[1]]@Polygons))
                
                bindingpolygons1<-list()
                for (q in 1:length(temp1@polygons[[1]]@Polygons)) {
                  workingprob1[q]<-prob[i]*(temp1@polygons[[1]]@Polygons[[q]]@area/area1)
                  bindingpolygons1[[q]]<-Polygons(list(temp1@polygons[[1]]@Polygons[[q]]),q)
                }
                bindingpolygons1<-SpatialPolygons(bindingpolygons1)
                proj4string(bindingpolygons1)<-proj
                
                bindingpolygons2<-list()
                for (q in 1:length(temp2@polygons[[1]]@Polygons)) {
                  workingprob2[q]<-prob[j]*(temp2@polygons[[1]]@Polygons[[q]]@area/area2)
                  bindingpolygons2[[q]]<-Polygons(list(temp2@polygons[[1]]@Polygons[[q]]),q)
                }
                bindingpolygons2<-SpatialPolygons(bindingpolygons2)
                proj4string(bindingpolygons2)<-proj
                
                bindingpolygons3<-list()
                for (q in 1:length(temp3@polygons[[1]]@Polygons)) {
                  workingprob3[q]<-prob[i]*(temp3@polygons[[1]]@Polygons[[q]]@area/area1)+prob[j]*(temp3@polygons[[1]]@Polygons[[q]]@area/area2)
                  bindingpolygons3[[q]]<-Polygons(list(temp3@polygons[[1]]@Polygons[[q]]),q)
                }
                bindingpolygons3<-SpatialPolygons(bindingpolygons3)
                proj4string(bindingpolygons3)<-proj
                
                prob<-c(prob,workingprob1,workingprob2,workingprob3)
                prob<-prob[-c(i,j)]
                
                sppolygons<-c(sppolygons,bindingpolygons1,bindingpolygons2,bindingpolygons3)
                sppolygons<-do.call(bind, sppolygons)
                sppolygons<-SpatialPolygonsDataFrame(sppolygons,data.frame(q=c(1:length(sppolygons)),match.ID = F))
                print(length(sppolygons))
                sppolygons<-sppolygons[-c(i,j),]
                print(length(sppolygons))
                sppolygons<-SpatialPolygons(sppolygons@polygons,proj4string=sppolygons@proj4string)
                ##if j inside i
              } else if (marker1==1&marker2==0) {
                print("A difference")
                temp1<-gDifference(sppolygons[i],sppolygons[j])
                temp2<-gIntersection(sppolygons[i],sppolygons[j])@polyobj
                area1<-area(sppolygons[i])
                area2<-area(sppolygons[j])
                workingprob1<-rep(NA,length(temp1@polygons[[1]]@Polygons))
                workingprob2<-rep(NA,length(temp2@polygons[[1]]@Polygons))
                
                bindingpolygons1<-list()
                for (q in 1:length(temp1@polygons[[1]]@Polygons)) {
                  workingprob1[q]<-prob[i]*(temp1@polygons[[1]]@Polygons[[q]]@area/area1)
                  bindingpolygons1[[q]]<-Polygons(list(temp1@polygons[[1]]@Polygons[[q]]),q)
                }
                bindingpolygons1<-SpatialPolygons(bindingpolygons1)
                proj4string(bindingpolygons1)<-proj
                
                bindingpolygons2<-list()
                for (q in 1:length(temp2@polygons[[1]]@Polygons)) {
                  workingprob2[q]<-prob[i]*(temp2@polygons[[1]]@Polygons[[q]]@area/area1)+prob[j]*(temp2@polygons[[1]]@Polygons[[q]]@area/area2)
                  bindingpolygons2[[q]]<-Polygons(list(temp2@polygons[[1]]@Polygons[[q]]),q)
                }
                bindingpolygons2<-SpatialPolygons(bindingpolygons2)
                proj4string(bindingpolygons2)<-proj
                
                prob<-c(prob,workingprob1,workingprob2)
                prob<-prob[-c(i,j)]
                
                sppolygons<-c(sppolygons,bindingpolygons1,bindingpolygons2)
                sppolygons<-do.call(bind, sppolygons)
                sppolygons<-SpatialPolygonsDataFrame(sppolygons,data.frame(q=c(1:length(sppolygons)),match.ID = F))
                print(length(sppolygons))
                sppolygons<-sppolygons[-c(i,j),]
                print(length(sppolygons))
                sppolygons<-SpatialPolygons(sppolygons@polygons,proj4string=sppolygons@proj4string)
                ##if i inside j
              } else if (marker1==0&marker2==1) {
                print("B difference")
                temp1<-gDifference(sppolygons[j],sppolygons[i])
                temp2<-gIntersection(sppolygons[i],sppolygons[j])@polyobj
                area1<-area(sppolygons[j])
                area2<-area(sppolygons[i])
                workingprob1<-rep(NA,length(temp1@polygons[[1]]@Polygons))
                workingprob2<-rep(NA,length(temp2@polygons[[1]]@Polygons))
                
                bindingpolygons1<-list()
                for (q in 1:length(temp1@polygons[[1]]@Polygons)) {
                  workingprob1[q]<-prob[j]*(temp1@polygons[[1]]@Polygons[[q]]@area/area1)
                  bindingpolygons1[[q]]<-Polygons(list(temp1@polygons[[1]]@Polygons[[q]]),q)
                }
                bindingpolygons1<-SpatialPolygons(bindingpolygons1)
                proj4string(bindingpolygons1)<-proj
                
                bindingpolygons2<-list()
                for (q in 1:length(temp2@polygons[[1]]@Polygons)) {
                  workingprob2[q]<-prob[j]*(temp2@polygons[[1]]@Polygons[[q]]@area/area1)+prob[i]*(temp2@polygons[[1]]@Polygons[[q]]@area/area2)
                  bindingpolygons2[[q]]<-Polygons(list(temp2@polygons[[1]]@Polygons[[q]]),q)
                }
                bindingpolygons2<-SpatialPolygons(bindingpolygons2)
                proj4string(bindingpolygons2)<-proj
                
                prob<-c(prob,workingprob1,workingprob2)
                prob<-prob[-c(i,j)]
                
                sppolygons<-c(sppolygons,bindingpolygons1,bindingpolygons2)
                sppolygons<-do.call(bind, sppolygons)
                sppolygons<-SpatialPolygonsDataFrame(sppolygons,data.frame(q=c(1:length(sppolygons)),match.ID = F))
                print(length(sppolygons))
                sppolygons<-sppolygons[-c(i,j),]
                print(length(sppolygons))
                sppolygons<-SpatialPolygons(sppolygons@polygons,proj4string=sppolygons@proj4string)
              } else {
                print("Equal")
                temp1<-gIntersection(sppolygons[i],sppolygons[j])@polyobj
                area1<-area(sppolygons[i])
                area2<-area(sppolygons[j])
                workingprob1<-rep(NA,length(temp1@polygons[[1]]@Polygons))
                
                bindingpolygons1<-list()
                for (q in 1:length(temp1@polygons[[1]]@Polygons)) {
                  workingprob1[q]<-prob[i]*(temp1@polygons[[1]]@Polygons[[q]]@area/area1)+prob[j]*(temp1@polygons[[1]]@Polygons[[q]]@area/area2)
                  bindingpolygons1[[q]]<-Polygons(list(temp1@polygons[[1]]@Polygons[[q]]),q)
                }
                bindingpolygons1<-SpatialPolygons(bindingpolygons1)
                proj4string(bindingpolygons1)<-proj
                
                prob<-c(prob,workingprob1)
                prob<-prob[-c(i,j)]
                
                sppolygons<-c(sppolygons,bindingpolygons1)
                sppolygons<-do.call(bind, sppolygons)
                sppolygons<-SpatialPolygonsDataFrame(sppolygons,data.frame(q=c(1:length(sppolygons)),match.ID = F))
                sppolygons<-sppolygons[-c(i,j),]
                sppolygons<-SpatialPolygons(sppolygons@polygons,proj4string=sppolygons@proj4string)
              }
              
              areas<-rep(NA,length(sppolygons))
              for (m in 1:length(sppolygons)) {
                areas[m]<-area(sppolygons[m])
              }
              
              sppolygons<-SpatialPolygonsDataFrame(sppolygons,data.frame(q=c(1:length(sppolygons))),match.ID = F)
              sppolygons<-sppolygons[(areas>1),]
              sppolygons<-SpatialPolygons(sppolygons@polygons,proj4string=sppolygons@proj4string)
              
              prob<-prob[(areas>1)]
              
              ######very very hacky - change
              
              valid<-rep(T,length(sppolygons))
              for (m in 1:length(sppolygons)) {
                if (!gIsValid(sppolygons[m])) {
                  valid[m]<-F
                  print("Invalid geometry detected")
                }
              }
              
              sppolygons<-SpatialPolygonsDataFrame(sppolygons,data.frame(q=c(1:length(sppolygons))),match.ID = F)
              sppolygons<-sppolygons[valid,]
              sppolygons<-SpatialPolygons(sppolygons@polygons,proj4string=sppolygons@proj4string)
              
              prob<-prob[valid]
              
              i<-i-1
              break()
            }
          }
        }
        if (class(gIntersection(sppolygons[i],sppolygons[j]))=="SpatialPolygons") {
          if (area(gIntersection(sppolygons[i],sppolygons[j]))!=0) {
            print(paste(i,j))
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
              print("Double difference")
              temp1<-gDifference(sppolygons[i],sppolygons[j])
              temp2<-gDifference(sppolygons[j],sppolygons[i])
              temp3<-gIntersection(sppolygons[i],sppolygons[j])
              area1<-area(sppolygons[i])
              area2<-area(sppolygons[j])
              workingprob1<-rep(NA,length(temp1@polygons[[1]]@Polygons))
              workingprob2<-rep(NA,length(temp2@polygons[[1]]@Polygons))
              workingprob3<-rep(NA,length(temp3@polygons[[1]]@Polygons))
              
              bindingpolygons1<-list()
              for (q in 1:length(temp1@polygons[[1]]@Polygons)) {
                workingprob1[q]<-prob[i]*(temp1@polygons[[1]]@Polygons[[q]]@area/area1)
                bindingpolygons1[[q]]<-Polygons(list(temp1@polygons[[1]]@Polygons[[q]]),q)
              }
              bindingpolygons1<-SpatialPolygons(bindingpolygons1)
              proj4string(bindingpolygons1)<-proj
              
              bindingpolygons2<-list()
              for (q in 1:length(temp2@polygons[[1]]@Polygons)) {
                workingprob2[q]<-prob[j]*(temp2@polygons[[1]]@Polygons[[q]]@area/area2)
                bindingpolygons2[[q]]<-Polygons(list(temp2@polygons[[1]]@Polygons[[q]]),q)
              }
              bindingpolygons2<-SpatialPolygons(bindingpolygons2)
              proj4string(bindingpolygons2)<-proj
              
              bindingpolygons3<-list()
              for (q in 1:length(temp3@polygons[[1]]@Polygons)) {
                workingprob3[q]<-prob[i]*(temp3@polygons[[1]]@Polygons[[q]]@area/area1)+prob[j]*(temp3@polygons[[1]]@Polygons[[q]]@area/area2)
                bindingpolygons3[[q]]<-Polygons(list(temp3@polygons[[1]]@Polygons[[q]]),q)
              }
              bindingpolygons3<-SpatialPolygons(bindingpolygons3)
              proj4string(bindingpolygons3)<-proj
              
              prob<-c(prob,workingprob1,workingprob2,workingprob3)
              prob<-prob[-c(i,j)]
              
              sppolygons<-c(sppolygons,bindingpolygons1,bindingpolygons2,bindingpolygons3)
              sppolygons<-do.call(bind, sppolygons)
              sppolygons<-SpatialPolygonsDataFrame(sppolygons,data.frame(q=c(1:length(sppolygons))),match.ID = F)
              print(length(sppolygons))
              sppolygons<-sppolygons[-c(i,j),]
              print(length(sppolygons))
              sppolygons<-SpatialPolygons(sppolygons@polygons,proj4string=sppolygons@proj4string)
              ##if j inside i
            } else if (marker1==1&marker2==0) {
              print("A difference")
              temp1<-gDifference(sppolygons[i],sppolygons[j])
              temp2<-gIntersection(sppolygons[i],sppolygons[j])
              area1<-area(sppolygons[i])
              area2<-area(sppolygons[j])
              workingprob1<-rep(NA,length(temp1@polygons[[1]]@Polygons))
              workingprob2<-rep(NA,length(temp2@polygons[[1]]@Polygons))
              
              bindingpolygons1<-list()
              for (q in 1:length(temp1@polygons[[1]]@Polygons)) {
                workingprob1[q]<-prob[i]*(temp1@polygons[[1]]@Polygons[[q]]@area/area1)
                bindingpolygons1[[q]]<-Polygons(list(temp1@polygons[[1]]@Polygons[[q]]),q)
              }
              bindingpolygons1<-SpatialPolygons(bindingpolygons1)
              proj4string(bindingpolygons1)<-proj
              
              bindingpolygons2<-list()
              for (q in 1:length(temp2@polygons[[1]]@Polygons)) {
                workingprob2[q]<-prob[i]*(temp2@polygons[[1]]@Polygons[[q]]@area/area1)+prob[j]*(temp2@polygons[[1]]@Polygons[[q]]@area/area2)
                bindingpolygons2[[q]]<-Polygons(list(temp2@polygons[[1]]@Polygons[[q]]),q)
              }
              bindingpolygons2<-SpatialPolygons(bindingpolygons2)
              proj4string(bindingpolygons2)<-proj
              
              prob<-c(prob,workingprob1,workingprob2)
              prob<-prob[-c(i,j)]
              
              sppolygons<-c(sppolygons,bindingpolygons1,bindingpolygons2)
              sppolygons<-do.call(bind, sppolygons)
              sppolygons<-SpatialPolygonsDataFrame(sppolygons,data.frame(q=c(1:length(sppolygons)),match.ID = F))
              print(length(sppolygons))
              sppolygons<-sppolygons[-c(i,j),]
              print(length(sppolygons))
              sppolygons<-SpatialPolygons(sppolygons@polygons,proj4string=sppolygons@proj4string)
              ##if i inside j
            } else if (marker1==0&marker2==1) {
              print("B difference")
              temp1<-gDifference(sppolygons[j],sppolygons[i])
              temp2<-gIntersection(sppolygons[i],sppolygons[j])
              area1<-area(sppolygons[j])
              area2<-area(sppolygons[i])
              workingprob1<-rep(NA,length(temp1@polygons[[1]]@Polygons))
              workingprob2<-rep(NA,length(temp2@polygons[[1]]@Polygons))
              
              bindingpolygons1<-list()
              for (q in 1:length(temp1@polygons[[1]]@Polygons)) {
                workingprob1[q]<-prob[j]*(temp1@polygons[[1]]@Polygons[[q]]@area/area1)
                bindingpolygons1[[q]]<-Polygons(list(temp1@polygons[[1]]@Polygons[[q]]),q)
              }
              bindingpolygons1<-SpatialPolygons(bindingpolygons1)
              proj4string(bindingpolygons1)<-proj
              
              bindingpolygons2<-list()
              for (q in 1:length(temp2@polygons[[1]]@Polygons)) {
                workingprob2[q]<-prob[j]*(temp2@polygons[[1]]@Polygons[[q]]@area/area1)+prob[i]*(temp2@polygons[[1]]@Polygons[[q]]@area/area2)
                bindingpolygons2[[q]]<-Polygons(list(temp2@polygons[[1]]@Polygons[[q]]),q)
              }
              bindingpolygons2<-SpatialPolygons(bindingpolygons2)
              proj4string(bindingpolygons2)<-proj
              
              prob<-c(prob,workingprob1,workingprob2)
              prob<-prob[-c(i,j)]
              
              sppolygons<-c(sppolygons,bindingpolygons1,bindingpolygons2)
              sppolygons<-do.call(bind, sppolygons)
              sppolygons<-SpatialPolygonsDataFrame(sppolygons,data.frame(q=c(1:length(sppolygons)),match.ID = F))
              print(length(sppolygons))
              sppolygons<-sppolygons[-c(i,j),]
              print(length(sppolygons))
              sppolygons<-SpatialPolygons(sppolygons@polygons,proj4string=sppolygons@proj4string)
            } else {
              print("Equal")
              temp1<-gIntersection(sppolygons[i],sppolygons[j])
              area1<-area(sppolygons[i])
              area2<-area(sppolygons[j])
              workingprob1<-rep(NA,length(temp1@polygons[[1]]@Polygons))
              
              bindingpolygons1<-list()
              for (q in 1:length(temp1@polygons[[1]]@Polygons)) {
                workingprob1[q]<-prob[i]*(temp1@polygons[[1]]@Polygons[[q]]@area/area1)+prob[j]*(temp1@polygons[[1]]@Polygons[[q]]@area/area2)
                bindingpolygons1[[q]]<-Polygons(list(temp1@polygons[[1]]@Polygons[[q]]),q)
              }
              bindingpolygons1<-SpatialPolygons(bindingpolygons1)
              proj4string(bindingpolygons1)<-proj
              
              prob<-c(prob,workingprob1)
              prob<-prob[-c(i,j)]
              
              sppolygons<-c(sppolygons,bindingpolygons1)
              sppolygons<-do.call(bind, sppolygons)
              sppolygons<-SpatialPolygonsDataFrame(sppolygons,data.frame(q=c(1:length(sppolygons)),match.ID = F))
              print(length(sppolygons))
              sppolygons<-sppolygons[-c(i,j),]
              print(length(sppolygons))
              sppolygons<-SpatialPolygons(sppolygons@polygons,proj4string=sppolygons@proj4string)
            }
            
            areas<-rep(NA,length(sppolygons))
            for (m in 1:length(sppolygons)) {
              areas[m]<-area(sppolygons[m])
            }
            
            sppolygons<-SpatialPolygonsDataFrame(sppolygons,data.frame(q=c(1:length(sppolygons))),match.ID = F)
            sppolygons<-sppolygons[(areas>1),]
            sppolygons<-SpatialPolygons(sppolygons@polygons,proj4string=sppolygons@proj4string)
            
            prob<-prob[(areas>1)]
            
            ######very very hacky - change
            
            valid<-rep(T,length(sppolygons))
            for (m in 1:length(sppolygons)) {
              if (!gIsValid(sppolygons[m])) {
                valid[m]<-F
                print("Invalid geometry detected")
              }
            }
            
            sppolygons<-SpatialPolygonsDataFrame(sppolygons,data.frame(q=c(1:length(sppolygons))),match.ID = F)
            sppolygons<-sppolygons[valid,]
            sppolygons<-SpatialPolygons(sppolygons@polygons,proj4string=sppolygons@proj4string)
            
            prob<-prob[valid]
            
            i<-i-1
            break()
          }
        }
      }
    }
    marker1<-0
    marker2<-0
    i<-i+1
  }
}
