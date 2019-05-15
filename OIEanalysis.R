###code for reading in OIE data on a per sequence level and generating 
library(sp)
library(raster)
library(rgdal)
library(rgeos)
library(maptools)
library(spatstat)
library(sf)
library(lwgeom)

recursivesimplify<-function (x) {
  ##start pointer
  b<-1
  ##save temporary variable
  target<-x
  while (b<=nrow(target)) {
    ##if the geometry if MULTIPOLYGON cast to polygon
    if (class(target$geometry[b])[1]=="sfc_MULTIPOLYGON") {
      ##join
      target<-rbind(target,st_cast(target[b,],"POLYGON"))
      ##remove original
      target<-target[-b,]
      b<-b-1
    }
    ##if the geometry if GEOMETRYCOLLECTION cast to separate can call function recursively again
    if (class(target$geometry[b])[1]=="sfc_GEOMETRYCOLLECTION") {
      temp<-st_cast(target[b,])
      ##remove area 0 area geometries
      temp<-temp[as.numeric(st_area(temp))!=0,]
      temp<-recursivesimplify(temp)
      #join
      target<-rbind(target,temp)
      #remove original
      target<-target[-b,]
      b<-b-1
    }
    b<-b+1
  }
  ##return temporary variable
  target<-target[as.numeric(st_area(target))!=0,]
  return(target)
}

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
seqinfo<-read.csv("~/Dropbox/BTV-8/KMLgenerator.csv",stringsAsFactors=F)

##iterate through sequences
for (i in 1:nrow(seqinfo)) {
  lostprobability<-0
  KML<-st_zm(read_sf(paste("~/Desktop/Glasgow Work/KMLs/",seqinfo$KML[i],sep=""),strsplit(seqinfo$KML[i],".kml")[[1]][1]))
  proj<-CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
  eqproj<-CRS(seqinfo$Projection[i])
  eqproj4<-seqinfo$Projection[i]
  
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
  working<-working[as.character(working$Code)%in%as.character(seqinfo$Code[i]),]
  
  ##filter on species where possible
  if (!is.na(seqinfo$Species[i])) {
    working<-working[grep(as.character(seqinfo$Species[i]),as.character(working$Species),ignore.case = T),]
  }
  
  print(working)
  
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
  
  ##calculate corrections for probabilities - because trimmed areas do not give the same probability that a case is in an area as is the case before trimming
  
  ##get raw probabilities
  prob<-working@data$Cases/sum(working@data$Cases)
  ##project generated polygons and KML to equal area projection
  probsfpolygons<-st_make_valid(st_as_sf(spTransform(sppolygons,eqproj)))
  probKML<-st_transform(KML,eqproj4)
  
  ##if not spatial polygons convert
  if (sum(!do.call("rbind",lapply(probKML$geometry, class))[,2]%in%c("POLYGON","MULTIPOLYGON"))!=0) {
    probKML<-st_polygonize(probKML)
    KML<-as(probKML, "Spatial")
    KML<-spTransform(KML,proj)
  }
  
  ##give polygons IDs for tracking purposes
  probsfpolygons$ID<-c(1:nrow(probsfpolygons))
  ##intersect polygons with KML, removing those outside of the region trimming those inside
  probsfpolygonsreduced<-st_intersection(probsfpolygons,probKML)
  ##if no intersecting polygons - break loop 
  if (nrow(probsfpolygonsreduced)==0) {
    sink(file=paste("~/Desktop/TestKMLs/",seqinfo$Sequence[i],".kml",sep=""))
    cat("No cases consistent with this metadata."); cat("\n")
    sink(NULL)
    next()
  }
  ##delete probabilies of polygons removed
  prob<-prob[probsfpolygonsreduced$ID]
  ##for each remaining polygon, test if area differs after intersection with KML and calculate correction factor
  for (u in 1:nrow(probsfpolygonsreduced)) {
    oriarea<-as.numeric(st_area(probsfpolygons[probsfpolygons$ID==probsfpolygonsreduced$ID[u],]))
    newarea<-as.numeric(st_area(probsfpolygonsreduced[u,]))
    if (oriarea>newarea) {
      prob[u]<-prob[u]*(newarea/oriarea)
    }
  }
  ##rescale probabilities post correction
  prob<-prob/sum(prob)
  
  #sppolygons<-gIntersection(sppolygons,KML,byid = T,checkValidity = T)
  sppolygons<-spTransform(as(probsfpolygonsreduced, "Spatial"),proj)
  
  ##initialise pointer and markers
  marker1<-0
  marker2<-0
  u<-1
  
  ##while less than the length of the polygons object
  while(u<length(sppolygons)) {
    #search from the target for intersecting polygons
    for (j in (u+1):length(sppolygons)) {
      if (!is.null(gIntersection(sppolygons[u],sppolygons[j]))) {
        ##if the intersection includes points and lines, as well as a polygon on the basis that two polygons that differ only 
        ##by lines and points are equal from the probability perspective
        if (class(gIntersection(sppolygons[u],sppolygons[j]))=="SpatialCollections") {
          ##check for ring objects if present error out
          if (!is.null(gIntersection(sppolygons[u],sppolygons[j])@polyobj)&!is.null(gIntersection(sppolygons[u],sppolygons[j])@ringobj)) {
            stop("Ring object")
          }
          ##if no ring objects
          if (!is.null(gIntersection(sppolygons[u],sppolygons[j])@polyobj)&is.null(gIntersection(sppolygons[u],sppolygons[j])@ringobj)) {
            ##check area of intersection is not zero for safety
            if (area(gIntersection(sppolygons[u],sppolygons[j])@polyobj)!=0) {
              print(paste(u,j))
              ##see if polygons are internal to one another
              if (!is.null(gDifference(sppolygons[u],sppolygons[j]))) {
                if (area(gDifference(sppolygons[u],sppolygons[j]))!=0) {
                  marker1<-1
                }
              }
              if (!is.null(gDifference(sppolygons[j],sppolygons[u]))) {
                if (area(gDifference(sppolygons[j],sppolygons[u]))!=0) {
                  marker2<-1
                }
              }
              ##if implied equal
              if (marker1==0&marker2==0) {
                ##take the intersection
                temp1<-gIntersection(sppolygons[u],sppolygons[j])@polyobj
                ##record the areas
                area1<-area(sppolygons[u])
                area2<-area(sppolygons[j])
                ##initialise object for new probabilities
                workingprob1<-rep(NA,length(temp1@polygons[[1]]@Polygons))
                
                ##initialise object for new polygons
                bindingpolygons1<-list()
                ##save new polygons and correct probability
                for (q in 1:length(temp1@polygons[[1]]@Polygons)) {
                  workingprob1[q]<-prob[u]*(temp1@polygons[[1]]@Polygons[[q]]@area/area1)+prob[j]*(temp1@polygons[[1]]@Polygons[[q]]@area/area2)
                  bindingpolygons1[[q]]<-Polygons(list(temp1@polygons[[1]]@Polygons[[q]]),q)
                }
                bindingpolygons1<-SpatialPolygons(bindingpolygons1)
                proj4string(bindingpolygons1)<-proj
                
                ##add probability to the list and remove original probabilities
                prob<-c(prob,workingprob1)
                prob<-prob[-c(u,j)]
                
                ##add polygon to the list and remove original polygons
                sppolygons<-c(sppolygons,bindingpolygons1)
                sppolygons<-do.call(bind, sppolygons)
                sppolygons<-SpatialPolygonsDataFrame(sppolygons,data.frame(q=c(1:length(sppolygons)),match.ID = F))
                print(nrow(sppolygons))
                sppolygons<-sppolygons[-c(u,j),]
                print(nrow(sppolygons))
                sppolygons<-SpatialPolygons(sppolygons@polygons,proj4string=sppolygons@proj4string)
                
                ##reset pointer as polygon at i now different
                u<-u-1
              }
              
              ##polygon validity checking
              
              ##find any polygons with invalid geometries
              valid<-rep(T,length(sppolygons))
              for (m in 1:length(sppolygons)) {
                if (!gIsValid(sppolygons[m])) {
                  ##if invalid, attempt to fix by using st_make_valid
                  sppolygons@polygons[m]<-as(st_make_valid(st_as_sf(sppolygons[m])), "Spatial")@polygons
                  ##if still invalid enter error routine
                  if (!gIsValid(sppolygons[m])) {
                    valid[m]<-F
                    print("Invalid geometry detected")
                  }
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
        if (class(gIntersection(sppolygons[u],sppolygons[j]))=="SpatialPolygons") {
          ##check area of intersection is not zero for safety
          if (area(gIntersection(sppolygons[u],sppolygons[j]))!=0) {
            print(paste(u,j))
            ##see if polygons are internal to one another
            if (!is.null(gDifference(sppolygons[u],sppolygons[j]))) {
              if (area(gDifference(sppolygons[u],sppolygons[j]))!=0) {
                marker1<-1
              }
            }  
            if (!is.null(gDifference(sppolygons[j],sppolygons[u]))) {
              if (area(gDifference(sppolygons[j],sppolygons[u]))!=0) {
                marker2<-1
              }
            }
            ##if implied equal
            if (marker1==0&marker2==0) {
              ##take the intersection
              temp1<-gIntersection(sppolygons[u],sppolygons[j])
              ##record the areas
              area1<-area(sppolygons[u])
              area2<-area(sppolygons[j])
              ##initialise object for new probabilities
              workingprob1<-rep(NA,length(temp1@polygons[[1]]@Polygons))
              
              ##initialise object for new polygons
              bindingpolygons1<-list()
              ##save new polygons and correct probability
              for (q in 1:length(temp1@polygons[[1]]@Polygons)) {
                workingprob1[q]<-prob[u]*(temp1@polygons[[1]]@Polygons[[q]]@area/area1)+prob[j]*(temp1@polygons[[1]]@Polygons[[q]]@area/area2)
                bindingpolygons1[[q]]<-Polygons(list(temp1@polygons[[1]]@Polygons[[q]]),q)
              }
              bindingpolygons1<-SpatialPolygons(bindingpolygons1)
              proj4string(bindingpolygons1)<-proj
              
              ##add probability to the list and remove original probabilities
              prob<-c(prob,workingprob1)
              prob<-prob[-c(u,j)]
              
              ##add polygon to the list and remove original polygons
              sppolygons<-c(sppolygons,bindingpolygons1)
              sppolygons<-do.call(bind, sppolygons)
              sppolygons<-SpatialPolygonsDataFrame(sppolygons,data.frame(q=c(1:length(sppolygons)),match.ID = F))
              print(nrow(sppolygons))
              sppolygons<-sppolygons[-c(u,j),]
              print(nrow(sppolygons))
              sppolygons<-SpatialPolygons(sppolygons@polygons,proj4string=sppolygons@proj4string)
              
              ##reset pointer as polygon at i now different
              u<-u-1
            }
            
            ##polygon validity checking
            
            ##find any polygons with invalid geometries
            valid<-rep(T,length(sppolygons))
            for (m in 1:length(sppolygons)) {
              if (!gIsValid(sppolygons[m])) {
                ##if invalid, attempt to fix by using st_make_valid
                sppolygons@polygons[m]<-as(st_make_valid(st_as_sf(sppolygons[m])), "Spatial")@polygons
                ##if still invalid enter error routine
                if (!gIsValid(sppolygons[m])) {
                  valid[m]<-F
                  print("Invalid geometry detected")
                }
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
    ##reset markers and update pointer
    marker1<-0
    marker2<-0
    u<-u+1
  }
  
  ##change to equal area CRS - depends on input region
  sppolygons<-spTransform(sppolygons,eqproj)
  
  #convert to sf object
  sfpolygons<-st_make_valid(st_as_sf(sppolygons))
  
  ##take intersection of all polygons with themselves and filter out results with 0 area
  sfintersection<-st_intersection(sfpolygons)
  sfintersection<-sfintersection[as.numeric(st_area(sfintersection))!=0,]
  sfintersection<-recursivesimplify(sfintersection)
  
  ##initialise object for final probabilities
  newprob<-rep(NA,nrow(sfintersection))
  
  ##record areas for probability calculation
  areasoriginal<-as.numeric(st_area(sfpolygons))
  areasfinal<-as.numeric(st_area(sfintersection))
  
  ##calculate correct final probabilities
  for (h in 1:nrow(sfintersection)) {
    origins<-sfintersection$origins[h][[1]]
    newprob[h]<-sum(prob[origins]*areasfinal[h]/areasoriginal[origins])
  }
  
  sfintersection$Prob<-newprob
  
  sfintersection<-st_transform(sfintersection,4326,check=T)[,c(3,4)]
  
  if (lostprobability>0) {
    stop(paste("Lost probability at iteration ",i,".",sep = ""))
  }
  ##generate KML - modification of Simon's code
  
  sink(file=paste("~/Desktop/TestKMLs/",seqinfo$Sequence[i],".kml",sep=""))
  cat("<?xml version=\"1.0\" encoding=\"UTF-8\"?>"); cat("\n")
  cat("<kml xmlns=\"http://earth.google.com/kml/2.2\">"); cat("\n")
  for (i in 1:nrow(sfintersection)) {
    cat(paste("\t<polygon id=\"",paste(seqinfo$Sequence[i],i,sep="_"),"\" samplingProbability=\"",sfintersection$Prob[i],"\">",sep="")); cat("\n")
    cat("\t\t<coordinates>"); cat("\n")
    for (j in 1:nrow(sfintersection[i,1]$geometry[[1]][1][[1]])) {
      cat(paste("\t\t\t",sfintersection[i,1]$geometry[[1]][1][[1]][j,2],",",sfintersection[i,1]$geometry[[1]][1][[1]][j,1],",0",sep="")); cat("\n")
    }
    cat("\t\t</coordinates>"); cat("\n")
    cat("\t</polygon>"); cat("\n")
  }
  cat("</kml>"); cat("\n")
  sink(NULL)
  } else {
    sink(file=paste("~/Desktop/TestKMLs/",seqinfo$Sequence[i],".kml",sep=""))
    cat("No cases consistent with this metadata."); cat("\n")
    sink(NULL)
  }
}
