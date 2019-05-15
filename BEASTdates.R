###code for reading in OIE data on a per sequence level and generating 
library(lubridate)

##read in sequence IDs, KMLs and dates
seqinfo<-read.csv("~/Dropbox/BTV-8/SamplesReducedDates.csv",stringsAsFactors=F)

decdates<-rep(NA,nrow(seqinfo))

##iterate through sequences
for (i in 1:nrow(seqinfo)) {
  ##filter on dates
  if (is.na(seqinfo$Month[i])) {
    decdates[i]<-as.numeric(seqinfo$Year[i])+0.5
  } else if (is.na(seqinfo$Day[i])) {
    if (nchar(seqinfo$Month[i])==1) {
      seqinfo$Month[i]<-paste0(0,seqinfo$Month[i])
    }
    if (as.numeric(seqinfo$Month[i])%in%c(4,6,9,11)) {
      decdates[i]<-(decimal_date(as.Date(paste(seqinfo$Year[i],"-",seqinfo$Month[i],"-15",sep="")))
                 +decimal_date(as.Date(paste(seqinfo$Year[i],"-",seqinfo$Month[i],"-16",sep=""))))/2
    } else if (as.numeric(seqinfo$Month[i])%in%c(2)) {
      if (leap_year(decimal_date(as.Date(paste(seqinfo$Year[i],"-",seqinfo$Month[i],"-01",sep=""))))==T) {
        decdates[i]<-decimal_date(as.Date(paste(seqinfo$Year[i],"-",seqinfo$Month[i],"-15",sep="")))
      } else {
        decdates[i]<-(decimal_date(as.Date(paste(seqinfo$Year[i],"-",seqinfo$Month[i],"-15",sep="")))
                      +decimal_date(as.Date(paste(seqinfo$Year[i],"-",seqinfo$Month[i],"-14",sep=""))))/2
      }
    } else {
      decdates[i]<-decimal_date(as.Date(paste(seqinfo$Year[i],"-",seqinfo$Month[i],"-16",sep="")))
    }
  } else {
    if (nchar(seqinfo$Day[i])==1) {
      seqinfo$Day[i]<-paste0(0,seqinfo$Day[i])
    }
    decdates[i]<-decimal_date(as.Date(paste(seqinfo$Year[i],"-",seqinfo$Month[i],"-",seqinfo$Day[i],sep="")))
  }
}

dates<-data.frame(seqinfo$Sequence,decdates)

write.table(dates,"~/Desktop/sequencedates.txt",col.names = FALSE,row.names = F,sep="\t")
