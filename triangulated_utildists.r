#R code for Utilization distributions (alpha hull, alpha shape, characteristic hull)
# make sure to put 'alpha-functions.R' file in your working directory

#load packages (make sure they are installed first)
library(adehabitatHR)
library(maptools)
library(ade4)
library(raster)
library(rgeos)
library(alphahull)

#######To run this code, you need to have a data file (in text format) with the following columns (using the exact names below (all uppercase)):
## ID - unique animal ID
## TID - unique track ID (to differentiate between multiple tracks for the same animal (for example, you can use the start date of the unique track as the TID))
## DATE - day of year in mm/dd/yyyy format
## TIME - minute of day (this can range from 0 - 1440)
## POINT_X - x spatial coordinate
## POINT_Y - y spatial coordiante
#######

######CHANGE SETTINGS SECTION########
#working directory - folder containing input data and where to export results. Set below
setwd("/directory")           #######CHANGE THIS#########
source("alpha-functions.R")

#load data table - this file should be in your working directory
data1<-read.table("file.txt",header=T)                          #######CHANGE THIS#########

######END CHANGE SETTINGS SECTION########

#######RUN rest of code########
dir.create("utilization_distributions")

head<-c("unique_id","num_total_locs","num_meandailylocs","ratio_xy","lscv_converge","lscv_h")
capture.output(head,file="utilization_distributions/ud_output.txt",append=T)

list.ani<-unique(data1$ID)
print(list.ani)

for (idnum in list.ani){
  
  #subset data for unique animal
  data<-subset(data1,data1$ID==as.character(idnum))
  
  for (u in unique(data$TID)){
    
    #subset data for unique animal track
    ani<-subset(data,data$TID==u)
    uniqid<-paste0(idnum,"_",gsub("/","_",u))
    print(uniqid)
    
    #make spatial data frame with all locations
    alllocs<-SpatialPoints(data.frame(x=ani$POINT_X,y=ani$POINT_Y))
    
    if (length(alllocs$x)>4){
      
      #Alpha hull
      #dv<-delvor(coordinates(alllocs))$mesh
      r <- tri.mesh(alllocs)
      k <- seq_len( r$tlnew - 1 )
      i <- r$tlist[k]          
      j <- r$tlist[r$tlptr[k]]
      keep <- i > 0
      i <- abs( i[ keep ] )
      j <- abs( j[ keep ] )
      plot(alllocs)
      segments( r$x[i], r$y[i], r$x[j], r$y[j], col="grey" )
      distances <- sqrt( ( r$x[i] - r$x[j] ) ^ 2 + ( r$y[i] - r$y[j] ) ^ 2 )
      m<-mean(distances)
      sdm<-sd(distances)/sqrt(length(distances))
      #defaults to 3 times standard deviation of the mean
      a.hull<-ahull(coordinates(alllocs),alpha=m+(3*sdm))
      a.shape<-ashape(coordinates(alllocs),alpha=m+(3*sdm))
      
      ah<-ahull_to_SPLDF(a.hull)
      ah$id<-uniqid
      as<-ashape_to_SPLDF(a.shape)
      as$id<-uniqid
      writeLinesShape(ah, paste0("utilization_distributions/",uniqid,"_ahull"))
      writeLinesShape(as, paste0("utilization_distributions/",uniqid,"_ashape"))
      
      ##charHull (Downs and Horner, 2009)
      ch<-CharHull(alllocs)
      numtris<-round(length(ch)*0.95,0)
      IDs<-rep(1,numtris)
      ch95a<-spCbind(ch[1:numtris,],IDs)
      area_m2<-sum(ch95a$area)
      ch95b<-unionSpatialPolygons(ch95a,IDs=IDs)
      ch95c<-SpatialPolygonsDataFrame(ch95b,data=as.data.frame(area_m2))
      ch95c$id<-uniqid
      writePolyShape(ch95c,paste0("utilization_distributions/",uniqid,"_charhull95"))
      
    } else {
      capture.output(c(uniqid,length(alllocs$x),file="utilization_distributions/ud_output.txt",append=T)
      print("Less than 5 total locations, no alpha shapes calculated.")}
  }
  
}

save.image("utilization_distributions/ud_rWorkspace.Rdata")
