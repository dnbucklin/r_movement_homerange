# This script takes a table of points and outputs a points shapefile with information about the trajectory steps (angles, distance, time)

#######To run this code, you need to have a data file (in text format) with the following columns (using the exact names below (all uppercase)):
## ID - unique animal ID
## DATE - day of year in mm/dd/yyyy format
## TIME - minute of day (this can range from 0 - 1440)
## POINT_X - x spatial coordinate
## POINT_Y - y spatial coordiante
#######

#change table input and output files, and then run rest of the code
input_table<-"N:/Kemps_ridley_internesting/3_randwalk_ud_allproj_6_2014/all_turtles_6_2014.txt"
output_shapefile<-"C:/folder/test2"
##########

#load libraries
library(adehabitatLT)
library(maptools)

#read table
dat<-read.table(input_table,header=T)

#add time
hour<-trunc(dat$TIME/60)
min<-trunc(((dat$TIME/60)-hour)*60)
time<-as.POSIXct(strptime(paste0(dat$DATE," ",hour,":",min),"%m/%d/%Y %H:%M"))
dat$fulltime<-time

#remove duplicates by ID,time
dat<-dat[!duplicated(dat[c("ID","fulltime")]),] 

#extract coordinates
coords<-coordinates(cbind(dat$POINT_X,dat$POINT_Y))

#convert to trajectory
lt<-as.ltraj(coords,date=dat$fulltime,id=dat$ID)

#output data
lt2<-ltraj2spdf(lt)
lt2@data$date<-as.character(lt2@data$date)

#add IDs
lt2@data$pkey<-as.character(lt2@data$pkey)
lt2@data$id<-unlist(lapply(lt2@data$pkey,function(x) {strsplit(x,".",fixed=T)[[1]][1]}))

#convert to degrees
lt2@data$a_ang_deg<- (lt2@data$abs.angle*180) / (pi)
lt2@data$r_ang_deg<- (lt2@data$rel.angle*180) / (pi)

#output shapefile
writePointsShape(lt2,output_shapefile)


