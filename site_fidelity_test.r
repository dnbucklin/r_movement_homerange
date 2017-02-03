## Site fidelity tests
## 2/3/2017

#######To run this code, you need to have a data file (in text format) with at least the following columns (using the exact names below):
## ID - unique animal ID
## TID - unique track ID (to differentiate between multiple tracks for the same animal - use the start date of the unique track)
## DATE - day of year in mm/dd/yyyy format
## TIME - minute of day (this can range from 0 - 1440)
## POINT_X - x spatial coordinate
## POINT_Y - y spatial coordiante
#######

## All animals in one text file, they must share the same projection

#load libraries
library(adehabitatLT)
library(ade4)
library(raster)
library(rgeos)
library(stats)
library(maptools)

############################################
####CHANGE SETTINGS IN THIS SECTION#########
#####RUN THIS SECTION FIRST#################

## working directory - make a new folder containing the input file and constraining shapefile. This is also where results are generated. Set below.
setwd("C:/folder")          #######CHANGE THIS#########

## load data table - this file should be in your working directory
data1<-read.table("file.txt",header=T)               #######CHANGE THIS#########

## How many random walks per animal track? Set below
rw.num<-100                                               #######CHANGE THIS#########

## Plot limit distance (in projection units - set to ~1 if using lat/lon coordinates, set to ~100000 if using meters)
d<-100000       #######CHANGE THIS#########

## can random walk trajectories cross outside constraining polygon? (points making up trajectory are not allowed outside)
## not allowing lines outside (FALSE) will cause tests to run slower
cross<-FALSE

## Constraint function shapefile (to restrict random walks from going outside certain areas)
## The shapefile should indicate where the animal CAN go, and all points MUST fall inside it (code checks this below)
## It must be in the same projection as your points
par<-readShapePoly("shapefile")     #######CHANGE THIS (last one)#########

############################################
#### END CHANGE SETTINGS SECTION ###########
############################################

##########################################################################
### THIS SECTION CHECKS IF POINTS OVERLAP CONSTRAINT POLYGON #############
##########################################################################

# add attribute to par
par@data<-data.frame(x = 1)

## run the checks below (until END CHANGE SETTINGS SECTION) prior to running rest of code
check<-over(SpatialPoints(cbind(data1$POINT_X,data1$POINT_Y)),par)
check.2 <- all(!is.na(check))
if (!check.2) {print(paste0("ERROR: ",sum(is.na(check))," of ",length(check[,1])," points in the tracks fall outside the bounding polygon. Remove them before running this code"))} else {print("Data are ready")}

plot(SpatialPoints(cbind(data1$POINT_X,data1$POINT_Y)))
plot(par,add=TRUE)

##########################################################################
#### END CHECK SECTION ###################################################
##########################################################################
## RUN REST OF CODE IF THERE ARE NO ERRORS ###############################

dir.create("site_fidelity")

col.names<-c("ani_id","num_relocations","r2","mean_r2_RW","r2_pval","linearity","mean_lin_RW","lin_pval")
capture.output(col.names,file="site_fidelity/random_walk_results.txt",append=T)

##constraint function for limiting area of movement on random walks
confun <- function(x, par){
  ## Define a SpatialPointsDataFrame from the trajectory
  coordinates(x) <- x[,1:2]
  
  ## overlap the relocations x to the map par (use for points)
  jo<-over(x,par)
  ## checks that there are no missing value
  res <- all(!is.na(jo))
  
  if (!cross & res) {
    #define spatiallines (use for lines)
    spl<-as(x,"SpatialLines")
    res <- gContains(par, spl) # gContains allows spl to intersect boundary of par
    # gContainsProperly would disallow boundary touching
  }
  
  ## return this check
  return(res)
}

### IF loop gets stuck and is taking too long for one animal (probably due to the number of points and/or the size of constraint polygon), hit escape to end process.
### You can restart the loop on the next animal by running all code below with the subset not yet run - set this by modifying the
### for (ani in list.ani){ ### line (5 lines down) to include just that subset 
### - for example, ## for (ani in list.ani[37:43]){ ## would run for animals 37-43 in the list.

#loop by animal ID
list.ani<-unique(data1$ID)
for (ani in list.ani){
  
  print(paste0("Animal #",which(list.ani==ani)," of ",length(list.ani)))
  data<-subset(data1,data1$ID==as.character(ani))
  print(ani)
  #date/time normalization - format is mm/dd/yyyy in DATE column, minute of day in TIME column
  da<-as.POSIXct(strptime(paste0(data$DATE," ",trunc(data$TIME/60)," ",(data$TIME/60-(trunc(data$TIME/60)))*60),"%m/%d/%Y %H %M"))
  data<-cbind(data,da)
  data<-data[!duplicated(data[,c("ID","da")]),]
  turt <- as.ltraj(xy = data[,c("POINT_X","POINT_Y")], date = data$da, id = as.character(data$ID),burst=paste0(data$ID,".",data$TID))
  
  #set up correlated random walks, with randomized angles, randomized distances, fixed starting point
  mos<-NMs.randomCRW(turt, rangles = TRUE, rdist = FALSE, fixedStart = TRUE,
                     x0 = NULL, rx = NULL, ry = NULL, treatment.func = NULL,
                     treatment.par = NULL, constraint.func = confun,
                     constraint.par = par, nrep = rw.num)
  
  ############### This step can take a lot of time if the track has many points
  
  ## Then apply the null model
  set.seed(125)   ## to make the calculation reproducible
  rand <- testNM(mos,count=TRUE)
  
  ##compare r2 and linearity of true path to random walks
  
  for (an in 1:length(turt)){
    
    #calculate true path r2 and linearity
    animal<-gsub("/","-",burst(turt)[an])
    
    track.t<-turt[[an]]
    relocs<-length(track.t[,1])
    
    #make plot
    png(file = paste0("site_fidelity/",animal,"_plot.png"),pointsize=12,width=1000,height=1000)
    
    #xlim<-c(bbox(par)[1,1],bbox(par)[1,2])
    #ylim<-c(bbox(par)[2,1],bbox(par)[2,2])
    
    xlim<-c(min(turt[[an]]$x)-d,max(turt[[an]]$x)+d)
    ylim<-c(min(turt[[an]]$y)-d,max(turt[[an]]$y)+d)
    
    plot(turt[an],spoldf=par, xlim=xlim, ylim=ylim,main=paste0(ani," tracks"))
    
    meanx.t<-mean(track.t$x)
    meany.t<-mean(track.t$y)
    r2.true<-(sum((track.t$x-meanx.t)^2)/length(track.t$x))+(sum((track.t$y-meany.t)^2)/length(track.t$y)) #r-squared of true track (Schoener 1981)
    
    leng<-sum(turt[[an]]$dist,na.rm=T)
    start<-track.t[1,][1:2]
    end<-track.t[length(track.t[,1]),][1:2]
    lin.true<-(dist(rbind(start,end)))/leng    #linearity of true track (distance from start to finish / total distance)
    foo.r2<-NULL
    foo.lin<-NULL
    
    #calculate r2 and linearity distributions for random walks
    for (i in 1:rw.num){
      
      track<-rand[[an]][[i]]
      meanx<-mean(track$x)
      meany<-mean(track$y)
      r2<-(sum((track$x-meanx)^2)/length(track$x))+(sum((track$y-meany)^2)/length(track$y)) #r-squared of simulated track (Schoener 1981)
      
      start<-track[1,][1:2]
      end<-track[length(track[,1]),][1:2]
      lin<-(dist(rbind(start,end)))/leng    #linearity of simulated track (distance from start to finish / total distance)
      lines(rand[[an]][[i]]$x,rand[[an]][[i]]$y,col="grey60") #plots the random walk line on the graph
      
      foo.r2[i]<-r2
      foo.lin[i]<-lin
    }
    
    #statistical tests (proportion of random r2 and linearity values less than true value)
    print(paste0("Animal ID: ",animal))
    
    mean.r2<-mean(foo.r2)
    mc.r2<-as.randtest(foo.r2,r2.true,alter="greater") #monte carlo test
    #print(mc.r2)
    
    mean.lin<-mean(foo.lin)
    mc.lin<-as.randtest(foo.lin,lin.true,alter="greater") #monte carlo test
    #print(mc.lin)
    
    #save result p-values
    vect<-c(animal,relocs,round(c(r2.true,mean.r2,mc.r2$pvalue,lin.true,mean.lin,mc.lin$pvalue),7))
    capture.output(vect,file="site_fidelity/random_walk_results.txt",append=T)
    lines(track.t$x,track.t$y,col="black", lwd = 3)
    
    dev.off()
    plot.new()
  }
  
}

save.image("site_fidelity/site_fidelity_rWorkspace.Rdata")
