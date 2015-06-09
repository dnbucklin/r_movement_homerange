#R code for Utilization distributions (MCP, KDE)

#load packages (make sure they are installed first)
library(ade4)
library(adehabitatHR)
library(maptools)
library(raster)
library(rgeos)

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

#load data table - this file should be in your working directory
data1<-read.table("file.txt",header=T)                          #######CHANGE THIS#########

#standard deviation ratio settings for rescaling KDEs
sd.ratio.min<-0.5     #######CHANGE THIS - OPTIONAL#########
sd.ratio.max<-1.5     #######CHANGE THIS - OPTIONAL#########

#MCP percentage (0 - 100)
mcp.per<-100			#######CHANGE THIS - OPTIONAL#########
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
    
    #calcuate mean daily locations
    mdl<-aggregate(cbind(ani$POINT_X,ani$POINT_Y),by=list(ani$DATE),FUN="mean")
    mdl.out<-cbind(uniqid,mdl)
    write.table(mdl.out,file="mean_daily_locs_all.txt",append=T,col.names=F,row.names=F)
    print(paste0(length(mdl[,1])," mean daily locations."))
    
    #make spatial data frame with mean daily locations
    locs<-SpatialPoints(data.frame(x=mdl$V1,y=mdl$V2))
    alllocs<-SpatialPoints(data.frame(x=ani$POINT_X,y=ani$POINT_Y))
    
    if (length(alllocs$x)>4){
      
      #MCP (uses all locations)
      cp<-mcp(alllocs,percent=mcp.per)
      cp$id<-uniqid
      writePolyShape(cp, paste0("utilization_distributions/",uniqid,"_mcp"))
      
      if(length(locs$x)>19){
        #KDE (uses mean daily locations)
        #calculate x/y standard deviation ratio for KDE
        sd.x<-sd(locs$x)
        sd.y<-sd(locs$y)
        rat.xy<-sd.x/sd.y
        print(rat.xy)
        
        #if ratio is >0.5 or <1.5, use regular locations, otherwise divide by standard devation
        if (rat.xy > sd.ratio.min & rat.xy < sd.ratio.max)  {locs.kde<-locs
        print("Using regular coordinates")
        
        #KDE *what grid number to choose (grid ?) *this is number of cells for longest direction*
        kde<-kernelUD(locs.kde,h = "LSCV",grid=1000)
        h.val<-c(uniqid,length(locs$x),rat.xy,kde@h$convergence,kde@h$h)
        print(h.val)
        vud <- getvolumeUD(kde)
        
        #export 50/95% UDs
        hr50 <- getverticeshr(kde,percent=50,ida="kde50")
        hr50$id<-uniqid
        hr95 <- getverticeshr(kde,percent=95,ida="kde95")
        hr95$id<-uniqid
        writePolyShape(spRbind(hr95,hr50), paste0("utilization_distributions/",uniqid,"_kde"))
        capture.output(h.val,file="utilization_distributions/ud_output.txt",append=T)}
        
        else {locs.kde<-SpatialPoints(data.frame(x=locs$x/sd.x,y=locs$y/sd.y)) 
        print("Using transformed coordinates")
        
        kde<-kernelUD(locs.kde,h = "LSCV",grid=1000)
        h.val<-c(uniqid,length(locs$x),rat.xy,kde@h$convergence,kde@h$h)
        print(h.val)
        vud <- getvolumeUD(kde)
        
        #get kernel density volume 
        fud <- vud[[1]]
        hr<-as.data.frame(fud)[,1]
        
        hr2<-data.frame(hr)
        ka<-data.frame(x=(coordinates(vud)[,1])*sd.x,y=(coordinates(vud)[,2])*sd.y,z=hr2)
        kb<-rasterFromXYZ(ka)
        kc<-rasterToContour(kb,maxpixels=3000000,levels=c(50,95))
        kc$id<-uniqid
        writeLinesShape(kc, paste0("utilization_distributions/",uniqid,"_kde"))
        
        print(u)
        capture.output(h.val,file="utilization_distributions/ud_output.txt",append=T)
        } 
      } else {
        capture.output(c(uniqid,length(alllocs$x),length(locs$x)),file="utilization_distributions/ud_output.txt",append=T)
        print("Less than 20 mean daily locations, no kernel UD calcualated.")}
      
    } else {
      capture.output(c(uniqid,length(alllocs$x),length(locs$x)),file="utilization_distributions/ud_output.txt",append=T)
      print("Less than 5 total locations, no MCP calculated.")}
  }
  
}

save.image("utilization_distributions/ud_rWorkspace.Rdata")
