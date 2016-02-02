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
setwd("C:/Users/dbucklin/Desktop/test")           #######CHANGE THIS#########

#load data table - this file should be in your working directory
data1<-read.table("all_3_red.txt",header=T)                          #######CHANGE THIS#########

#standard deviation ratio settings for rescaling KDEs
sd.ratio.min<-0.5     #######CHANGE THIS - OPTIONAL#########
sd.ratio.max<-2       #######CHANGE THIS - OPTIONAL#########

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
        h.val<-c(uniqid,length(alllocs$x),length(locs$x),round(rat.xy,3),kde@h$convergence,round(kde@h$h,4))
        print(h.val)
        vud <- getvolumeUD(kde)
        
        #get kernel density volume 
        fud <- vud[[1]]
        hr<-as.data.frame(fud)[,1]
        
        hr2<-data.frame(hr)
        ka<-data.frame(x=(coordinates(vud)[,1]),y=(coordinates(vud)[,2]),z=hr2)
        kb<-rasterFromXYZ(ka,digits=8)
        kc<-rasterToContour(kb,maxpixels=3000000,levels=c(50,95))
        kc$id<-uniqid
      
        #writeLinesShape(kc, paste0("utilization_distributions/",uniqid,"_kde"))
        #for KDE line output, uncomment above
        capture.output(h.val,file="utilization_distributions/ud_output.txt",append=T)}
        
        else {locs.kde<-SpatialPoints(data.frame(x=locs$x/sd.x,y=locs$y/sd.y)) 
        print("Using transformed coordinates")
        
        kde<-kernelUD(locs.kde,h = "LSCV",grid=1000)
        h.val<-c(uniqid,length(alllocs$x),length(locs$x),round(rat.xy,3),kde@h$convergence,round(kde@h$h,4))
        print(h.val)
        vud <- getvolumeUD(kde)
        
        #get kernel density volume 
        fud <- vud[[1]]
        hr<-as.data.frame(fud)[,1]
        
        hr2<-data.frame(hr)
        ka<-data.frame(x=(coordinates(vud)[,1])*sd.x,y=(coordinates(vud)[,2])*sd.y,z=hr2)
        kb<-rasterFromXYZ(ka,digits=8)
        kc<-rasterToContour(kb,maxpixels=3000000,levels=c(50,95))
        kc$id<-uniqid
        
        #writeLinesShape(kc, paste0("utilization_distributions/",uniqid,"_kde"))
        #for KDEline output, uncomment above
        print(u)
        capture.output(h.val,file="utilization_distributions/ud_output.txt",append=T)
        } 
        
        #begin hole removal
        ct<-0
        vertsort<-c(95,50)
        #loop through each KDE level
        for (lev in vertsort) {
          ct<-ct+1
          t<-subset(kc,kc@data$level==lev)
          
          #convert to poly
          t1<-SpatialPolygonsDataFrame(gPolygonize(t),data=data.frame(level=seq(1,length(gPolygonize(t)),1)))
          t2<-unionSpatialPolygons(t1,IDs=rep(as.character(ct),length(t1)))
          
          holetab<-data.frame(id=NA,index=NA)
          
          for (i in t1@data$level)
          {
            l<-length(t1@polygons[[i]]@Polygons)
            if (l == 1) {next} else {
              for (p in 1:length(t1@polygons[[i]]@Polygons)) {
                polys<-t1@polygons[[i]]@Polygons
                poly<-t1@polygons[[i]]@Polygons[[p]]
                if (isTRUE(poly@hole)) {
                  holetab<-rbind(holetab,c(i,p))
                }
              }
            }
          }
          
          holetab<-holetab[complete.cases(holetab),]
          if (length(holetab[,1]) == 0) {
            out<-SpatialPolygonsDataFrame(t2,data=data.frame(row.names=ct,level=lev))
            assign(paste0('kernel',lev),out)
            next
          }
          
          for (m in 1:length(holetab[,1]))
          {
            ind<-holetab[m,]
            hole<-t1@polygons[[ind$id]]@Polygons[[ind$index]]
            hole@hole<-FALSE
            hole<-Polygons(list(hole),'1')
            hole2<-SpatialPolygons(list(hole))
            t2<-gDifference(t2,hole2,id=as.character(ct))
          }
          
          out<-SpatialPolygonsDataFrame(t2,data=data.frame(row.names=ct,level=lev))
          assign(paste0('kernel',lev),out)
        }
        
        #browser()
        if (length(vertsort) > 1) {
          str<-paste('kde.out<-rbind(',paste(paste('kernel',vertsort,sep=""),collapse=','),')',sep="") 
        } else {
          str<-paste('kde.out<-',paste('kernel',vertsort,sep=""),sep="")
        }
        
        eval(parse(text=str))
        #end hole removal
        
        writePolyShape(kde.out, paste0("utilization_distributions/",uniqid,"_kdepoly"))
        
      } else {
        capture.output(c(uniqid,length(alllocs$x),length(locs$x)),file="utilization_distributions/ud_output.txt",append=T)
        print("Less than 20 mean daily locations, no kernel UD calcualated.")}
      
    } else {
      capture.output(c(uniqid,length(alllocs$x),length(locs$x)),file="utilization_distributions/ud_output.txt",append=T)
      print("Less than 5 total locations, no MCP calculated.")}
  }
  
}

save.image("utilization_distributions/ud_rWorkspace.Rdata")