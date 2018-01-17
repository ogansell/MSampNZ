library(rgdal)
library(mapview)
library(geosphere)
library(sp)

#Function for generating transects from a start point provided by the master sample
#' @name maketransect
#' @title Make transects of a specified distance from SpatialPointsDataFrames
#' @export
maketransect <-function(x,tran.length){

  nztm <-"+proj=tmerc +lat_0=0 +lon_0=173 +k=0.9996 +x_0=1600000 +y_0=10000000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"
  wgs84 <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")#Set CRS to convert to WGS84

df1 <- x

df1$deg <-sample(1:360,nrow(df1), replace = TRUE)#Random bearing for transect

proj4string(df1) = CRS(nztm)
end.pt<- as.data.frame(destPoint(coord(df1), b=df1$deg, d=tran.length))#Calculate endpoint of transect given mastersample point is start point
end.pt.sp<-SpatialPoints(end.pt[,c(1,2)])
end.pt <-SpatialPointsDataFrame(end.pt.sp,end.pt)
proj4string(end.pt) <- CRS(wgs84) # specify CRS
end.pt <-as.data.frame(end.pt)
end.pt <-end.pt[,c(1:2)] #Extract coords
df1 <- spTransform(df1, wgs84)
st.pt <-as.data.frame(df1) #Make a dataframe of start point
st.pt <-st.pt[,c(4:5)] #Extract coords
names(st.pt)<-c("lon","lat") #Fix names

degrees <-round(as.data.frame(finalBearing(st.pt,end.pt, sphere = TRUE))+180,0)#Calculate true bearing between start and end point
names(degrees)<-"degrees"
degrees$ID <-1:nrow(df1)

## Create a raw list to store Lines objects
l <- vector("list", nrow(st.pt))
library(sp)
for (k in seq_along(l)) {
  l[[k]] <- Lines(list(Line(rbind(st.pt[k,], end.pt[k,]))), as.character(k))
}
#Create SpatialLines object
sline<-SpatialLines(l)
sline.df = SpatialLinesDataFrame(sline,data.frame(ID = 1:nrow(df1)))

df1$ID <-1:nrow(df1) #Create ID to merge on
df1.lines <-merge(sline.df,as.data.frame(df1[,c(1:2,4)]),by = "ID")
df1.lines <-merge(df1.lines,degrees,by = "ID")
df1.lines <-merge(df1.lines,as.data.frame(x),by = c("SiteID","SiteOrder"))
# df1.lines@data$coords.x1.x<-NULL;df1.lines@data$coords.x2.x<-NULL;df1.lines@data$ID<-NULL
# names(df1.lines@data)[4:5]<-c("EastingNZTM","NorthingNZTM")
# df1.lines@data$EastingNZTM<-round(df1.lines@data$EastingNZTM,0);df1.lines@data$NorthingNZTM<-round(df1.lines@data$NorthingNZTM,0)

proj4string(df1.lines) <- CRS(wgs84) # specify CRS
df1.lines <- spTransform(df1.lines, nztm)#transform to NZTM
}
