########################
#Master Sample Code
#Paul van Dam-Bates
########################
#Take in a polygon and generate points from the master sample
#Step 1: Determine Island that the polygon falls in
#Step 2: Use random seed from that island to start Halton Sequence
#Step 3: Ouptut number of points required clipped for that region.
#--------------------------------------------------------------------

library(sp)
library(raster)
library(rgeos)
library(sp)

#' @export
#Halton Sequence:
RSHalton <- function(n = 10,seeds = c(0,0),bases = c(2,3)) {
  ##
  ## Generate n points from a random start d dimensional Halton sequence.
  ##
  ## Inputs:
  ##
  ## n 			sample size
  ## bases    	coprime bases e.g. c(2,3) (Halton Sequence)
  ## seeds  		random seeds  e.g. c(0,0) (Halton Sequence)


  ########### Initialize #########################################
  d <- length(bases);	pts <- mat.or.vec(n, d)
  if (length(seeds) != d){
    seeds <- rep(seeds[1],d)
  }

  ########### Main Loop #########################################
  for (i in 1:d) {
    b <- bases[i];   	u <- seeds[i];	k <- u:(u+n-1);
    xk <- (k %% b)/b;
    for (j in 1:(ceiling(logb(u+n,b)) + 2)) {
      xk <- xk + (floor(k/(b^j)) %% b)/(b^(j+1));
    }
    pts[,i] <- cbind(xk)
  }
  pts <- cbind(1:nrow(pts), pts)
  return(pts)
}

#' @export
masterSample <- function(island = "South", shp, N = 100){
  #Master Sample seed for South Island, chosen as first random start that fell into SI
  #seed.si <-  c(4887260, 18041662)
  #seed.ni <- c(5137598, 8906854)

  #Define CRS
  nztm <-"+proj=tmerc +lat_0=0 +lon_0=173 +k=0.9996 +x_0=1600000 +y_0=10000000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"

  #If CRS does not match fix that!
  if(nztm != proj4string(shp)){
    cat("ERROR: CRS does not match, please project for:", nztm, "\n")
    return()
  }

  if(island == "South")
  {
    bb <- data.frame(min = c(1089354,4747979), max = c(1721164,5516919), row.names = c("x","y"))
    seed <-  c(4887260, 18041662)
  }else if(island == "North")
  {
    bb <- data.frame(min = c(1510593,5390569), max = c(2092000,6223164), row.names = c("x","y"))
    seed <- c(5137598, 8906854)
  }else{
    cat("ERROR: Define Island for MS \n")
    return()
  }

  #Scale and shift Halton to fit into bounding box
  scale.bas <- bb[,2] - bb[,1]
  shift.bas <- bb[,1]

  #Area proportional sampling we can use this to make this code quicker.
  prop.area <- sum(area(shp))/(scale.bas[1]*scale.bas[2])
  bb.shp <- bbox(shp)

  draw <- ceiling(N/prop.area) + 1000

  getSample <- function(k = 0){
    if(k == 0){ seedshift <- seed
    }else seedshift <- k*draw + seed
    pts <- RSHalton(n = draw, seeds = seedshift, bases = c(2,3))
    pts[,2] <- pts[,2]*scale.bas[1] + shift.bas[1]
    pts[,3] <- pts[,3]*scale.bas[2] + shift.bas[2]

    #Give points a projection, clip them as needed.
    tmp.order <- (k*draw + 1):((k+1)*draw)
    #Clip the points that are not within our shape's bounding box, keep order That's Important!
    in.shp <- which((pts[,2] >= bb.shp[1,1]) & (pts[,2] <= bb.shp[1,2]) & (pts[,3] >= bb.shp[2,1]) & (pts[,3] <= bb.shp[2,2]))
    pts.coord <- SpatialPointsDataFrame(cbind(pts[in.shp,2],pts[in.shp,3]),proj4string=CRS(nztm), data.frame(SiteID = paste0(island, tmp.order[in.shp])))
    pts.coord <- pts.coord[shp,]
    return(pts.coord)
  }

  pts.sample <- getSample()
  if(nrow(pts.sample) == 0) {
    draw <- draw * 2
    pts.sample <- getSample()
  }

  di <- 1
  while(nrow(pts.sample) < N){
    new.pts <- getSample(k = di)
    if(nrow(new.pts) > 0) pts.sample <- rbind(pts.sample, new.pts)
    di <- di + 1
  }

  pts.sample$SiteOrder <- 1:nrow(pts.sample)

  return(pts.sample[1:N,])
}
