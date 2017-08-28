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
library(rgdal)
library(raster)

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
master.sample <- function(island = "South", shp, N = 100){
	#Master Sample seed for South Island, chosen as first random start that fell into SI
	#seed.si <-  c(4887260, 18041662)
	#seed.ni <- c(5137598, 8906854)

	#Start with by defining bounding box in NZTM
	bb.si <- data.frame(min = c(1089354,4747979), max = c(1721164,5516919), row.names = c("x","y"))
	bb.ni <- data.frame(min = c(1510593,5390569), max = c(2092000,6223164), row.names = c("x","y"))
	
	#Define CRS
	nztm <-"+proj=tmerc +lat_0=0 +lon_0=173 +k=0.9996 +x_0=1600000 +y_0=10000000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"
	
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
	d.b <- bb[,2] - bb[,1]
	scale.bas <- max(d.b)
	shift.bas <- bb[,1]

	pts <- RSHalton(n = 10000000, seeds = seed, bases = c(2,3))
	pts[,2] <- pts[,2]*scale.bas + shift.bas[1]
	pts[,3] <- pts[,3]*scale.bas + shift.bas[2]

	draw <- 1000000
	
	getSample <- function(k = 0){
		seedshift <- ifelse(k == 0, 0, k*draw + seed + 1)
		pts <- RSHalton(n = draw, seeds = seed + seedshift, bases = c(2,3))
		pts[,2] <- pts[,2]*scale.bas + shift.bas[1]
		pts[,3] <- pts[,3]*scale.bas + shift.bas[2]
		
		#Give points a projection, clip them as needed.
		tmp.order <- (k*draw + 1):((k+1)*draw)
		pts.coord <- SpatialPointsDataFrame(cbind(pts[,2],pts[,3]),proj4string=CRS("+proj=tmerc +lat_0=0 +lon_0=173 +k=0.9996 +x_0=1600000 +y_0=10000000 +ellps=GRS80 +units=m +no_defs"), data.frame(SiteOrder = tmp.order))
		return(pts.coord)
	}
	
	pts.sample <- getSample()
	pts.sample <- pts.sample[shp, ]

	
	if(nrow(pts.sample) < N){
		di <- 1
		while(nrow(pts.sample) < N){
			new.pts <- getSample(k = di)
			new.pts <- new.pts[shp, ]
			pts.sample <- rbind(pts.sample, new.pts)
			di <- di + 1
		}
		return(pts.sample[1:N,])
	} else{
		return(pts.sample[1:N,])
	}
}
