########################
#Master Sample Code
#Paul van Dam-Bates
########################
#Take in a polygon and generate points from the master sample
#Step 1: Determine Island that the polygon falls in
#Step 2: Use random seed from that island to start Halton Sequence
#Step 3: Ouptut number of points required clipped for that region.
#--------------------------------------------------------------------

library(raster)
library(rgeos)
library(sp)

#' @export
#Halton Sequence:
RSHalton <- function(n = 10, seeds = c(0,0),bases = c(2,3), boxes = 0, J = c(0,0)) {
  ##
  ## Generate n points from a random start d dimensional Halton sequence.
  ##
  ## Inputs:
  ##
  ## n 			sample size
  ## bases    	coprime bases e.g. c(2,3) (Halton Sequence)
  ## seeds  	random seeds  e.g. c(0,0) (Halton Sequence)
  ## boxes 		Index of the Halton Sequence that the Box falls in
  ## B  		Number of boxes to divide Halton Sequence into


  ########### Initialize #########################################
  d <- length(bases);	pts <- mat.or.vec(n, d)
  if (length(seeds) != d){
    seeds <- rep(seeds[1],d)
  }

  boxes <- where2Start(boxes = boxes, J = J, seeds = seeds, bases = bases)
  B <- prod(bases^J)

  ########### Main Loop #########################################
  for (i in 1:d) {
    b <- bases[i];   	u <- seeds[i];
    k <- (rep(u + boxes, ceiling(n/length(boxes))) + rep(seq(0, ceiling(n/length(boxes)) - 1)*B, each = length(boxes)))[1:n]
    xk <- (k %% b)/b;
    for (j in 1:(ceiling(logb(u+n,b)) + 2)) {
      xk <- xk + (floor(k/(b^j)) %% b)/(b^(j+1));
    }
    pts[,i] <- cbind(xk)
  }
  pts <- cbind(k+1-u, pts)
  return(pts)
}

#' @export
#Solve Congruence to get order of boxes:
systemCong <- function(L = c(1/4, 1/3), J = c(2,2), base = c(2,3))
{
  x <- 0:(base[1]^J[1]-1)/1/base[1]^J[1]	#Fix rounding errors!
  y <- 0:(base[2]^J[2]-1)/1/base[2]^J[2]

  L <- c(x[which.min(abs(x - L[1]))], y [which.min(abs(y - L[2]))])

  a1 <- sum( ( floor(L[1]*base[1]^(1:J[1])) %% base[1]) * base[1]^( 1:J[1]-1))
  a2 <- sum( ( floor(L[2]*base[2]^(1:J[2])) %% base[2]) * base[2]^( 1:J[2]-1))
  mod <- base^J
  B <- prod(mod)
  possible <- 0:(B-1)
  sol1 <- possible %% mod[1] == a1
  sol2 <- possible %% mod[2] == a2
  return(possible[which(sol1 + sol2 == 2)])
}

#' @name makeFrame
#' @title Make a Halton grid over the bounding box
#' @export
#Create a Halton Grid over the Bounding Box
makeFrame <- function(base = c(2,3), J = c(2,2), bb)
{
  B <- prod(base^J)
  halt.grid <- raster(extent(as.matrix( bb )), nrow=base[2]^J[2], ncol=base[1]^J[1])
  halt.grid <- rasterToPolygons(halt.grid)
  return(halt.grid)
}

# Wrap a Halton Frame over the sample Shape.
#' @name shape2Frame
#' @title Wrap a Halton Frame over the sample shape
#' @export
shape2Frame <- function(shp, bb = NULL, base = c(2,3), J = c(2,2), projstring = NULL)
{
  if( !is.null( bb))
  {
    scale.bas <- bb[,2] - bb[,1]
    shift.bas <- bb[,1]
  }else{ return("Define Bounding Box Please.")}

  if( is.null( projstring)) {
	projstring <- getProj(island = "South")
	cat("Assuming NZTM Projection\n")
	}
  if(proj4string(shp) != projstring) shp <- spTransform(shp, projstring)

  #Stretch bounding box to Halton Frame Size:
  bb2 <- bbox(shp)
  xy <- (bb2 - shift.bas)/scale.bas
  lx <- floor(xy[1,1] / (1/base[1]^J[1]))/(base[1]^J[1])
  ly <- floor(xy[2,1] / (1/base[2]^J[2]))/(base[2]^J[2])
  ux <- ceiling(xy[1,2] /(1/base[1]^J[1]))/(base[1]^J[1])
  uy <- ceiling(xy[2,2] /(1/base[2]^J[2]))/(base[2]^J[2])
  nx <- (ux-lx)*base[1]^J[1]
  ny <- (uy-ly)*base[2]^J[2]

  bb.new <- data.frame(min = c(lx, ly), max = c(ux, uy), row.names = c("x","y"))
  halt.frame <- raster(extent(as.matrix( bb.new*scale.bas + shift.bas )), nrow=ny, ncol=nx)
  projection(halt.frame) <- projstring
  halt.poly <- rasterToPolygons(halt.frame)
}


#' @export
#Where to start the Halton Sequence
where2Start <- function(J = c(1,1), seeds = c(0,0), bases = c(2,3), boxes = NULL)
{
  B <- prod(bases^J)
  possible <- 0:(B-1)
  sol1 <- possible %% bases[1]^J[1] == (seeds[1] %% bases[1]^J[1])
  sol2 <- possible %% bases[2]^J[2] == (seeds[2] %% bases[2]^J[2])
  boxInit <- possible[which(sol1 + sol2 == 2)]

  if(is.null(boxes)) return()
  boxes <- ifelse(boxes < boxInit, B + (boxes - boxInit), boxes - boxInit)
  return(sort(boxes))
}

#' @name getProj
#' @title Define spatial objects in NZTM projection
#' @export
getProj <- function(island = "South")
{	#NZTM
  if(island != "AucklandIslands")
	return("+proj=tmerc +lat_0=0 +lon_0=173 +k=0.9996 +x_0=1600000 +y_0=10000000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")
	#NZTM for Auckland Islands
	return("+proj=tmerc +lat_0=0 +lon_0=166 +k=1 +x_0=3500000 +y_0=10000000 +ellps=GRS80 +units=m +no_defs")
 }

#' @export
#' @name getBB
#' @title Get the bounding box for other functions
getBB <- function(island = "South")
{
  if(island == "South"){
    bb <- data.frame(min = c(1089354,4747979), max = c(1721164,5516919), row.names = c("x","y"))
  }
  if(island == "North"){
    bb <- data.frame(min = c(1510593,5390569), max = c(2092000,6223164), row.names = c("x","y"))
  }
  if(island == "AucklandIslands"){
    bb <- data.frame(min = c(3490761,4355830), max = c(3525662,4405196), row.names = c("x","y"))
  }
  return(bb)
}


getSeed <- function(island = "South")
{
  if(island == "South"){
    seed <-  c(4887260, 18041662)
  }
  if(island == "North"){
    seed <- c(5137598, 8906854)
  }
  if(island == "AucklandIslands"){
    seed <- c(1925006, 6242772)
  }
  return(seed)
}

#' @name masterSample
#' @title Generate sample points in New Zealand using BAS master sample
#' @description Generates BAS sample points in a specified sample frame based on New Zealand terrestrial mastersample. Users need to specify 'island' the sample frame is in and how many points are required.
#' @export
masterSample <- function(island = "South", shp, N = 100, J = c(0,0)){
  #Define CRS

  if(!island %in% c("South", "North","AucklandIslands")) return("Define the island please.")
  nztm <- getProj(island)
  if(proj4string(shp) != nztm) shp <- spTransform(shp, nztm)
  
  bb <- getBB(island)
  seed <- getSeed(island)

  #Scale and shift Halton to fit into bounding box
  scale.bas <- bb[,2] - bb[,1]
  shift.bas <- bb[,1]

  #We can use Halton Boxes to speed up code when the polygons are small and all over the place.
  #Kind of like magic!
  draw <- N + 5000
  
  if(sum(J) != 0){
    hal.frame <- shape2Frame(shp, J = J, bb = bb, projstring = nztm)
    boxes <- which(rowSums(gIntersects(shp, hal.frame, byid = TRUE)) > 0)
    hal.polys <- hal.frame[boxes,]@polygons
    box.lower <- do.call("rbind", lapply(hal.polys, FUN = function(x){data.frame(x = min(x@Polygons[[1]]@coords[,1]), y = min(x@Polygons[[1]]@coords[,2]))}))
    box.lower <- (t(box.lower) - shift.bas)/scale.bas
    halt.rep <- apply(box.lower, 2, FUN = systemCong, J = J, base = c(2,3))
    B <- prod(c(2,3)^J)
  }else{
    halt.rep = 0; B = 1;
  }

  getSample <- function(k = 0, endPoint = 0){
    if(k == 0){ seedshift <- seed
    }else seedshift <- endPoint + seed
    pts <- RSHalton(n = draw, seeds = seedshift, bases = c(2,3), boxes = halt.rep, J = J)
    pts[,2] <- pts[,2]*scale.bas[1] + shift.bas[1]
    pts[,3] <- pts[,3]*scale.bas[2] + shift.bas[2]
    pts.coord <- SpatialPointsDataFrame(cbind(pts[,2],pts[,3]),proj4string=CRS(nztm), data.frame(SiteID = paste0(island, pts[,1] + endPoint), Count = pts[,1] + endPoint))
    indx <- gIntersects(shp, pts.coord, byid = TRUE)
    pts.coord <- pts.coord[rowSums(indx) > 0,]
    return(pts.coord)
  }

  pts.sample <- getSample()
  if(nrow(pts.sample) == 0) {
    draw <- draw * 2
    pts.sample <- getSample()
  }

  di <- 1
  while(nrow(pts.sample) < N){
    last.pt <- pts.sample@data$Count[nrow(pts.sample)]
    new.pts <- getSample(k = di, endPoint = last.pt)
    if(nrow(new.pts) > 0) pts.sample <- rbind(pts.sample, new.pts)
    di <- di + 1
  }

  return(pts.sample[1:N,])
}

#############
# Take BAS Point and Make Halton Frame around it.
#############
#' @name point2Frame
#' @title Make a halton frame around a BAS point
#' @export
point2Frame <- function(pt, bb = NULL, base = c(2,3), J = c(2,2), projstring = NULL)
{
  if(!is.null(bb))
  {
    scale.bas <- bb[,2] - bb[,1]
    shift.bas <- bb[,1]
  }else{ return("Define Bounding Box Please.")}

  if(is.null(projstring)) projstring <- proj4string(pt)

  xy <- coordinates(pt)
  xy <- (xy - shift.bas)/scale.bas
  lx <- floor(xy[1] / (1/base[1]^J[1]))/(base[1]^J[1])
  ly <- floor(xy[2] / (1/base[2]^J[2]))/(base[2]^J[2])
  frame.order <- systemCong(c(lx, ly), base = base, J = J)
  lx <- lx*scale.bas[1] + shift.bas[1]
  ly <- ly*scale.bas[2] + shift.bas[2]
  framei <- Polygons(list(Polygon(cbind(c(lx,lx,lx + scale.bas[1]/base[1]^J[1], lx + scale.bas[1]/base[1]^J[1],lx), c(ly,ly + scale.bas[2]/base[2]^J[2],ly + scale.bas[2]/base[2]^J[2], ly, ly)))),
                     ID = frame.order)
  framei <- SpatialPolygonsDataFrame(SpatialPolygons(list(framei), proj4string = CRS(projstring)), data = data.frame(Order = frame.order, row.names = frame.order))
  return(framei)
}

#' @name lineSamp
#' @title Generate samples in linear features based on BAS mastersample
#' @export
# Sampling a linear feature:
lineSamp <- function (n = 10, x, seed = 0, halt = TRUE) 
{	
	cc <- do.call("c", coordinates(x))
	# Trick here is to figure out where the "breaks" in the lines occur
	# By index cumulative sum we can identify them
	brks <- sapply(cc, nrow)	
	brk <- cumsum(brks[-length(brks)])	
    cc.df <- do.call("rbind", cc)
    cc.mat <- as.matrix(cc.df)
    lengths = LineLength(cc.mat, longlat = FALSE, sum = FALSE)
	lengths[brk] <- 0	# Remove the length for discontinuities.
    csl = c(0, cumsum(lengths))
    maxl = csl[length(csl)]
    if (halt == TRUE) {
        pts = lineHalton(n, u = seed) * maxl
    }
    else {
        pts = runif(n) * maxl
    }
    int = findInterval(pts, csl, all.inside = TRUE)
    where = (pts - csl[int])/diff(csl)[int]
    xy = cc.mat[int, , drop = FALSE] + where * (cc.mat[int + 
        1, , drop = FALSE] - cc.mat[int, , drop = FALSE])
    samp <- SpatialPoints(xy, proj4string = CRS(proj4string(x)))
	return(samp)
}

#' @export
lineHalton <- function(n = 10, u = 0, b = 5) {
  k <- u:(u+n-1);    xk <- (k %% b)/b;
  for (j in 1:(ceiling(logb(u+n,b)) + 2)) {
    xk <- xk + (floor(k/(b^j)) %% b)/(b^(j+1));
  }
  return(xk)
}
