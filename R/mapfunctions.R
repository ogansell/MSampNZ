library(rgdal)
library(mapview)

#Function for convert crs and then outputting a mapview map for easy viewing of spatial data generated in NZTM
#' @export
toMap <- function(x){
  nztm <-"+proj=tmerc +lat_0=0 +lon_0=173 +k=0.9996 +x_0=1600000 +y_0=10000000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"
  prj.LatLong <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")#Set CRS to convert to WGS84
  proj4string(x) = CRS(nztm)
  x <-spTransform(x, prj.LatLong)#Convert CRS to WGS84
  mapview::mapview(x)
}

#Function to make it easy to assign NZTM CRS to spatial data and then convert to WGS for use in leaflet
#' @export
coOrd <- function(x){
  nztm <-"+proj=tmerc +lat_0=0 +lon_0=173 +k=0.9996 +x_0=1600000 +y_0=10000000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"
  prj.LatLong <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")#Set CRS to convert to WGS84
  proj4string(x) = CRS(nztm)
  x <-spTransform(x, prj.LatLong)#Convert CRS to WGS84
}

#Quick conversion to NZTM
#' @export
proj.nztm <- function(x){
  nztm <-"+proj=tmerc +lat_0=0 +lon_0=173 +k=0.9996 +x_0=1600000 +y_0=10000000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"
  proj4string(x) = CRS(nztm)
}

#Function to make it easy make leaflet map to view shapefile and mastersample points
#' @export
viewMap <- function(points,rad,shp){
  leaflet(options = leafletOptions(zoomControl = FALSE))%>%
    addTiles("http://tiles-a.data-cdn.linz.govt.nz/services;key=d27d21709f324848b2d1ffc5e2220036/tiles/v4/layer=767/EPSG:3857/{z}/{x}/{y}.png",
             group = "Topo") %>% #Add plain topo 50 imagery
    addProviderTiles("Esri.WorldImagery", group = "Satellite")%>%
    addPolygons(data = coOrd(shp), popup = popupTable(shp), group = "shp")%>%
    addCircles(data = coOrd(points), popup = popupTable(points),group = "points", color = "red",
               radius = rad)%>%
    addLayersControl(baseGroups = c("Topo","Satellite"),overlayGroups = c("shp", "points"),
                     options = layersControlOptions(collapsed = FALSE))%>%
    addMeasure(primaryLengthUnit = "meters", primaryAreaUnit = "hectares", secondaryLengthUnit = NULL )
}

#Function to make it easy make leaflet map to view shapefile and mastersample points and transects generated from mastersample
#' @export
viewLinesandpoints <-function (points, rad, shp, polyline) {
  leaflet(options = leafletOptions(zoomControl = FALSE)) %>%
    addTiles("http://tiles-a.data-cdn.linz.govt.nz/services;key=d27d21709f324848b2d1ffc5e2220036/tiles/v4/layer=767/EPSG:3857/{z}/{x}/{y}.png",
             group = "Topo") %>%
    addProviderTiles("Esri.WorldImagery", group = "Satellite") %>%
    addPolygons(data = coord(shp), group = "shp") %>%
    addPolylines(data = coord(polyline), popup = popupTable(polyline), color = "green",group = "polyline")%>%
    addCircles(data = coord(points), popup = popupTable(points), group = "points", color = "red", radius = rad) %>%
    addLayersControl(baseGroups = c("Topo","Satellite"), overlayGroups = c("shp", "points","polyline"), options = layersControlOptions(collapsed = FALSE)) %>%
    addMeasure(primaryLengthUnit = "meters", primaryAreaUnit = "hectares",
               secondaryLengthUnit = NULL)
}
