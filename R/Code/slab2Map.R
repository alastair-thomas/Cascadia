library(maps)
library(mapproj)
library(mapdata)
library(rgeos)
library(maptools)
library(sp)
library(raster)
library(rgdal)
library(geodata)
library(dplyr)
library(ggplot2)
library(sf)
library(plotly)
library(ncdf4)
library(stringr)

# reads in the Canada country data from GADM
if (file.exists("Data/Canada/gadm/gadm41_CAN_1_pk.rds")){
  canada <- readRDS("Data/Canada/gadm/gadm41_CAN_1_pk.rds")
} else {
  canada <- gadm("CAN", level=1, path="Data/Canada", resolution=2) # provinces
}

# reads in the US country data from GADM
if (file.exists("Data/US/gadm/gadm41_USA_1_pk.rds")){
  us <- readRDS("Data/US/gadm/gadm41_USA_1_pk.rds")
} else {
  us <- gadm("USA", level=1, path="Data/US", resolution=2) # states
}

us_sf <- st_as_sf(us)
canada_sf <- st_as_sf(canada)

us_states <- us_sf[is.element(us_sf$NAME_1,
                              c("California", "Oregon", "Washington",
                                "Idaho", "Nevada", "Arizona")),]
canada_prov <- canada_sf[is.element(canada_sf$NAME_1,
                                    c("British Columbia", "Alberta")),]

state_names <- data.frame(lon = c(-123, 121),
                          lat = c(40, 43),
                          text = c("California", "Oregon"))

readGrid <- function(varName){
  Grid <- nc_open(paste("Data/Slab2/", varName, ".grd", sep=""))
  x <- ncvar_get(Grid,"x")
  y <- ncvar_get(Grid,"y")
  z <- ncvar_get(Grid,"z")
  
  allLat <- rep(y, each=length(x))
  allLon <- rep(x, length(y))
  allLon <- allLon - 360
  
  z = as.vector(z)
  
  df <- data.frame(lon = allLon, lat = allLat, z = z)
  df = na.omit(df)
  rownames(df) <- 1:nrow(df) 
  return(df)
}

dip_df <- readGrid("dip")
depth_df <- readGrid("depth")
strike_df <- readGrid("strike")
thick_df <- readGrid("thickness")
unc_df <- readGrid("uncertainty")

depth_df$z <- abs(depth_df$z)

plotGrid <- function(df, varName){
  
  state_names <- data.frame(Lon = c(-121, -121, -121, -121.5),
                            Lat = c(41, 44, 47.5, 50),
                            State = c("California", "Oregon", "Washington", "British Columbia"))
  
  country_names <- data.frame(Lon=c(-119.5, -119.5),
                              Lat=c(49.2, 48.8),
                              Country=c("Canada", "USA"))
  
  border <- read_sf('Data/Border/Canada_and_US_Border.shp')
  
  
  
  ggplot() +
    geom_sf(data = us_states, fill = "white", colour="black") +
    geom_sf(data = canada_prov, fill = "white", colour="black") +
    geom_sf(data=border, size=10, colour="black") +
    geom_contour_filled(data=df[depth_df$z < 50,], aes(x=lon, y=lat, z = z), alpha=0.6) +
    scale_fill_brewer(palette = "Spectral") +
    geom_text(data=state_names, aes(x=Lon, y=Lat, label=State), size=2.5) +
    geom_text(data=country_names, aes(x=Lon, y=Lat, label=Country), size=3, hjust=1) +
    coord_sf(xlim=-c(130, 120), ylim=c(40, 50))+
    labs(x = "Longitude", y = "Latitude", fill=varName) +
    theme_bw() +
    ggtitle(paste("Cascadia Slab2 Data -", varName))
}


plotGrid(unc_df, "Uncertainty (km)")

ggsave(file="slab2-uncertainty.png", path="Plots/", height=6, width=6)

# project into easting and northing



