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
library(GADMTools)
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
  x = x - 360
  z = as.vector(z)
  df <- data.frame(lon = x, lat = y, z = z)
  df = na.omit(df)
  rownames(df) <- 1:nrow(df) 
  return(df)
}

dip_df <- readGrid("dip")
depth_df <- readGrid("depth")
strike_df <- readGrid("strike")
thick_df <- readGrid("thickness")
unc_df <- readGrid("uncertainty")

plotGrid <- function(df, varName){
  ggplot() +
    geom_contour_filled(data=df, aes(x=lon, y=lat, z = z)) +
    geom_sf(data = us_states, fill = NA, colour="black") +
    geom_sf(data = canada_prov, fill = NA, colour="black") +
    coord_sf(xlim=-c(130, 117), ylim=c(35, 55))+
    labs(x = "Longitude", y = "Latitude", fill=varName) +
    theme_bw() +
    ggtitle(paste("Slab2", varName, "- Cascadia"))
}

plotGrid(depth_df, "Depth")

load("Data/Subsidence/DR1.RData")

ggplot() +
  geom_sf(data = us_states, fill = "white", colour="black") +
  geom_sf(data = canada_prov, fill = "white", colour="black") +
  geom_point(data=dr1, aes(x=Lon, y=Lat, colour=quality), size=1) +
  coord_sf(xlim=-c(130, 120), ylim=c(40, 50))+
  labs(x = "Longitude", y = "Latitude", colour="Quality") +
  theme_bw() +
  ggtitle(paste("Subsidence Samples", "- Cascadia"))
