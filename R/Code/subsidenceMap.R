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
library(RColorBrewer)
library(viridis)
library(ggrepel)

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

load("Data/Subsidence/DR1.RData")

state_names <- data.frame(Lon = c(-122, -121.5, -121, -121.5),
                          Lat = c(41, 44, 47.5, 50),
                          State = c("California", "Oregon", "Washington", "British Columbia"))

Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

Sites <- data.frame(Site=dr1$Site, Lat=dr1$Lat, Lon=dr1$Lon)
Sites <- aggregate(.~Site, data=Sites, mean)
Sites2 <- data.frame(Site=dr1$Site, Quality=dr1$quality)
Sites2 <- aggregate(.~Site, data=Sites2, Mode)
Sites$Quality <- factor(Sites2$Quality, levels=c(1,2,3), labels=c("Bad","OK", "Good"))

border <- read_sf('Data/Border/Canada_and_US_Border.shp')

country_names <- data.frame(Lon=c(-119.5, -119.5),
                            Lat=c(49.2, 48.8),
                            Country=c("Canada", "USA"))
ggplot() +
  geom_sf(data = us_states, fill = "white", colour="black") +
  geom_sf(data = canada_prov, fill = "white", colour="black") +
  geom_sf(data=border, size=10, colour="black") +
  geom_point(data=Sites, aes(x=Lon, y=Lat, colour=Quality), size=2) +
  scale_colour_brewer(palette="Spectral") +
  geom_text(data=state_names, aes(x=Lon, y=Lat, label=State), size=2.5) +
  geom_text(data=country_names, aes(x=Lon, y=Lat, label=Country), size=3, hjust=1) +
  geom_text_repel(data=Sites, aes(x=Lon, y=Lat, label=Site),
                  size=2.5, hjust=0,
                  xlim = c(-130, -127), force=5,
                  direction="y", box.padding=0.15)+
  coord_sf(xlim=-c(131, 120), ylim=c(40, 50))+
  labs(x = "Longitude", y = "Latitude", colour="Quality:") +
  theme_bw() +
  ggtitle("Cascadia Subsidence Data - Locations and Quality") +
  theme(legend.position="top", legend.margin=margin(t = 0, unit='cm'),
        plot.title = element_text(hjust = 0.5))

ggsave(file="Subsidence.png", path="Plots/", height=6, width=6)
