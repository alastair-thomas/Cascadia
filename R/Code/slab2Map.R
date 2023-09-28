
# gets the border data for the country and required states / provinces
readCountries <- function(){
  
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
  
  usSF <- st_as_sf(us)
  canadaSF <- st_as_sf(canada)
  
  usStates <- usSF[is.element(usSF$NAME_1,
                                c("California", "Oregon", "Washington",
                                  "Idaho", "Nevada", "Arizona")),]
  canadaProv <- canadaSF[is.element(canadaSF$NAME_1,
                                      c("British Columbia", "Alberta")),]
  
  border <- read_sf('Data/Border/Canada_and_US_Border.shp')
  
  return(list(Canada=canadaProv, US=usStates, canUSBorder=border))
}

# varName ∈ ["depth", "dip", "strike", "thickness", "uncertainty"]
# returns a data frame of Lat, Lon, Variable
#
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


plotGridGG <- function(varName, title){

  # read in the border data
  borders <- readCountries()
  canadaBorders <- ms_simplify(borders$Canada, weighting=0.7) #, tol=0.05, topologyPreserve=TRUE)
  usBorders <- ms_simplify(borders$US, weighting=0.7)
  canUSBorder <- borders$canUSBorder
  
  # read in the slab2 data
  grid <- readGrid(varName)
  
  # create dataframes for labels
  placeNames <- data.frame(Lon = c(-125.5),
                           Lat = c(49.7),
                           Place = c("Vancouver Island"))
  
  stateNames <- data.frame(Lon = c(-121.75, -121.5, -121, -121.5),
                            Lat = c(40.9, 44, 47.5, 50),
                            State = c("California", "Oregon", "Washington", "British Columbia"))
  
  countryNames <- data.frame(Lon=c(-119.5, -119.5),
                              Lat=c(49.2, 48.8),
                              Country=c("Canada", "USA"))
  
  myBreaks <- c(seq(-100,0,100/9))
  
  ggplot() +
    geom_sf(data = usBorders, fill = "white", colour="black") +
    geom_sf(data = canadaBorders, fill = "white", colour="black") +
    geom_sf(data=canUSBorder, linewidth=0.6, colour="black", linetype=11) +
    geom_contour_filled(data=grid[grid$z > -100,], aes(x=lon, y=lat, z = z), alpha=0.6, bins=9) +
    scale_fill_brewer(palette = "OrRd", direction=1) +
    geom_text(data=stateNames, aes(x=Lon, y=Lat, label=State), size=4) +
    geom_text(data=countryNames, aes(x=Lon, y=Lat, label=Country), size=4.5, hjust=1) +
    geom_text(data=placeNames, aes(x=Lon, y=Lat, label=Place), size=3.5, angle=-40) +
    coord_sf(xlim=-c(130, 120), ylim=c(40, 50))+
    labs(x = "Longitude", y = "Latitude", fill=varName) +
    theme_bw() +
    ggtitle(paste("Cascadia Slab2 Data -", title)) +
    theme(panel.background = element_rect('#f0feff')) +
    theme(plot.title = element_text(size=25)) +
    theme(legend.position = "bottom")
}

plotGridGG("depth", "Depth (km)")

ggsave(file="slab2-depth.png", path="Plots/", height=9, width=9)

# project into easting and northing

plotGridPlotly <- function(varName, title){
  
}

