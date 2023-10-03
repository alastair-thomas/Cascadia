
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
  border2 <- border[is.element(border$SectionEng, c('Straits of Georgia and Juan de Fuca',
                                                    'The 49th Parallel Boundary, Columbia Valley to Pacific Ocean',
                                                    'The 49th Parallel Boundary, Similkameen River to Columbia Valley',
                                                    'The 49th Parallel Boundary, West Kootenay to the Similkameen River')), ]
  
  return(list(Canada=canadaProv, US=usStates, canUSBorder=border2))
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
  
  if (varName == "depth"){
    z = -1 * z
  }
  
  df <- data.frame(lon = allLon, lat = allLat, z = z)
  df = na.omit(df)
  rownames(df) <- 1:nrow(df) 
  return(df)
}

# project into UTM co-ordinate system
toUTM <- function(x, isDataFrame=FALSE){
  
  # determine the "to" projection systems
  utmProj = st_crs("EPSG:32610")
  #lonLatProj = st_crs("EPSG:4326")
  
  toProj = utmProj
  #fromProj = lonLatProj
  
  # convert to an sf object if input is dataframe
  if (isDataFrame){
    xConvert = st_as_sf(x, coords = c("Lon", "Lat"))
  }
  else{
    xConvert = x
  }
  
  # transform coordinates and convert back to a matrix
  out = st_transform(xConvert, toProj)
  
  if (isDataFrame){
    out = st_coordinates(out)[,1:2]
  }
  
  return(out)
}

fromUTM <- function(x){
  
}

testPlot <- function(){
  # read in the border data
  borders <- readCountries()
  canadaBorders = ms_simplify(borders$Canada, weighting=0.7) #, tol=0.05, topologyPreserve=TRUE)
  usBorders = ms_simplify(borders$US, weighting=0.7)
  canUSBorder = borders$canUSBorder
  
  canadaBorders = toUTM(canadaBorders)
  usBorders = toUTM(usBorders)
  canUSBorder = toUTM(canUSBorder)
  
  #browser()
  
  utmProj = st_crs("EPSG:32610")
  
  ggplot() +
    geom_sf(data = usBorders, fill = "white", colour="black") +
    geom_sf(data = canadaBorders, fill = "white", colour="black") +
    geom_sf(data=canUSBorder, linewidth=0.6, colour="black", linetype=11)
    #coord_sf(xlim=-c(130, 120), ylim=c(40, 50), datum = utmProj)
}

testPlot()

plotGridGG <- function(varName, title, limitZ = NA){
  lonLatProj = st_crs("EPSG:4326")
  
  # read in the border data
  borders <- readCountries()
  canadaBorders <- ms_simplify(borders$Canada, weighting=0.7) #, tol=0.05, topologyPreserve=TRUE)
  usBorders <- ms_simplify(borders$US, weighting=0.7)
  canUSBorder <- borders$canUSBorder
  
  #canadaBorders = toUTM(canadaBorders)
  #usBorders = toUTM(usBorders)
  # = toUTM(canUSBorder)
  
  # read in the slab2 data
  grid <- readGrid(varName)
  if (is.na(limitZ) == FALSE){
    grid2 <- grid[grid$z < limitZ,]
  }
  
  # create dataframes for labels
  placeNames <- data.frame(Lon = c(-125.5),
                           Lat = c(49.7),
                           Place = c("Vancouver Island"))
  placeNames = st_as_sf(placeNames, coords = c("Lon", "Lat"), crs=lonLatProj)
  
  stateNames <- data.frame(Lon = c(-121.75, -121.5, -121, -121.5),
                            Lat = c(40.9, 44, 47.5, 50),
                            State = c("California", "Oregon", "Washington", "British Columbia"))
  stateNames = st_as_sf(stateNames, coords = c("Lon", "Lat"), crs=lonLatProj)
  
  countryNames <- data.frame(Lon=c(-119.5, -119.5),
                              Lat=c(49.2, 48.8),
                              Country=c("Canada", "USA"))
  countryNames = st_as_sf(countryNames, coords = c("Lon", "Lat"), crs=lonLatProj)
  
  #browser()
  utmProj = st_crs("EPSG:32610")
  
  ggplot() +
    geom_sf(data = usBorders, fill = "white", colour="black") +
    geom_sf(data = canadaBorders, fill = "white", colour="black") +
    geom_sf(data=canUSBorder, linewidth=0.6, colour="black", linetype=11) +
    metR::geom_contour_fill(data=grid2, aes(x=lon, y=lat, z = z), alpha=0.65, bins=100) +
    scale_fill_distiller(palette = "OrRd", direction=1,
                         name = title) +
                         #labels = c('0 km', '10 km', '20 km', '30 km', '40 km', '50 km')) +
    #scale_fill_gradient2(name="Depth (Km)") +
    geom_sf_text(data=stateNames, aes(x=Lon, y=Lat, label=State), size=4) +
    geom_sf_text(data=countryNames, aes(x=Lon, y=Lat, label=Country), size=4.5, hjust=1) +
    geom_sf_text(data=placeNames, aes(x=Lon, y=Lat, label=Place), size=3.5, angle=-40) +
    coord_sf(xlim=-c(130, 120), ylim=c(40, 50), datum = utmProj)+
    labs(x = "Easting", y = "Northing", fill=varName) +
    theme_bw() +
    theme(panel.background = element_rect('#f0feff')) +
    theme(legend.position = "right",
          legend.key.height = unit(3, 'cm'))
}

plotGridGG("depth", "Depth (km)", 50)

ggsave(file="slab2-depth.png", path="Plots/", height=9, width=9)

