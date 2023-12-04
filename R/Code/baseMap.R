
# gets the border data for the country and required states / provinces
readCountries = function(){
  
  setwd("~/Uni/NTNU/Masters Project/CSZ/R")
  # reads in the Canada country data from GADM
  if (file.exists("Data/Canada/gadm/gadm41_CAN_1_pk.rds")){
    canada = readRDS("Data/Canada/gadm/gadm41_CAN_1_pk.rds")
  } else {
    canada = gadm("CAN", level=1, path="Data/Canada", resolution=2) # provinces
  }
  
  # reads in the US country data from GADM
  if (file.exists("Data/US/gadm/gadm41_USA_1_pk.rds")){
    us = readRDS("Data/US/gadm/gadm41_USA_1_pk.rds")
  } else {
    us = gadm("USA", level=1, path="Data/US", resolution=2) # states
  }
  
  usSF = st_as_sf(us, crs=st_crs("EPSG:4326"))
  canadaSF = st_as_sf(canada, crs=st_crs("EPSG:4326"))
  
  usStates = usSF[is.element(usSF$NAME_1,
                              c("California", "Oregon", "Washington",
                                "Idaho", "Nevada", "Arizona")),]
  canadaProv = canadaSF[is.element(canadaSF$NAME_1,
                                    c("British Columbia", "Alberta")),]
  
  border = read_sf('Data/Border/Canada_and_US_Border.shp')
  border2 = border[is.element(border$SectionEng, c('Straits of Georgia and Juan de Fuca',
                                                    'The 49th Parallel Boundary, Columbia Valley to Pacific Ocean',
                                                    'The 49th Parallel Boundary, Similkameen River to Columbia Valley',
                                                    'The 49th Parallel Boundary, West Kootenay to the Similkameen River')), ]
  
  return(list(Canada=canadaProv, US=usStates, canUSBorder=border2))
}

# project into UTM co-ordinate system
# requires an sf object as x
toUTM <- function(x){
  
  # determine the "to" projection systems
  utmProj = st_crs("EPSG:32610")
  
  # transform coordinates
  out = st_transform(x, utmProj)
  
  return(out)
}

# project to the Lon/Lat co-ordinate system
# requires an sf object
fromUTM <- function(x){
  
  # the correct projection system
  lonLatProj = st_crs("EPSG:4326")
  
  # transform coordinates
  out = st_transform(x, lonLatProj)
  
  return(out)
  
}

plotBase = function(scale=1, labels=TRUE){
  
  # read in the border data
  borders = readCountries()
  
  # extract each country and simplfy the geometry so it plots nicely
  canadaBorders = ms_simplify(borders$Canada, weighting=0.7, keep_shapes=TRUE)
  usBorders = ms_simplify(borders$US, weighting=0.7, keep_shapes=TRUE)
  canUSBorder = ms_simplify(borders$canUSBorder, weighting=0.7, keep_shapes=TRUE)
  
  placeNames = data.frame(Lon = c(-126),
                           Lat = c(49.9),
                           Place = c("Vancouver Island"))
  placeNames = st_as_sf(x=placeNames, coords = c("Lon", "Lat"), crs=st_crs("EPSG:4326"))
  
  stateNames = data.frame(Lon = c(-121.75, -121.5, -121, -121.5),
                           Lat = c(40.9, 44, 47.5, 50),
                           State = c("California", "Oregon", "Washington", "British Columbia"))
  stateNames = st_as_sf(x=stateNames, coords = c("Lon", "Lat"), crs=st_crs("EPSG:4326"))
  
  countryNames = data.frame(Lon=c(-119.5, -119.5),
                             Lat=c(49.2, 48.8),
                             Country=c("Canada", "USA"))
  countryNames = st_as_sf(x=countryNames, coords = c("Lon", "Lat"), crs=st_crs("EPSG:4326"))
  
  # add land
  g = ggplot() +
      geom_sf(data = usBorders, fill = "white", colour="black") +
      geom_sf(data = canadaBorders, fill = "white", colour="black") +
      geom_sf(data=canUSBorder, linewidth=scale*0.6, colour="black", linetype=11)
  
  if (labels){
    g = g +
        geom_sf_text(data=stateNames, aes(label=State), size=scale*2.5) +
        geom_sf_text(data=countryNames, aes(label=Country), size=scale*3, hjust=1) +
        geom_sf_text(data=placeNames, aes(label=Place), size=scale*2, angle=-40)
  }
      
    # control the appearance
  g = g +
      theme_bw() +
      theme(panel.background = element_rect('#f0feff')) +
      labs(x = "",
           y = "",
           title = "")
  
  return(g)
}
