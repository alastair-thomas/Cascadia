meshPlot <- function(){
  borders <- readCountries()
  canadaBorders <- ms_simplify(borders$Canada, weighting=0.7) #, tol=0.05, topologyPreserve=TRUE)
  usBorders <- ms_simplify(borders$US, weighting=0.7)
  canUSBorder <- borders$canUSBorder
  
  
  # create dataframes for labels
  placeNames <- data.frame(Lon = c(-125.5),
                           Lat = c(49.7),
                           Place = c("Vancouver Island"))
  #placeNames = st_as_sf(placeNames, coords = c("Lon", "Lat"), crs=lonLatProj)
  
  stateNames <- data.frame(Lon = c(-121.75, -121.5, -121, -121.5),
                           Lat = c(40.9, 44, 47.5, 50),
                           State = c("California", "Oregon", "Washington", "British Columbia"))
  #stateNames = st_as_sf(stateNames, coords = c("Lon", "Lat"), crs=lonLatProj)
  
  countryNames <- data.frame(Lon=c(-119.5, -119.5),
                             Lat=c(49.2, 48.8),
                             Country=c("Canada", "USA"))
  
  ggplot() +
    geom_sf(data = usBorders, fill = "white", colour="black") +
    geom_sf(data = canadaBorders, fill = "white", colour="black") +
    geom_sf(data=canUSBorder, linewidth=0.6, colour="black", linetype=11) +
    # metR::geom_contour_fill(data=grid2, aes(x=lon, y=lat, z = z), alpha=0.65, bins=100) +
    # scale_fill_distiller(palette = "OrRd", direction=1,
    #                      name = title) +
    #labels = c('0 km', '10 km', '20 km', '30 km', '40 km', '50 km')) +
    #scale_fill_gradient2(name="Depth (Km)") +
    geom_text(data=stateNames, aes(x=Lon, y=Lat, label=State), size=4) +
    geom_text(data=countryNames, aes(x=Lon, y=Lat, label=Country), size=4.5, hjust=1) +
    geom_text(data=placeNames, aes(x=Lon, y=Lat, label=Place), size=3.5, angle=-40) +
    coord_sf(xlim=-c(130, 120), ylim=c(40, 50))+
    theme_bw() +
    theme(panel.background = element_rect('#f0feff')) +
    theme(legend.position = "right",
          legend.key.height = unit(3, 'cm'))
}

faultGeom = discretizeSlab2()


