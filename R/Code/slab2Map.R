
# varName ∈ ["depth", "dip", "strike", "thickness", "uncertainty"]
# returns a data frame of Lat, Lon, Variable
#
readGrid = function(varName){
  Grid = nc_open(paste("Data/Slab2/", varName, ".grd", sep=""))
  x = ncvar_get(Grid,"x")
  y = ncvar_get(Grid,"y")
  z = ncvar_get(Grid,"z")
  
  allLat = rep(y, each=length(x))
  allLon = rep(x, length(y))
  allLon = allLon - 360
  
  z = as.vector(z)
  
  grid = data.frame(Lon = allLon, Lat = allLat, z = z)
  grid = na.omit(grid)
  rownames(grid) = 1:nrow(grid)
  grid = st_as_sf(x=grid, coords = c("Lon", "Lat"), crs=st_crs("EPSG:4326"))
  return(grid)
}


plotGrid = function(scale=2, varName = "depth", projection="NULL"){
  # gets the base map in given projection
  g = plotBase(scale=scale, labels=FALSE)
  
  # read in the slab2 data for depth
  grid = readGrid(varName)
  
  # Extract point data
  point_data = st_coordinates(grid)
  
  # Convert to a data frame
  grid2 = as.data.frame(point_data)
  grid2$Z = -grid$z
  
  # limit the depth
  
  #grid$Discrete <- cut(grid$z, seq(0,max(grid$depth),1), include.lowest=T)
  
  g = g +
    geom_tile(data=grid2, aes(x=X, y=Y, fill=Z), alpha=0.75) +
    scale_fill_gradient(low = "red", high = "blue", name = "Depth (km)", trans="reverse") +
    theme(legend.position = "right", legend.key.height = unit(3, 'cm'))

  # set the map in the correct projection
  # works because everything is an sf object
  if (projection == "UTM"){
    
    limits = data.frame(x = -c(130, 120), y = c(40, 50))
    limits = st_as_sf(x=limits, coords = c("x", "y"), crs=st_crs("EPSG:4326"))
    limits = toUTM(limits)
    
    g = g + coord_sf(xlim = c(st_coordinates(limits)[,1]), ylim = c(st_coordinates(limits)[,2]), crs=st_crs("EPSG:32610"))
  }
  else{
    g = g + coord_sf(xlim = -c(129, 119), ylim = c(39, 50), crs=st_crs("EPSG:4326"))
  }
  
  return(g)
}

#plotGrid()
#ggsave(file="Slab2.png", path="Plots/", height=9, width=6)

