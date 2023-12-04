
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


plotGrid = function(varName = "depth", projection="NULL"){
  # gets the base map in given projection
  g = plotBase(scale=2)
  
  # read in the slab2 data for depth
  grid = readGrid(varName)
  # limit the depth
  if (varName == "depth"){
    grid = grid[grid$z > -30,]
  }
  
  g = g +
    geom_sf(data=grid, aes(color=z)) +
    #geom_sf(data=grid, aes(color=z), stat="contour_filled") +
    scale_fill_distiller(palette = "OrRd", direction=1, name = "Depth (km)") +
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
    g = g + coord_sf(xlim = -c(130, 120), ylim = c(40, 50), crs=st_crs("EPSG:4326"))
  }
  
  return(g)
}

#plotGridGG("depth", "Depth (km)", 50)

#ggsave(file="slab2-depth.png", path="Plots/", height=9, width=9)

