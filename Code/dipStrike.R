varName = "dip"
scale=2
projection="UTM"
legendTitle="Dip Angle (°)"

# gets the base map in given projection
g = plotBase(scale=scale, proj="UTM", labels=FALSE)

# read in the slab2 data for depth
grid = readGrid(varName)

# make strike centered around zero
if (varName == "strike"){
  grid$Z[grid$Z > 100] = grid$Z[grid$Z > 100] - 360
}

# limit to 30km depths
if (varName == "depth"){
  depths = -grid$Z
} else{
  tempGrid = readGrid("depth")
  depths = -tempGrid$Z
}
depthMask = which(depths <= 30)
grid2 = grid[depthMask,]

hull = concaveman(as.matrix(grid2[,1:2]))

# set the map in the correct projection
# works because everything is an sf object
if (projection == "UTM"){
  
  # convert the grid data to easting, northing
  xy = data.frame(x = grid2$X, y = grid2$Y)
  xy = st_as_sf(x=xy, coords = c("x", "y"), crs=st_crs("EPSG:4326"))
  xy = toUTM(xy)
  xy$Z = grid2$Z
  
  # extract as a dataframe
  grid3 = data.frame(X = c(st_coordinates(xy)[,1]),
                     Y = c(st_coordinates(xy)[,2]),
                     Z = grid2$Z)
  # convex hull of grid3
  hull = concaveman(as.matrix(grid3[,1:2]))
  
  # Interpolate onto a regular grid
  xs = seq(min(grid3$X), max(grid3$X), length.out = 200)
  ys = seq(min(grid3$Y), max(grid3$Y), length.out = 200)
  mesh = akima::interp(grid3$X, grid3$Y, grid3$Z, xo=xs, yo=ys)
  grid4 = akima::interp2xyz(mesh, data.frame=T)
  grid4 = na.omit(grid4)
  
  # now take only points inside the convex hull
  inside = point.in.polygon(point.x = grid4$x,
                            point.y = grid4$y,
                            pol.x   = hull[,1],
                            pol.y   = hull[,2])
  # 1 == inside, 0 == outside
  grid4 = grid4[inside == 1,]
  
  # create the limits
  limits = data.frame(x = -c(128.5, 122),
                      y = c(39.8, 50.2))
  limits = st_as_sf(x=limits,
                    coords = c("x", "y"),
                    crs=st_crs("EPSG:4326"))
  limits = toUTM(limits)
  
  # plot
  g = g +
    geom_raster(data=grid4, aes(x=x, y=y, fill=z)) +
    scale_fill_viridis(alpha=0.75, option="mako") +
    labs(fill=legendTitle) +
    theme(legend.key.height = unit(1, "cm")) +
    coord_sf(xlim = c(st_coordinates(limits)[,1]),
             ylim = c(st_coordinates(limits)[,2]),
             crs=st_crs("EPSG:32610"))
}

plot(g)

