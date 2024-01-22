mesh = inla_mesh
z = xSum[,1]
proj="northing"

if (proj=="northing"){
  xy = mesh$loc[,1:2]
  latLon = projCSZ(xy, inverse=TRUE, units="km")
  mesh$loc[,1:2] = latLon
  mesh$crs = st_crs("EPSG:4326")
} else{
  xy = projCSZ(mesh$loc[,1:2], units="km")
}


NX = 500
NY = 500
plotData = akima::interp(x=xy[,1], y=xy[,2], z=z, nx=NX, ny=NY)

meshGrid = expand.grid(x = plotData$x, y = plotData$y)
meshGridLatLon = projCSZ(as.matrix(meshGrid), inverse=TRUE, units="km")

plotDF = data.frame(xp=meshGridLatLon[,1],
                    yp=meshGridLatLon[,2],
                    zp=c(plotData$z))

plotDF = na.omit(plotDF)

tileWidth = 3*abs(max(plotDF$xp) - min(plotDF$xp)) / NX
tileHeight = 3*abs(max(plotDF$yp) - min(plotDF$yp)) / NY

g = plotBase(scale=1.5, labels=FALSE, countryBoundary=FALSE)
g = g +
  geom_tile(data=plotDF, aes(x=xp, y=yp, fill=zp),
            width=tileWidth, height=tileHeight) +
  scale_fill_viridis_c(alpha=0.1, option = "viridis", name="legendTitle")

g = g +
  coord_sf(xlim=-c(131, 121), ylim=c(36, 54)) +
  xlab("Longitude") +
  ylab("Latitude")

plot(g)
