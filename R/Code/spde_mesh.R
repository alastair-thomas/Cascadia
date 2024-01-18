## Sets INLA mesh for the spde model
require(splancs)
require(rgl)
require(INLA)
require(lattice)
require(concaveman)

setwd("~/Uni/NTNU/Masters Project/CSZ/R/Data/Fault Geometry")
load("simpleFault.RData")
setwd("~/Uni/NTNU/Masters Project/CSZ/R/Code")

nsf = length(simpleFault)

x = rep(0, nsf)
y = rep(0, nsf)

for (i in 1:nsf){
  x[i] = simpleFault[[i]]$lon
  y[i] = simpleFault[[i]]$lat
}

xy = cbind(x, y)

xy = projCSZ(xy, units="km")

concaveHull = concaveman(xy, concavity = 2, length_threshold=0)

concaveInt = inla.mesh.segment(concaveHull, is.bnd=FALSE)

# make sure concaveInt is in a format expected by inla.mesh.2d
concaveInt$idx = rbind(concaveInt$idx, c(nrow(concaveInt$loc), 1))
concaveInt$idx = matrix(as.integer(concaveInt$idx), ncol=2)
concaveInt$grp = matrix(rep(as.integer(1), nrow(concaveInt$loc)), ncol=1)
concaveInt$loc = matrix(concaveInt$loc[,1:2], ncol=2) # I changed this line slightly
concaveInt$loc = concaveInt$loc[nrow(concaveInt$loc):1, ]

hullExt = inla.nonconvex.hull.basic(xy, resolution=150, convex=-.4)

# construct mesh with INLA
inla_mesh = inla.mesh.2d(n=2000,
                         loc=xy,
                         boundary=list(concaveInt, hullExt),
                         max.edge=c(40, 1000), # larger max edge for outerior
                         cutoff = 10)

# INLA SPDE object
# don't need to change
# corresponds to a smoothness of 1

inla_spde = inla.spde2.matern(inla_mesh, alpha=2)


# create the projection matrix
#A = inla.spde.make.A(inla_mesh, loc = xy)

# save
setwd("~/Uni/NTNU/Masters Project/CSZ/R/Data/SPDE")
save(inla_mesh, inla_spde, file="spdeMesh.RData")
setwd("~/Uni/NTNU/Masters Project/CSZ/R/Code")


# A function to plot the mesh onto the base map
#
# mesh - inla.mesh object
# pts - options co-ordinates to plot
# scale - doesn't do anything atm
# proj - which projection the given data is in (should all be in same format)
plotMesh = function(mesh, pts=NULL, scale=1.5, proj="northing"){
  
  if (proj == "northing"){
    xy = mesh$loc[,1:2]
    latLon = projCSZ(xy, inverse=TRUE, units="km")
    mesh$loc[,1:2] = latLon
  }
  
  g = plotBase(scale=scale, labels=FALSE, countryBoundary=FALSE)
  
  g = g +
    gg(mesh,
       edge.color="red", edge.linewidth=0.15,
       int.linewidth=0.75, int.color="blue",
       ext.linewidth = 0.75, ext.color = "blue") +
    coord_sf(xlim=-c(131, 120), ylim=c(36, 54))
  
  # if want to plot the centers of the subfaults
  if (is.null(pts) == FALSE){
    if (proj == "northing"){
      pts = projCSZ(pts, inverse=TRUE, units="km")
    }
    
    print(data.frame(pts))
    
    g = g +
      geom_point(data=data.frame(pts), aes(x=X, y=Y), size=0.1, color="red")
  }
  
  return(g)
}

plotBothMesh = function(fault, mesh, meshProj="northing", scale=1.5){
  if (meshProj == "northing"){
    xy = mesh$loc[,1:2]
    latLon = projCSZ(xy, inverse=TRUE, units="km")
    mesh$loc[,1:2] = latLon
  }
  
  # add the underlying spatial field plot SPDE
  g = plotBase(scale=scale, labels=FALSE)
  
  g = g +
    gg(mesh, edge.color="red",
       edge.linewidth=0.15,
       interior=FALSE,
       exterior=FALSE)
  
  # add the sub faults
  nsf = length(fault)
  
  ids = factor(1:nsf)
  
  depths = rep(0, nsf)
  lons = rep(0, 3*nsf)
  lats = rep(0, 3*nsf)
  
  for (i in 1:nsf){
    lons[(((i-1)*3)+1):((i*3))] = fault[[i]]$corners[,1]
    lats[(((i-1)*3)+1):((i*3))] = fault[[i]]$corners[,2]
  }
  
  values = data.frame(
    id = ids,
    depth = -depths
  )
  
  positions = data.frame(
    id = rep(ids, each = 3),
    x = lons,
    y = lats
  )
  
  # merge together
  datapoly = merge(values, positions, by = c("id"))
  
  # add subfault mesh
  g = g +
      geom_polygon(data=datapoly, aes(x=x, y=y, group=id),
                   color="blue", fill=NA, size=0.15)
  
  
  # Finally add the centroids of subfaults
  x = rep(0, nsf)
  y = rep(0, nsf)
  
  for (i in 1:nsf){
    x[i] = fault[[i]]$lon
    y[i] = fault[[i]]$lat
  }
  
  centers = data.frame(x=x, y=y)
  
  g = g +
    geom_point(data=centers, aes(x=x, y=y),
               colour="green", size=2)
  
  colors = c("Sub Faults" = "blue", "SPDE Mesh" = "red", "Sub Fault Centroids" = "green")
  g = g +
      labs(colour = "Legend") +
      scale_colour_manual(values = colors) +
      coord_sf(xlim=-c(124, 126), ylim=c(41, 43))
  
  return(g)
}
