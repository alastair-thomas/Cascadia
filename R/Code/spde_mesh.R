## Sets INLA mesh for the spde model
require(splancs)
require(rgl)
require(INLA)
require(lattice)
require(concaveman)

# load in the fault geometry
setwd("~/Uni/NTNU/Masters Project/CSZ/R/Data/Fault Geometry")
load("fault.RData") # loads in triangular fault geometry as triGeomFault
setwd("~/Uni/NTNU/Masters Project/CSZ/R/Code")

K = length(triGeomFull) # number of subfaults
x = rep(0, K)
y = rep(0, K)

# extract centers of subfaults
for (i in 1:K){
  x[i] = triGeomFull[[i]]$lon
  y[i] = triGeomFull[[i]]$lat
}

# This should use the midpoint of our subfaults
locations = cbind(x, y)

concaveHull = concaveman(locations)

concaveInt = inla.mesh.segment(concaveHull, is.bnd=FALSE)
# make sure concaveInt is in a format expected by inla.mesh.2d
concaveInt$idx = rbind(concaveInt$idx, c(nrow(concaveInt$loc), 1))
concaveInt$idx = matrix(as.integer(concaveInt$idx), ncol=2)
concaveInt$grp = matrix(rep(as.integer(1), nrow(concaveInt$loc)), ncol=1)
concaveInt$loc = matrix(concaveInt$loc, ncol=2)
concaveInt$loc = concaveInt$loc[nrow(concaveInt$loc):1, ]

hullExt = inla.nonconvex.hull.basic(locations, resolution=150, convex=-.4)

# construct mesh with INLA
inla_mesh = inla.mesh.2d(n=2000,
                         loc=locations,
                         boundary=list(concaveInt, hullExt),
                         cutoff = 0.01)

# plot inla mesh object to check it looks okay.
plot(inla_mesh)

# INLA SPDE object
# don't need to change
# corresponds to a smoothness of 1
inla_spde = inla.spde2.matern(inla_mesh, alpha=2)

# create the projection matrix
A = inla.spde.make.A(inla_mesh, loc = locations)

# save
setwd("~/Uni/NTNU/Masters Project/CSZ/R/Data/SPDE")
save(inla_mesh, inla_spde, A, file="spde_mesh.RData")
setwd("~/Uni/NTNU/Masters Project/CSZ/R/Code")

plotMesh = function(mesh, scale=1.5){
  g = plotBase(scale=scale, labels=FALSE)
  
  g = g +
    gg(mesh, edge.color="#D35400", ext.linewidth = 0.75, ext.color = "#4d4c4b") +
    coord_sf(xlim=-c(131, 120), ylim=c(36.8, 52.5))
  
  return(g)
}
