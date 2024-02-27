## Sets INLA mesh for the spde model
require(splancs)
require(rgl)
require(INLA)
require(lattice)
require(concaveman)

fault = fullFault

nsf = length(fault)

x = rep(0, nsf)
y = rep(0, nsf)

for (i in 1:nsf){
  x[i] = fault[[i]]$lon
  y[i] = fault[[i]]$lat
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
# Some parameters added to try and get the mesh to hit all centroids
# hits about 99%...
inla_mesh = inla.mesh.2d(n=2000,
                         loc=xy,
                         boundary=list(concaveInt, hullExt),
                         max.edge=c(12, 1000), # larger max edge for outerior
                         cutoff = 2,
                         min.angle = 10)

spdeMesh = inla_mesh
