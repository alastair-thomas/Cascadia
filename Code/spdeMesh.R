getSPDEMesh = function(fault, max.edge=c(15, 1000), cutoff=3){
  require(splancs)
  require(rgl)
  require(INLA)
  require(lattice)
  require(concaveman)
  
  S = length(fault)
  
  x = rep(0, S)
  y = rep(0, S)
  
  for (i in 1:S){
    x[i] = fault[[i]]$lon
    y[i] = fault[[i]]$lat
  }
  
  xy = cbind(x, y)
  xy = projCSZ(xy, units="km")
  
  concaveHull = concaveman(xy)
  
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
                           max.edge=max.edge, # larger max edge for outerior
                           cutoff = cutoff)
  
  A = inla.spde.make.A(inla_mesh,
                       loc=xy)
  
  return(list(mesh=inla_mesh, A=A))
}
