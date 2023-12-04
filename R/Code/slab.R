# script for loading and discretizing Slab 2.0 data

# Construct the full triangulated geometry of the fault for use in the Okada 
# model.
# n: number of nodes to form triangulation from (each node is the corner of 
#     several triangles)
# max.n: The maximum number of vertices allowed, overriding max.edge only 
#        (default=-1, meaning no limit). One or two values, where the second 
#        value gives the number of additional vertices allowed for the 
#        extension.
# max.edge: maximum triangle edge length in km (for making triangle regular)
# maxDepth: maximum depth of the fault geometry
# ...: other inputs passed to inla.mesh.2d
getFullFaultGeom = function(triangulatedGeom=NULL, n=2000, max.n=-1, max.edge=c(15, 100), 
                            maxDepth=30, ...) {
  
  # construct the triangulated geometry if need be
  if(is.null(triangulatedGeom)) {
    triangulatedGeom = discretizeSlab2(n=n, max.n=max.n, max.edge=max.edge, 
                                       maxDepth=maxDepth, ...)
    
    # convert coordinates back to lon/lat. Only corners is used
    for(i in 1:length(triangulatedGeom$corners)) {
      triangulatedGeom$corners[[i]][,1:2] = projCSZ(as.matrix(triangulatedGeom$corners[[i]][,1:2]), inverse=TRUE, units="km")
    }
  }
  cornerList = triangulatedGeom$corners
  
  # now convert to the format used by okadaTri
  faultGeom = list()
  # browser()
  for(i in 1:length(cornerList)) {
    thisCorners = as.matrix(cornerList[[i]])
    thisCorners[,3] = thisCorners[,3]*10^3 # km to m
    
    # calculate the geometry of the subfault
    triGeom = calculate_geometry_triangles(cbind(projCSZ(thisCorners[,1:2], units="m"), thisCorners[,3]))
    
    # reset to unprojected coordinates as expect by okadaTri
    triGeom$corners = thisCorners
    lonLat = projCSZ(cbind(triGeom$lon, triGeom$lat), inverse=TRUE, units="m")
    triGeom$lon = lonLat[1]
    triGeom$lat = lonLat[2]
    
    # append to list of subfaults
    faultGeom = c(faultGeom, list(triGeom))
  }
  
  # return results
  faultGeom
}

# loads the Slab 2.0 model output points and depths
loadSlab2 = function() {
  # need to change to where slab2 data is stored
  slabDir = "C://users/alast/OneDrive/Documents/Uni/NTNU/Masters Project/CSZ/R/Data/Slab2/"
  
  setwd(slabDir)
  # doesn't work/bad depth values due to NAs:
  # slab = read.csv(paste0(slabDir, "cas_slab2_nod_09.04.23.csv"))
  # slab
  
  require(ncdf4)
  ncdat <- nc_open("depth.grd")
  depth = ncvar_get(ncdat, "z")
  lon = ncvar_get(ncdat, "x")
  lat = ncvar_get(ncdat, "y")
  
  allLats = rep(lat, each=length(lon))
  allLons = rep(lon, length(lat))
  allDepths = c(depth)
  
  ncdat <- nc_open("dip.grd")
  dip = ncvar_get(ncdat, "z")
  allDips = c(dip)
  
  ncdat <- nc_open("strike.grd")
  strike = ncvar_get(ncdat, "z")
  allStrikes = c(strike)
  
  out = data.frame(list(lon=allLons, lat=allLats, depth=allDepths, dip=allDips, strike=allStrikes))
  
  extentI = is.finite(allDepths)
  
  setwd("C://users/alast/OneDrive/Documents/Uni/NTNU/Masters Project/CSZ/R/Code")
  
  out[extentI,]
}

# Discretizes the Slab 2.0 geometry into a triangulation mesh.
# NOTE: this function does not currently vary the resolution of the mesh as a 
#     function of the depth, since this will likely not substantially impact the 
#     subsidence along the coast.
# n: number of nodes to form triangulation from (each node is the corner of 
#     several triangles)
# max.n: The maximum number of vertices allowed, overriding max.edge only 
#        (default=-1, meaning no limit). One or two values, where the second 
#        value gives the number of additional vertices allowed for the 
#        extension.
# max.edge: maximum triangle edge length in km (for making triangle regular)
# maxDepth: maximum depth of the fault geometry
# cutoff: minimum allowable edge length in km
# method: either linear, nearest neighbor, r Gaussian kernel smoother 
#         interpolation of the fault geometry. Defaults to linear
# ...: other inputs passed to inla.mesh.2d
discretizeSlab2 = function(n=2000, max.n=-1, max.edge=c(15, 100), maxDepth=30, 
                           cutoff=3, method=c("linear", "NN", "kernel"), ...) {
  method = match.arg(method)
  
  # first load in the Slab 2.0 geometry
  slab = loadSlab2()
  lonLat = cbind(slab$lon, slab$lat)
  depths = slab$depth
  
  # get projected coordinates
  xy = projCSZ(lonLat)
  
  # if(FALSE) {
  #   plotWithColor(xy[,1], xy[,2], depths, pch=19, cex=.3, xlab="Easting (km)", 
  #                 ylab="Northing (km)")
  # }
  
  # keep only points with appropriate depth
  goodI = abs(depths) <= maxDepth
  lonLat = lonLat[goodI,]
  xy = xy[goodI,]
  depths = depths[goodI]
  
  # compute concave hull of slab geometry
  # hullI = chull(xy)
  # xyHull = xy[hullI,]
  
  require(concaveman)
  
  concaveHull = concaveman(xy)
  
  #plot(xy, pch=".", asp=1)
  #polygon(concaveHull, border="blue")
  
  concaveInt = inla.mesh.segment(concaveHull, is.bnd=FALSE)
  # make sure concaveInt is in a format expected by inla.mesh.2d
  concaveInt$idx = rbind(concaveInt$idx, c(nrow(concaveInt$loc), 1))
  concaveInt$idx = matrix(as.integer(concaveInt$idx), ncol=2)
  concaveInt$grp = matrix(rep(as.integer(1), nrow(concaveInt$loc)), ncol=1)
  concaveInt$loc = matrix(concaveInt$loc, ncol=2)
  concaveInt$loc = concaveInt$loc[nrow(concaveInt$loc):1, ]
  
  hullExt = inla.nonconvex.hull.basic(xy, resolution=150, convex=-.4)
  
  # construct mesh with INLA
  mesh = inla.mesh.2d(n=n, boundary=list(concaveInt, hullExt), max.n=max.n, 
                      max.edge=max.edge, cutoff=cutoff, ...)
  
  # if(FALSE) {
  #   plotWithColor(xy[,1], xy[,2], depths, pch=19, cex=.3, xlab="Easting (km)", 
  #                 ylab="Northing (km)")
  #   polygon(xyHull[,1], xyHull[,2], border="green")
  # }
  # 
  # if(FALSE) {
  #   # old way to construct interior boundary:
  #   # hullInt = inla.nonconvex.hull.basic(xy, resolution=350, convex=-.012)
  #   # xyHull = hullInt$loc
  #   # 
  #   # mesh = inla.mesh.2d(n=n, boundary=list(hullInt, hullExt), max.n=max.n, 
  #   #                     max.edge=max.edge, cutoff=cutoff, ...)
  #   
  #   # plot the mesh
  #   plot(mesh, asp=1)
  #   points(xy[,1], xy[,2], col="red", pch=".")
  #   plot(xy[,1], xy[,2], pch=".", col="blue")
  #   polygon(concaveHull)
  # }
  
  faultGeom = getGeomFromMesh(mesh, extent=concaveHull, maxDepth=maxDepth, 
                              method=method)
  
  
  return(list(mesh=mesh, geom=faultGeom, extent=concaveHull, maxDepth=maxDepth))
}

# given a triangulated mesh constructed from discretizeSlab2, constructs 
# a list of triangle centers and corners (also lists of external triangles 
# that are later ignored). Calculated the depths of each point in the geometry.
# inputs:
# mesh: triangulated mesh object from inla.mesh.2d
# extent: points of a polygon describing the extent of the fault geometry
# maxDepth: maximum depth of the fault geometry to use for depth interpolation. 
#           Should be slightly deeper than the extent, e.g. 1km deeper
# method: method used in getPointDepths for depth interpolation
getGeomFromMesh = function(mesh, extent, maxDepth=30, method=c("linear", "NN", "kernel")) {
  method = match.arg(method)
  
  corners = mesh$loc[,1:2] # corners of the triangles, i.e. vertices
  
  # t: triangles, v: vertices
  tv = mesh$graph$tv
  vt = mesh$graph$vt
  tt = mesh$graph$tt
  tti = mesh$graph$tti
  vv = mesh$graph$vv # sparse matrix, 1 if connected
  
  nV = nrow(corners)
  nT = nrow(tv)
  
  inds = vt
  
  # for each triangle, calculate its center (for calculating spatial 
  # covariances)
  centers = matrix(nrow=nT, ncol=2)
  triCorners = list()
  for(ti in 1:nT) {
    vInds = tv[ti,]
    thisCoords = corners[vInds,]
    
    centers[ti,] = colMeans(thisCoords)
    triCorners = c(triCorners, list(thisCoords))
  }
  
  # get depths at all points
  allDepths = getPointDepths(rbind(centers, corners), maxDepth=maxDepth, method=method)
  centerDepths = allDepths[1:nrow(centers)]
  allCornerDepths = allDepths[-(1:nrow(centers))]
  
  # add depths to coordinates
  centers = data.frame(list(lon=centers[,1], lat=centers[,2], depth=centerDepths))
  goodTri = rep(TRUE, nrow(centers)) # good if has some dip
  internalI = rep(TRUE, nrow(centers)) # all corners within extent
  for(i in 1:nT) {
    vInds = tv[i,]
    thisDepths = allCornerDepths[vInds]
    
    if(any(!is.finite(thisDepths)) || all(thisDepths == thisDepths[1])) {
      goodTri[i] = FALSE
    }
    if(any(!fields::in.poly(triCorners[[i]][,1:2], extent))) {
      internalI[i] = FALSE
    }
    triCorners[[i]] = data.frame(lon=triCorners[[i]][,1], lat=triCorners[[i]][,2], 
                                 depth=thisDepths)
  }
  
  # determine which triangles are in the fault extent
  # internalI = fields::in.poly(centers, extent)
  
  internalCenters = centers[internalI,]
  externalCenters = centers[!internalI,]
  internalTriCorners = list()
  externalTriCorners = list()
  for(i in 1:nrow(centers)) {
    thisTriCorners = triCorners[[i]]
    if(internalI[i] && goodTri[i]) {
      internalTriCorners = c(internalTriCorners, list(thisTriCorners))
    } else {
      externalTriCorners = c(externalTriCorners, list(thisTriCorners))
    }
  }
  # for(i in 1:nrow(centers)) {
  #   thisTriCorners = triCorners[[i]]
  #   if(internalI[i]) {
  #     internalTriCorners = c(internalTriCorners, list(thisTriCorners))
  #   } else {
  #     externalTriCorners = c(externalTriCorners, list(thisTriCorners))
  #   }
  # }
  
  if(FALSE) {
    centers = internalCenters
    corners = internalTriCorners
    
    weirdCenters = centers[internalI & !goodTri,]
    weirdCorners = triCorners[which(internalI & !goodTri)]
    
    for(i in 1:length(weirdCorners)) {
      thisCorners = weirdCorners[[i]]
      
      inExt = fields::in.poly(as.matrix(thisCorners), extent, inflation=-1e-06)
      if(all(inExt)) {
        print(i)
      }
    }
    
    plot(weirdCenters[,1], weirdCenters[,2], pch=19)
    polygon(extent)
    
    plotPolyDat(corners, asp=1)
    plotPolyDat(weirdCorners, col="red", border=rgb(1,0,0), new=FALSE)
    
    getDepths = function(x) {mean(x[,3])}
    plotPolyDat(corners, sapply(corners, getDepths), asp=1, borders=rgb(1,1,1,0), 
                leaveRoomForLegend=TRUE)
    plotPolyDat(weirdCorners, sapply(weirdCorners, getDepths), border=rgb(1,0,0,0), 
                asp=1, new=FALSE, leaveRoomForLegend=TRUE)
  }
  
  list(centers=internalCenters, corners=internalTriCorners,
       externalTriangulation=list(centers=externalCenters, corners=externalTriCorners))
}

# simple nearest neighbor algorithm
# assume pts is in utm10
# Inputs:
# pts: nx2 matrix of longitudes and latitudes at which to calculate fault depths
# maxDepth: maximum depth of fault geometry points used to interpolate between. 
#           should be slightly larger (e.g. 1km) larger than the depths of pts 
#           so calculated depths at pts is interpolated rather than extrapolated
# method: either linear interpolation (default), nearest neighbor estimate, or Gaussian kernel smoother.
#         Just use a linear interpolator here, others were for experiments...
# # res: if method=="kernel", resolution of grid approximating points
getPointDepths = function(pts, maxDepth=Inf, method=c("linear", "NN", "kernel"), res=1) {
  method = match.arg(method)
  
  # first load in the Slab 2.0 geometry
  slab = loadSlab2()
  slab$lon = slab$lon - 360
  lonLat = cbind(slab$lon, slab$lat)
  depths = slab$depth
  
  # get projected Slab 2.0 geometry coordinates
  xy = projCSZ(lonLat)
  
  # keep only points with appropriate depth
  goodI = abs(depths) <= (maxDepth + 1)
  xy = xy[goodI,]
  depths = depths[goodI]
  
  # Things to try:
  #   compactly supported covariance function
  #   kernel smoother
  #   construct my own linear smoother
  
  if(method == "kernel") {
    # require(fdapace)
    
    # calculate minimum bandwidth for 2 nearest neighbors:
    # distMat = rdist(pts, xy)
    # sortedMat = t(apply(distMat, 1, sort))
    # minDists = apply(sortedMat, 2, min)
    # minDists[1:3]
    # [1] 0.02807917 1.84774971 3.32015467
    
    # round coordinates of pts to grid. res resolution or better
    lx = ceiling(diff(range(pts[,1]))/res)
    ly = ceiling(diff(range(pts[,2]))/res)
    xgrid = seq(min(pts[,1]), max(pts[,1]), l=lx)
    ygrid = seq(min(pts[,2]), max(pts[,2]), l=ly)
    xinds = roundToGrid(pts[,1], xgrid, returnInds=TRUE)
    yinds = roundToGrid(pts[,2], ygrid, returnInds=TRUE)
    
    # kernel smoother (Epanechnikov/parabolic kernel)
    # bw=5km
    # totTime = system.time(out <- Lwls2D(bw=5, kern="quart", xin=xy, yin=depths, 
    #                                     xout=pts))[3]
    # totTime = system.time(out <- Lwls2D(bw=4, kern="epan", xin=xy, yin=depths, 
    #                                     xout=pts, crosscov=TRUE))[3]
    # totTime = system.time(out <- Lwls2D(bw=5, kern="epan", xin=xy, yin=depths, 
    #                                     xout1=xgrid, xout2=ygrid))[3]
    # totTime = system.time(out <- Lwls2D(bw=4, kern="gauss", xin=xy, yin=depths,
    #                                     xout1=xgrid, xout2=ygrid, crosscov=TRUE))[3]
    # ptsDepths = out[cbind(xinds, yinds)]
    # totTime/60
    # 6.678383
    
    # smooth.2d(Y, ind = NULL, weight.obj = NULL, setup = FALSE, grid = NULL,
    #           x = NULL, nrow = 64, ncol = 64, surface = TRUE, cov.function =
    #             gauss.cov, Mwidth = NULL, Nwidth = NULL, ...)
    
    totTime = system.time(out <- smooth.2d(Y=depths, x=xy, grid=list(x=xgrid, y=ygrid), 
                                           aRange=4))[3]
    # totTime/60
    # 0.1844333 # for res=1
    # 0.61845 # for res=.5
    # 3.374683 # for res=.25 (but actually more like 2 minutes on laptop...)
    ptsDepths = out$z[cbind(xinds, yinds)]
    
    if(FALSE) {
      
      plotWithColor(pts[,1], pts[,2], ptsDepths, zlim=c(-40, 0), 
                    asp=1, forceColorsInRange=TRUE, cex=.2, pch=19, 
                    xlim=range(xgrid), ylim=range(ygrid))
      plotWithColor(xy[,1], xy[,2], depths, zlim=c(-40, 0), 
                    asp=1, forceColorsInRange=TRUE, cex=.2, pch=19, 
                    xlim=range(xgrid), ylim=range(ygrid))
    }
    
    return(ptsDepths)
  } else if(method == "NN") {
    distMat = rdist(pts, xy)
    nearestInds = apply(distMat, 1, which.min)
    return(depths[nearestInds])
  } else if(method == "linear") {
    totTime = system.time(out <- interp::interp(x=xy[,1], y=xy[,2], z=depths,
                                      xo=pts[,1], yo=pts[,2], output="points",
                                      method="linear"))[3]
    # totTime
    # 2.085 
    
    return(out$z)
  }
}


# old fault geometry code ----
# no longer used except for rectangular faults

# rows: table where each row is a rectangular subfault
# nDown: num division in dip direction
# nStrike: num divisions in strike direction
divideFault = function(rows, nDown=3, nStrike=4) {
  subFaults = apply(rows, 1, divideSubfault, nDown=nDown, nStrike=nStrike)
  do.call("rbind", subFaults)
}

# row: row of a table. Described a rectangular subfault
# nDown: num division in dip direction
# nStrike: num divisions in strike direction
divideSubfault = function(row, nDown=3, nStrike=4) {
  if(!is.list(row))
    row = as.list(row)
  
  geom = calcGeom(row)
  corners = geom$corners
  
  # generate set of nDown points mid strike and going down dip with the correct depths
  centers = matrix(nrow=nDown, ncol=3) # rows are 1,2,...,nDown, and cols are lon,lat,depth
  
  # Simple conversion factors
  #lat2meter = util.dist_latlong2meters(0.0, 1.0)[1]
  #LAT2METER = 110.574 #* 10^3
  LAT2METER = 111133.84012073894 #/10^3
  lat2meter = LAT2METER
  DEG2RAD = 2*pi/360
  
  # Set depths
  centers[,3] = row$depth + (0:(nDown-1))/nDown * row$width * sin(row$dip * DEG2RAD)
  
  # Vector *up_dip* goes from bottom edge to top edge, in meters,
  # from point 2 to point 0 in the figure in the class docstring.
  # Vector *up_strike* goes along the top edge from point d to point a
  # in the figure in the class docstring. (this is different from in calcGeom)
  # up_depth is the depth difference from top to bottom of fault
  up_dip = c(-row$width * cos(row$dip * DEG2RAD) * cos(row$strike * DEG2RAD) 
             / (LAT2METER * cos(row$latitude * DEG2RAD)),
             row$width * cos(row$dip * DEG2RAD) 
             * sin(row$strike * DEG2RAD) / LAT2METER)
  up_strike = c(row$length * sin(row$strike * DEG2RAD) 
                / (lat2meter * cos(geom$centers[3,2] * DEG2RAD)),
                row$length * cos(row$strike * DEG2RAD) / lat2meter)
  
  # Set lon and lat of centers
  centers[,1] = row$longitude - (0:(nDown-1))/nDown * up_dip[1]
  centers[,2] = row$latitude - (0:(nDown-1))/nDown * up_dip[2]
  
  # get points along down strike edge of subfault by moving in down strike direction 
  # from centers
  starts = centers
  starts[,1:2] = starts[,1:2] + (- 0.5 + 1/(2*nStrike))*matrix(rep(up_strike, nrow(centers)), ncol=2, byrow = TRUE)
  
  # compute complete set of subfault coordinates by shifting starts by various vectors in
  # the up strike direction
  coords = cbind(rep(starts[,1], nStrike), rep(starts[,2], nStrike), rep(starts[,3], nStrike))
  shifts = matrix(rep(up_strike, nStrike*nDown), ncol=2, byrow = TRUE)
  mult = rep((0:(nStrike-1))/nStrike, rep(nDown, nStrike))
  shifts = sweep(shifts, MARGIN=1, STATS = mult, FUN = "*")
  coords[,1:2] = coords[,1:2] + shifts
  
  # construct subfaults
  rows = data.frame(matrix(unlist(rep(row, nStrike*nDown)), ncol=length(row), byrow = TRUE))
  names(rows) = names(row)
  rows$Fault = row$Fault
  rows$length = rows$length/nStrike
  rows$width = rows$width/nDown
  rows$depth = coords[,3]
  rows$longitude = coords[,1]
  rows$latitude = coords[,2]
  
  return(rows)
}




