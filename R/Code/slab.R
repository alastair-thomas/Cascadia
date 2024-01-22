
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
getFullFaultGeom = function(triangulatedGeom=NULL, n=2000, max.n=-1, max.edge=c(15, 1000), cutoff=3, 
                            maxDepth=30, ...) {
  
  # construct the triangulated geometry if need be
  if(is.null(triangulatedGeom)) {
    discr = discretizeSlab2(n=n, max.n=max.n, max.edge=max.edge, cutoff=cutoff, 
                                       maxDepth=maxDepth, ...)
    
    triangulatedGeom = discr$geom
    
    boundary = discr$extent
    
    # find the number of subfaults
    nSubFaults = length(triangulatedGeom$corners)
    
    # convert all coordinates back to lon/lat.
    for(i in 1:nSubFaults) {
      triangulatedGeom$corners[[i]][,1:2] = projCSZ(as.matrix(triangulatedGeom$corners[[i]][,1:2]), inverse=TRUE, units="km")
      triangulatedGeom$centers[[i]][,1:2] = projCSZ(as.matrix(triangulatedGeom$centers[[i]][,1:2]), inverse=TRUE, units="km")
    }
  }
  
  cornerList = triangulatedGeom$corners
  centerList = triangulatedGeom$centers
  
  # now convert to the format used by okadaTri
  faultGeom = list()
  
  # loop over each subfault
  for(i in 1:nSubFaults) {
    
    # get a 3x3 matrix, each row is a subfault corner with lon, lat, depth
    thisCorners = as.matrix(cornerList[[i]])
    thisCorners[,3] = thisCorners[,3]*(10^3) # depth from km to m
    
    thisCenter = centerList[[i]]
    thisCenter[3] = thisCenter[3]*(10^3)
    
    # calculate the geometry of the subfault
    triGeom = calculate_geometry_triangles(cbind(projCSZ(thisCorners[,1:2], units="m"), thisCorners[,3]))
    
    # reset already known values
    triGeom$corners = thisCorners
    
    # I changed this to already calculated centers
    # not sure if should use output from calculate_geometry_triangles
    # lon, lat are the same but depth slightly different
    
    triGeom$lon = as.numeric(thisCenter[1])
    triGeom$lat = as.numeric(thisCenter[2])
    triGeom$depth = as.numeric(thisCenter[3])
    
    # append to list of subfaults
    faultGeom = c(faultGeom, list(triGeom))
  }
  
  # return results
  return(faultGeom)
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
discretizeSlab2 = function(n=2000, max.n=-1, max.edge=c(15, 1000), maxDepth=30, 
                           cutoff=3, method=c("linear", "NN", "kernel"), ...) {
  method = match.arg(method)
  
  # first load in the Slab 2.0 geometry
  slab = loadSlab2()
  lonLat = cbind(slab$lon - 360, slab$lat) # maybe need a -360 here
  depths = slab$depth
  
  # get projected coordinates
  xy = projCSZ(lonLat, units="km")
  
  # keep only points with appropriate depth
  goodI = abs(depths) <= maxDepth
  lonLat = lonLat[goodI,]
  xy = xy[goodI,]
  depths = depths[goodI]
  
  # create the interior hull
  require(concaveman)
  
  # the concavity and length_threshold are set at default values.
  # don't think they need changed.
  concaveHull = concaveman(xy, concavity=2, length_threshold=0)
  
  # this sections slightly shrinks the hull
  # this guarantees that all the points for the subfaults will have a depth
  # otherwise some depths become NA when depth is interpolated form Slab2
  
  v = vect(concaveHull, "polygons") # convert to format needed for buffer function
  v = buffer(v, -0.1) # shrink the hull by 0.1km
  v = as.points(v) # convert to points
  concaveHull = crds(v) # get the hull in same format as before
  
  #plot(projCSZ(xy, unit="km", inverse=TRUE))
  #plot(projCSZ(concaveHull, units="km", inverse=TRUE))
  
  # testing is.bnd=TRUE since this is a boundary?
  # The shape being made is boundary points, thus I think is.bnd=TRUE is fine?
  concaveInt = inla.mesh.segment(concaveHull, is.bnd=TRUE)
  
  # make sure concaveInt is in a format expected by inla.mesh.2d
  concaveInt$idx = rbind(concaveInt$idx, c(nrow(concaveInt$loc), 1))
  concaveInt$idx = matrix(as.integer(concaveInt$idx), ncol=2)
  concaveInt$grp = matrix(rep(as.integer(1), nrow(concaveInt$loc)), ncol=1)
  concaveInt$loc = matrix(concaveInt$loc[,1:2], ncol=2) # This line got changed
  concaveInt$loc = concaveInt$loc[nrow(concaveInt$loc):1, ]
  
  # create the exterior hull
  hullExt = inla.nonconvex.hull.basic(xy, resolution=150, convex=-.4)
  
  # construct mesh with INLA
  mesh = inla.mesh.2d(n=n, boundary=list(concaveInt, hullExt), max.n=max.n, 
                      max.edge=max.edge, cutoff=cutoff, ...)
  
  # The geometry produced here is a bit weird
  # number of slabs doesn't equal number of centroids!
  # I fixed it, mismatch in deciding which slabs were interior/exterior
   
  faultGeom = getGeomFromMesh(mesh, extent=concaveHull, maxDepth=maxDepth, 
                              method=method)
  
  
  return(list(mesh=mesh, geom=faultGeom, extent=projCSZ(concaveHull, units="km", inverse=TRUE), maxDepth=maxDepth))
}

getTrianglesFromMesh = function(mesh) {
  
  corners = mesh$loc[,1:2] # corners of the triangles, i.e. vertices
  
  # t: triangles, v: vertices
  tv = mesh$graph$tv
  vt = mesh$graph$vt
  tt = mesh$graph$tt
  tti = mesh$graph$tti
  vv = mesh$graph$vv # sparse matrix, 1 if connected
  
  nV = nrow(corners) # number of points in mesh
  nT = nrow(tv) # number of triangles in mesh
  
  # for each triangle, calculate its center
  centers = matrix(nrow=nT, ncol=2) # (number of subfaults, 2)
  triCorners = list()
  for(ti in 1:nT) { # loop over each subfault
    vInds = tv[ti,] # get the indicies for corners of this subfault
    thisCoords = corners[vInds,] # get the corners of the subfault
    
    centers[ti,] = colMeans(thisCoords) # center is just the mean of Lon, Lat, Depth
    triCorners = c(triCorners, list(thisCoords))
  }
  
  return(triCorners)
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
  
  nV = nrow(corners) # number of points in mesh
  nT = nrow(tv) # number of triangles in mesh
  
  # for each triangle, calculate its center
  centers = matrix(nrow=nT, ncol=2) # (number of subfaults, 2)
  triCorners = list()
  for(ti in 1:nT) { # loop over each subfault
    vInds = tv[ti,] # get the indicies for corners of this subfault
    thisCoords = corners[vInds,] # get the corners of the subfault
    
    centers[ti,] = colMeans(thisCoords) # center is just the mean of Lon, Lat, Depth
    triCorners = c(triCorners, list(thisCoords))
  }
  
  # centers - (number of subfaults, 2) containing each subfault centroid
  # triCorners - vector containings list of corners
  
  cornerDepths = getPointDepths(corners, maxDepth=30, method=method)
  centerDepths  = getPointDepths(centers, maxDepth=30, method=method)
  
  # add depths to center coordinates
  # centers becomes a dataframe
  centers = data.frame(list(lon=centers[,1], lat=centers[,2], depth=centerDepths))
  
  # indicator variables
  goodTri = rep(TRUE, nT) # good if has some dip
  internalI = rep(TRUE, nT) # all corners within extent
  
  
  # add depth to corners
  for(ti in 1:nT) {
    vInds = tv[ti,] # get each row of corner indicies
    thisDepths = cornerDepths[vInds]
    
    # if depth NA or all depths the same
    if(any(!is.finite(thisDepths)) || all(thisDepths == thisDepths[1])) {
      # print(thisDepths)
      # Everything excluded was because of NA values
      goodTri[ti] = FALSE
    }
    
    # if any vertex lies outside the interior hull
    
    if(any(!fields::in.poly(triCorners[[ti]][,1:2], extent))) {
      internalI[ti] = FALSE
    }
    
    triCorners[[ti]] = data.frame(lon=triCorners[[ti]][,1], lat=triCorners[[ti]][,2], depth=thisDepths)
  }
  
  
  
  # determine which triangles are in the fault extent
  # internalI = fields::in.poly(centers, extent)
  
  # This didn't take into account the goodTri condition
  # created a mis-match between the number of internal triangles and internal centers
  #internalCenters = centers[internalI,]
  #externalCenters = centers[!internalI,]
  
  internalCenters = list()
  externalCenters = list()
  
  internalTriCorners = list()
  externalTriCorners = list()
  
  # loop over each subfault
  # add it as an interior or exterior subfault
  for(ti in 1:nT) {
    
    thisTriCorners = triCorners[[ti]]
    thisCenter = centers[ti,]
    row.names(thisCenter) = NULL
    
    if(internalI[ti]) {
       if (goodTri[ti]){
         internalTriCorners = c(internalTriCorners, list(thisTriCorners))
         internalCenters    = c(internalCenters, list(thisCenter))
       }
      else{
        xy = cbind(thisTriCorners$lon, thisTriCorners$lat)
        xy = projCSZ(xy, inverse=TRUE, units="km")
        
        print(cbind(xy, thisTriCorners$depth))
      }
      
    } else {
      
      externalTriCorners = c(externalTriCorners, list(thisTriCorners))
      externalCenters    = c(externalCenters, list(thisCenter))
    }
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
    
    totTime = system.time(out <- smooth.2d(Y=depths, x=xy, grid=list(x=xgrid, y=ygrid), 
                                           aRange=4))[3]
    # totTime/60
    # 0.1844333 # for res=1
    # 0.61845 # for res=.5
    # 3.374683 # for res=.25 (but actually more like 2 minutes on laptop...)
    
    ptsDepths = out$z[cbind(xinds, yinds)]
    
    return(ptsDepths)
    
  } else if(method == "NN") {
    distMat = rdist(pts, xy)
    nearestInds = apply(distMat, 1, which.min)
    return(depths[nearestInds])
    
  } else if(method == "linear") {
    
    out = interp::interp(x=xy[,1], y=xy[,2], z=depths,
                         xo=pts[,1], yo=pts[,2],
                         output="points", method="linear")
    
    #totTime = system.time(out <- interp::interp(x=xy[,1], y=xy[,2], z=depths,
     #                                 xo=pts[,1], yo=pts[,2], output="points",
     #                                 method="linear"))[3]
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

# AT - a function I found that Finn had made
# Maybe useful for us
# Link: https://groups.google.com/g/r-inla-discussion-group/c/z1n1exlZrKM
#
#' Convert inla.mesh to sp objects
#'
#' @param mesh An \code{\link{inla.mesh}} object
#' @return A list with \code{sp} objects for triangles and vertices:
#' \describe{
#' \item{triangles}{\code{SpatialPolygonsDataFrame} object with the triangles in
#' the same order as in the original mesh, but each triangle looping through
#' the vertices in clockwise order (\code{sp} standard) instead of
#' counterclockwise order (\code{inla.mesh} standard). The \code{data.frame}
#' contains the vertex indices for each triangle, which is needed to link to
#' functions defined on the vertices of the triangulation.
#' \item{vertices}{\code{SpatialPoints} object with the vertex coordinates,
#' in the same order as in the original mesh.}
#' }
#' @export
inla.mesh2sp <- function(mesh) {
  crs <- inla.CRS(inla.CRSargs(mesh$crs))
  isgeocentric <- identical(inla.as.list.CRS(crs)[["proj"]], "geocent")
  if (isgeocentric || (mesh$manifold == "S2")) {
    stop(paste0(
      "'sp' doesn't support storing polygons in geocentric coordinates.\n",
      "Convert to a map projection with inla.spTransform() before
calling inla.mesh2sp()."))
    }

  triangles <- SpatialPolygonsDataFrame(
    Sr = SpatialPolygons(lapply(
      1:nrow(mesh$graph$tv),
      function(x) {
        tv <- mesh$graph$tv[x, , drop = TRUE]
        Polygons(list(Polygon(mesh$loc[tv[c(1, 3, 2, 1)],
                                       1:2,
                                       drop = FALSE])),
                 ID = x)
      }
    ),
    proj4string = crs
    ),
    data = as.data.frame(mesh$graph$tv[, c(1, 3, 2), drop = FALSE]),
    match.ID = FALSE
  )
  vertices <- SpatialPoints(mesh$loc[, 1:2, drop = FALSE], proj4string = crs)
  
  list(triangles = triangles, vertices = vertices)
}





