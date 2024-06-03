# utility functions (miscellaneous but useful)

getMomentFromSlip = function(slips, rigidity=4*10^10, doTaper=FALSE, lambda=1, fault=csz, depths=getFaultCenters(fault)[,3], 
                             normalizeTaper=FALSE, dStar=28000) {
  # get fault total area
  areas = fault$length*fault$width
  totArea = sum(areas)
  
  if(doTaper)
    slips = slips*taper(depths, lambda, dStar=dStar, normalize=normalizeTaper)
  
  Mo = sum(slips*rigidity*areas)
  
  if(Mo <= 0) {
    warning("negative seismic moment observed")
    return(0)
  }
  
  (log10(Mo) - 9.05)/1.5
}

# Projects from lon/lat (EPSG:4326) to utm coordinates in km based on UTM 10 
# (EPSG:32610) or back if inverse==TRUE.
# x: a matrix of points, where each row is the coordinates of a single point
# inverse: if FALSE, converts from lon/lat to utm, else does the reverse
# units: either "km" (the default) or "m". Sets the input units if 
#        inverse==TRUE or the output units otherwise
projCSZ = function(x, inverse=FALSE, units=c("km", "m")) {
  units = match.arg(units)
  
  # determine the "to" projection systems
  utmProj = st_crs("EPSG:32610")
  lonLatProj = st_crs("EPSG:4326")
  
  if(!inverse) {
    toProj = utmProj
    fromProj = lonLatProj
  } else {
    toProj = lonLatProj
    fromProj = utmProj
  }
  
  # convert to an sf object
  x = st_multipoint(x, dim="XY")
  if(inverse && (units == "km")) {
    # make sure we convert from km back to m as the projection is expecting
    x = x*1000
  }
  x = st_sfc(x, crs=fromProj)
  
  # we already know it is utm zone 10
  if(FALSE) {
    # calculate UTM zone based on slab 2 geometry
    slab = loadSlab2()
    lonRange = range(slab$lon)
    latRange = range(slab$lat)
    midLon = mean(lonRange)
    
    utmZone = ceiling((midLon-180)/6) # UTM zone 10
  }
  
  # transform coordinates and convert back to a matrix
  out = st_transform(x, toProj)
  out = st_coordinates(out)[,1:2]
  
  # convert from m to km if necessary
  if((toProj == utmProj) && (units == "km")) {
    out = out/1000
  }
  
  out
}

# same as projCSZ, but under old projection system used by geoclaw
projCSZ2 = function(x, inverse=FALSE, units=c("km", "m")) {
  units = match.arg(units)
  
  DEG2RAD = 2*pi/360
  LAT2METER = 111133.84012073894 #/10^3
  if(units == "km") {
    LAT2METER = LAT2METER / 10^3
  }
  
  # do the projection
  if(!inverse) {
    # from lon/lat to m
    
    LON2METER = LAT2METER * cos( DEG2RAD*x[,2] ) 
    
    LONLAT2METER = cbind(LON2METER, LAT2METER)
    scales = LONLAT2METER
  } else {
    # from m to lon/lat
    lats = x[,2] * (1/LAT2METER)
    LON2METER = LAT2METER * cos( DEG2RAD*lats ) 
    
    LONLAT2METER = cbind(LON2METER, LAT2METER)
    scales = 1/LONLAT2METER
  }
  
  # return scaled coordinates
  x * scales
}

# round vector of coordinates to the given grid of equally spaced points
roundToGrid = function(coordVec, coordGrid, returnInds=FALSE) {
  inds = (coordVec - min(coordGrid))/(max(coordGrid) - min(coordGrid))
  inds = round(inds*(length(coordGrid)-1)) + 1
  
  if(returnInds) {
    return(inds)
  } else {
    return(coordGrid[inds])
  }
}


