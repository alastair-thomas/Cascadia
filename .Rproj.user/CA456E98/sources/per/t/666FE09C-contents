# Util functions from Python ----

# function to create grids of 2d arrays forming mesh based on 1d vectors
# (output format should be the same as Python's numpy.meshgrid)
meshgrid = function(x, y, transpose=FALSE) {
  if(transpose) {
    X = matrix(rep(x, length(y)), nrow=length(x))
    Y = matrix(rep(y, length(x)), ncol=length(y), byrow = TRUE)
  }
  else { # numpy format
    X = matrix(rep(x, length(y)), nrow=length(y), byrow=TRUE)
    Y = matrix(rep(y, length(x)), ncol=length(x))
  }
  return(list(X=X, Y=Y))
}

# Rectangles ----

okadaRect = function(rows, x, y, slips=rep(0, nrow(rows)), rakes=rep(90, nrow(rows)), 
                     poisson=0.25, inds=NULL) {
  # allow user to send a constant for slips and rakes instead of vector for ease of use
  if(length(slips) == 1)
    slips = rep(slips, nrow(rows))
  if(length(rakes) == 1)
    rakes = rep(rakes, nrow(rows))
  
  fullOkada = matrix(0, nrow=length(y), ncol=length(x))
  for(i in 1:nrow(rows)) {
    fullOkada = fullOkada + okadaSubfaultRect(rows[i,], x, y, slips[i], rakes[i], poisson, inds)
  }
  return(fullOkada)
}

# same as Okada, but return a matrix of dimension nrow(datCoords) x nrow(rows)
# with a decomposition of the seaDef at each data location induced by the slip 
# at each subfault (G matrix from model)
okadaAllRect = function(rows, lonGrid, latGrid, datCoords, slips=rep(0, nrow(rows)), rakes=rep(90, nrow(rows)), 
                        poisson=0.25) {
  # allow user to send a constant for slips and rakes instead of vector for ease of use
  if(length(slips) == 1)
    slips = rep(slips, nrow(rows))
  if(length(rakes) == 1)
    rakes = rep(rakes, nrow(rows))
  
  # calculate the simulated subsidence at the data locations
  # round the data locations to the lon lat grid to make this function faster
  # and to have consistent grid across different subfaults
  roundToRange = function(coordVec, coordRange) {
    inds = (coordVec - min(coordRange))/(max(coordRange) - min(coordRange))
    inds = round(inds*(length(coordRange)-1)) + 1
    roundedCoordVec = coordRange[inds]
    return(roundedCoordVec)
  }
  roundLon = roundToRange(datCoords[,1], lonGrid)
  roundLat = roundToRange(datCoords[,2], latGrid)
  roundCoords = cbind(roundLon, roundLat)
  
  # find indices of grid coords corresponding to rounded data coords
  findIndex = function(rCoords, gCoords) {
    return(match(data.frame(t(rCoords)), data.frame(t(gCoords))))
  }
  coordGrid = make.surface.grid(list(lonGrid, latGrid))
  inds = findIndex(roundCoords, coordGrid)
  
  # compute the result of the Okada model for each subfault at each data location
  fullOkada = matrix(0, nrow=nrow(datCoords), ncol=nrow(rows))
  for(i in 1:nrow(rows)) {
    # simulate the seaDef
    # dZ = okadaSubfaultRect(rows[i,], x, y, slips[i], rakes[i], poisson)
    fullOkada[,i] = okadaSubfaultRect(rows[i,], lonGrid, latGrid, slips[i], rakes[i], poisson, inds)
    
    # get the simulated seaDef at the data locations
    # fullOkada[i,] = c(t(dZ))[inds]
  }
  return(fullOkada)
}

# function based on dtopotools.Subfault.okada function in python for rectangular 
# fault geometries.
#
# x and y are 1D vectors forming the coordinates of a grid (as in 
# make.surface.grid(list(x, y))   )
# ind is a set of indices such that c(t(dZ))[inds] gives the desired seaDefs
# where dZ is the result of okadaSubfault
okadaSubfaultRect = function(row, x, y, slip=0, rake=90, poisson=0.25, inds=NULL) {
  # Apply Okada to this subfault and return a DTopography object.
  # 
  # :Input:
  # - x,y are 1d arrays
  # :Output:
  # - DTopography object with dZ array of shape (1,len(x),len(y))
  # with single static displacement and times = [0.].
  # 
  # Currently only calculates the vertical displacement.
  # 
  # Okada model is a mapping from several fault parameters
  # to a surface deformation.
  # See Okada 1985 [Okada85]_, or Okada 1992, Bull. Seism. Soc. Am.
  # 
  # okadamap function riginally written in Python by Dave George for
  # Clawpack 4.6 okada.py routine, with some routines adapted
  # from fortran routines written by Xiaoming Wang.
  # 
  # Rewritten and made more flexible by Randy LeVeque
  # 
  # **Note:** *self.coordinate_specification* (str) specifies the location on 
  # each subfault that corresponds to the (longitude,latitude) and depth 
  # subfault.
  # 
  # See the documentation for *SubFault.calculate_geometry* for dicussion of the 
  # possible values *self.coordinate_specification* can take.
  
  if(!is.list(row))
    row = as.list(row)
  
  # Setup coordinate arrays
  #   corners = [[None, None, None], # a 
  #              [None, None, None], # b
  #              [None, None, None], # c
  #              [None, None, None]] # d
  #   centers = [[None, None, None], # 1
  #              [None, None, None], # 2 
  #              [None, None, None]] # 3
  geom = calcGeom(row)
  centers = geom$centers
  corners = geom$corners
  
  # Okada model assumes x,y are at bottom center:
  x_bottom = centers[3,1]
  y_bottom = centers[3,2]
  depth_bottom = centers[3,3]
  
  length = row$length
  width = row$width
  depth = row$depth
  
  halfL = 0.5*length
  w  =  width
  
  # convert angles to radians:
  DEG2RAD = 2*pi/360
  ang_dip = DEG2RAD * row$dip
  ang_rake = DEG2RAD * rake
  ang_strike = DEG2RAD * row$strike
  
  # this format should be the same as numpy.meshgrid
  mesh = meshgrid(x, y)   # use convention of upper case for 2d
  X = mesh$X
  Y = mesh$Y
  
  # Convert distance from (X,Y) to (x_bottom,y_bottom) from degrees to
  # meters:
  LAT2METER = 111133.84012073894 #/10^3
  xx = LAT2METER * cos(DEG2RAD * Y) * (X - x_bottom)   
  yy = LAT2METER * (Y - y_bottom)
  
  # if user only wants a subset of seaDefs given by inds, subset now to save time
  if(!is.null(inds)) {
    xx = c(t(xx))[inds]
    yy = c(t(yy))[inds]
  }
  
  # Convert to distance along strike (x1) and dip (x2):
  x1 = xx * sin(ang_strike) + yy * cos(ang_strike) 
  x2 = xx * cos(ang_strike) - yy * sin(ang_strike) 
  
  # In Okada's paper, x2 is distance up the fault plane, not down dip:
  x2 = -x2
  
  p = x2 * cos(ang_dip) + depth_bottom * sin(ang_dip)
  q = x2 * sin(ang_dip) - depth_bottom * cos(ang_dip)
  
  # to save computation time, set strike slip to 0 if rake = 90
  if(rake != 90) {
    f1 = strike_slip(x1 + halfL, p,     ang_dip, q, poisson)
    f2 = strike_slip(x1 + halfL, p - w, ang_dip, q, poisson)
    f3 = strike_slip(x1 - halfL, p,     ang_dip, q, poisson)
    f4 = strike_slip(x1 - halfL, p - w, ang_dip, q, poisson)
  }
  else {
    f1 = f2 = f3 = f4 = rep(0, length(x1))
  }
  
  g1=dip_slip(x1 + halfL, p,     ang_dip, q, poisson)
  g2=dip_slip(x1 + halfL, p - w, ang_dip, q, poisson)
  g3=dip_slip(x1 - halfL, p,     ang_dip, q, poisson)
  g4=dip_slip(x1 - halfL, p - w, ang_dip, q, poisson)
  
  # Displacement in direction of strike and dip:
  ds = slip * cos(ang_rake)
  dd = slip * sin(ang_rake)
  
  us = (f1 - f2 - f3 + f4) * ds
  ud = (g1 - g2 - g3 + g4) * dd
  
  dz = (us+ud)
  
  # I opted to only return the dZ instead of the list of objects for simplicity
  #dtopo = list()
  #dtopo$X = X
  #dtopo$Y = Y
  #dtopo$dZ = dz
  #dtopo$times = 0.
  return(dz)
}

strike_slip = function(y1, y2, ang_dip, q, poisson=0.25) {
  # Used for Okada's model
  # Methods from Yoshimitsu Okada (1985)
  
  sn = sin(ang_dip)
  cs = cos(ang_dip)
  d_bar = y2*sn - q*cs
  r = sqrt(y1^2 + y2^2 + q^2)
  xx = sqrt(y1^2 + q^2)
  a4 = 2.0*poisson/cs*(log(r+d_bar) - sn*log(r+y2))
  # rewritten to use multiplication instead of division when possible (faster)
  # f = -(d_bar*q/r/(r+y2) + q*sn/(r+y2) + a4*sn)/(2.0*pi)
  f = -(d_bar*q/(r*(r+y2)) + q*sn/(r+y2) + a4*sn)*(1/(2.0*pi))
  
  return(f)
}

dip_slip = function(y1, y2, ang_dip, q, poisson=0.25) {
  # Based on dtopotools.SubFault._strike_slip and Okada's paper (1985)
  # Added by Xiaoming Wang
  
  sn = sin(ang_dip)
  cs = cos(ang_dip)
  
  d_bar = y2*sn - q*cs;
  r = sqrt(y1^2 + y2^2 + q^2)
  xx = sqrt(y1^2 + q^2)
  #   a5 = 4.*poisson/cs*atan((y2*(xx+q*cs)+xx*(r+xx)*sn)/y1/(r+xx)/cs)
  #   f = -(d_bar*q/r/(r+y1) + sn*atan(y1*y2/q/r) - a5*sn*cs)/(2.0*pi)
  # rewritten to use multiplication instead of division when possible (faster)
  a5 = 4.*poisson/cs*atan((y2*(xx+q*cs)+xx*(r+xx)*sn)/(y1*(r+xx)*cs))
  f = -(d_bar*q/(r*(r+y1)) + sn*atan(y1*y2/(q*r)) - a5*sn*cs)*(1/(2.0*pi))
  
  return(f)
}

calcGeom = function(subfault) {
  #   Calculate the fault geometry.
  #   
  #   Routine calculates the class attributes *corners* and 
  #   *centers* which are the corners of the fault plane and 
  #   points along the centerline respecitvely in 3D space.
  #   
  #   **Note:** *self.coordinate_specification*  specifies the location on each
  #   subfault that corresponds to the (longitude,latitude) and depth 
  #   of the subfault.
  #   Currently must be one of these strings:
  #     
  #   - "bottom center": (longitude,latitude) and depth at bottom center
  #   - "top center": (longitude,latitude) and depth at top center
  #   - "centroid": (longitude,latitude) and depth at centroid of plane
  #   - "noaa sift": (longitude,latitude) at bottom center, depth at top,  
  #   This mixed convention is used by the NOAA SIFT
  #   database and "unit sources", see:
  #     http://nctr.pmel.noaa.gov/propagation-database.html
  #   
  #   The Okada model is expressed assuming (longitude,latitude) and depth
  #   are at the bottom center of the fault plane, so values must be
  #   shifted or other specifications.
  
  row = subfault
  if(!is.list(row))
    row = data.frame(row)
  
  # Simple conversion factors
  #lat2meter = util.dist_latlong2meters(0.0, 1.0)[1]
  #LAT2METER = 110.574 #* 10^3
  LAT2METER = 111133.84012073894 #/10^3
  lat2meter = LAT2METER
  DEG2RAD = 2*pi/360
  
  #   Setup coordinate arrays
  #   Python format:
  #   Top edge    Bottom edge
  #     a ----------- b          ^ 
  #     |             |          |         ^
  #     |             |          |         |
  #     |             |          |         | along-strike direction
  #     |             |          |         |
  #     0------1------2          | length  |
  #     |             |          |
  #     |             |          |
  #     |             |          |
  #     |             |          |
  #     d ----------- c          v
  #     <------------->
  #     width
  #   
  #     <-- up dip direction
  
  #   corners = [[None, None, None], # a 
  #              [None, None, None], # b
  #              [None, None, None], # c
  #              [None, None, None]] # d
  #   centers = [[None, None, None], # 1
  #              [None, None, None], # 2 
  #              [None, None, None]] # 3
  corners = matrix(nrow=4, ncol=3) # rows are a,b,c,d, and cols are lon,lat,depth
  centers = matrix(nrow=3, ncol=3) # rows are 1,2,3 (0,1,2), and cols are lon,lat,depth
  
  # Set depths
  centers[1,3] = row$depth
  centers[2,3] = row$depth + 0.5 * row$width * sin(row$dip * DEG2RAD)
  centers[3,3] = row$depth + row$width * sin(row$dip * DEG2RAD)
  
  corners[1,3] = centers[1,3]
  corners[4,3] = centers[1,3]
  corners[2,3] = centers[3,3]
  corners[3,3] = centers[3,3]
  
  # Locate fault plane in 3D space:  
  # See the class docstring for a guide to labeling of corners/centers.
  
  # Vector *up_dip* goes from bottom edge to top edge, in meters,
  # from point 2 to point 0 in the figure in the class docstring.
  
  up_dip = c(-row$width * cos(row$dip * DEG2RAD) * cos(row$strike * DEG2RAD) 
             / (LAT2METER * cos(row$latitude * DEG2RAD)),
             row$width * cos(row$dip * DEG2RAD) 
             * sin(row$strike * DEG2RAD) / LAT2METER)
  
  centers[1,1:2] = c(row$longitude, row$latitude)
  centers[2,1:2] = c(row$longitude - 0.5 * up_dip[1],
                     row$latitude - 0.5 * up_dip[2])
  centers[3,1:2] = c(row$longitude - up_dip[1],
                     row$latitude - up_dip[2])
  
  # Calculate coordinates of corners:
  # Vector *strike* goes along the top edge from point 0 to point a
  # in the figure in the class docstring.
  
  up_strike = c(0.5 * row$length * sin(row$strike * DEG2RAD) 
                / (lat2meter * cos(centers[3,2] * DEG2RAD)),
                0.5 * row$length * cos(row$strike * DEG2RAD) / lat2meter)
  
  corners[1,1:2] = c(centers[1,1] + up_strike[1],
                     centers[1,2] + up_strike[2])
  corners[2,1:2] = c(centers[3,1] + up_strike[1], 
                     centers[3,2] + up_strike[2])
  corners[3,1:2] = c(centers[3,1] - up_strike[1],
                     centers[3,2] - up_strike[2])
  corners[4,1:2] = c(centers[1,1] - up_strike[1],
                     centers[1,2] - up_strike[2])
  
  return(list(corners=corners, centers=centers))
}

# Triangles ----

# Apply Okada to this subfault and return a DTopography object.
# 
# :Input:
#   - x,y are 1d arrays
# :Output:
#   - DTopography object with dZ array of shape (1,len(x),len(y))
# with single static displacement and times = [0.].
# 
# Currently only calculates the vertical displacement.
# 
# Okada model is a mapping from several fault parameters
# to a surface deformation.
# See Okada 1985 [Okada85]_, or Okada 1992, Bull. Seism. Soc. Am.
# 
# okadamap function riginally written in Python by Dave George for
# Clawpack 4.6 okada.py routine, with some routines adapted
# from fortran routines written by Xiaoming Wang.
# 
# Rewritten and made more flexible by Randy LeVeque
# 
# **Note:** *self.coordinate_specification* (str) specifies the location on 
# each subfault that corresponds to the (longitude,latitude) and depth 
# subfault.
# 
# See the documentation for *SubFault.calculate_geometry* for discussion of the 
# possible values *self.coordinate_specification* can take.

# JP's documentation:
# fault: a list of subfaults. See okadaSubfaultTri for details
# x, y: x, y coordinates of locations at which to calculate surface deformation
# slip: a vector of coseismic slips, 1 for each subfault
# rake: a vector of rakes, 1 for each subfault
okadaTri = function(fault, x, y, 
                    slip=rep(1, length(fault)), 
                    rake=rep(90, length(fault))) {
  
  if(length(slip) == 1) {
    slip = rep(slip, length(fault))
  }
  if(length(rake) == 1) {
    rake = rep(rake, length(fault))
  }
  
  # add slip, rake to subfaults
  for(i in 1:length(fault)) {
    fault[[i]]$rake = rake[i]
    fault[[i]]$slip = slip[i]
  }
  
  # get topographic deformations for all subfaults
  allTopos = lapply(fault, okadaSubfaultTri, x=x, y=y)
  
  # now add them together
  dtopo = allTopos[[1]]
  for(i in 2:length(fault)) {
    thisTopo = allTopos[[i]]
    
    dtopo$dX = dtopo$dX + thisTopo$dX
    dtopo$dY = dtopo$dY + thisTopo$dY
    dtopo$dZ = dtopo$dZ + thisTopo$dZ
  }
  
  dtopo
}

# subfault: a list with elements:
#     fix_orientation: whether the orientation is fixed. Set by 
#                      calculate_geometry_triangles
#     corners: 3x3 matrix where each row is (lon, lat, depth) of a single corner
#     centers: 
#     lon: longitude TODO
#     lat: latitude TODO
#     depth: depth TODO
#     strike: strike TODO
#     dip: dip TODO
#     slip: average coseismic slip for the subfault
#     rake: rake of the subfault
# x, y: x, y coordinates of locations at which to calculate surface deformation
# NOTE: if subfault$rake not set, defaults to 90
okadaSubfaultTri = function(subfault, x, y) {
  
  if(is.null(subfault$rake)) {
    subfault$rake = 90
  }
  if(is.null(subfault$slip)) {
    subfault$slip = 1
  }
  
  out = meshgrid(x, y)   # uppercase
  X1 = out$X
  X2 = out$Y
  # X3 = numpy.zeros(X1.shape)   # depth zero
  X3 = matrix(0, nrow=nrow(X1), ncol=ncol(X1))   # depth zero
  
  # compute burgers vector
  slipv = get_unit_slip_vector(subfault) 
  # browser() # check slipv
  burgersv = slipv * subfault$slip
  
  # get beta angles
  # out = self._get_leg_angles()
  out = get_leg_angles(subfault)
  reverse_list = out$reverse_list
  O1_list = out$O1_list
  O2_list = out$O2_list
  alpha_list = out$alpha_list
  beta_list = out$beta_list
  
  #
  # v11 = numpy.zeros(X1.shape)
  # v21 = numpy.zeros(X1.shape)
  # v31 = numpy.zeros(X1.shape)
  # 
  # v12 = numpy.zeros(X1.shape)
  # v22 = numpy.zeros(X1.shape)
  # v32 = numpy.zeros(X1.shape)
  # 
  # v13 = numpy.zeros(X1.shape)
  # v23 = numpy.zeros(X1.shape)
  # v33 = numpy.zeros(X1.shape)
  dims = dim(X1)
  v11 = matrix(0, nrow=dims[1], ncol=dims[2])
  v21 = matrix(0, nrow=dims[1], ncol=dims[2])
  v31 = matrix(0, nrow=dims[1], ncol=dims[2])
  
  v12 = matrix(0, nrow=dims[1], ncol=dims[2])
  v22 = matrix(0, nrow=dims[1], ncol=dims[2])
  v32 = matrix(0, nrow=dims[1], ncol=dims[2])
  
  v13 = matrix(0, nrow=dims[1], ncol=dims[2])
  v23 = matrix(0, nrow=dims[1], ncol=dims[2])
  v33 = matrix(0, nrow=dims[1], ncol=dims[2])
  
  # for(j in range(6)) {
  for(j in 0:5) {
    # k = j%3
    k = (j%%3)+1
    alpha = alpha_list[[k]]
    beta = beta_list[[k]]
    
    if (floor(j/3) == 0) {
      # Olong = O1_list[k][0]
      # Olat = O1_list[k][1]
      # Odepth = abs(O1_list[k][2])
      Olong = O1_list[[k]][1]
      Olat = O1_list[[k]][2]
      Odepth = abs(O1_list[[k]][3])
    } else if (floor(j/3) == 1) {
      # Olong = O2_list[k][0]
      # Olat = O2_list[k][1]
      # Odepth = abs(O2_list[k][2])
      Olong = O2_list[[k]][1]
      Olat = O2_list[[k]][2]
      Odepth = abs(O2_list[[k]][3])
    }
    
    if(reverse_list[k]) {
      sgn = (-1.)^(floor(j/3))
    } else {
      sgn = (-1.)^floor(j/3+1)
    }
    
    # fix orientation 
    if(subfault$fix_orientation) {
      sgn = -1 * sgn
    }
    
    out = get_halfspace_coords(subfault, X1,X2,X3,alpha,beta,Olong,Olat,Odepth)
    Y1 = out$Y1
    Y2 = out$Y2
    Y3 = out$Y3
    Z1 = out$Z1
    Z2 = out$Z2
    Z3 = out$Z3
    Yb1 = out$Yb1
    Yb2 = out$Yb2
    Yb3 = out$Yb3
    Zb1 = out$Zb1
    Zb2 = out$Zb2
    Zb3 = out$Zb3
    
    out = get_angular_dislocations_surface(subfault, Y1,Y2,Y3,beta,Odepth)
    w11 = out$v11
    w12 = out$v12
    w13 = out$v13
    w21 = out$v21
    w22 = out$v22
    w23 = out$v23
    w31 = out$v31
    w32 = out$v32
    w33 = out$v33
    
    if(any(is.na(w11))) {
      browser()
    }
    
    out = coord_transform(subfault, w11,w12,w13,w21,w22,w23,w31,w32,w33,alpha)
    w11 = out$v11
    w12 = out$v12
    w13 = out$v13
    w21 = out$v21
    w22 = out$v22
    w23 = out$v23
    w31 = out$v31
    w32 = out$v32
    w33 = out$v33
    
    if(any(is.na(w11))) {
      browser()
    }
    
    v11 = v11 + sgn*w11
    v21 = v21 + sgn*w21
    v31 = v31 + sgn*w31
    
    v12 = v12 + sgn*w12
    v22 = v22 + sgn*w22
    v32 = v32 + sgn*w32
    
    v13 = v13 + sgn*w13
    v23 = v23 + sgn*w23
    v33 = v33 + sgn*w33
  }
  
  # linear combination for each component of Burgers vectors
  # dX = -v11*burgersv[0] - v12*burgersv[1] + v13*burgersv[2]
  # dY = -v21*burgersv[0] - v22*burgersv[1] + v23*burgersv[2]
  # dZ = -v31*burgersv[0] - v32*burgersv[1] + v33*burgersv[2]
  dX = -v11*burgersv[1] - v12*burgersv[2] + v13*burgersv[3]
  dY = -v21*burgersv[1] - v22*burgersv[2] + v23*burgersv[3]
  dZ = -v31*burgersv[1] - v32*burgersv[2] + v33*burgersv[3]
  
  dtopo = list()
  dtopo$X = X1
  dtopo$Y = X2
  
  dtopo$dX = dX
  dtopo$dY = dY
  dtopo$dZ = dZ
  
  dtopo
}

# compute beta in radians
# (the angle between vertical depth axis and 
#   the tangent vector of side of the triangular subfault)
# 
# ordering: x2-x1, x3-x2, x1-x3
# 
# requires self.corners to have been computed.
# JP's notes: 
# lon/lat: 
#     inputs: subfault$corners, subfault$lat
#     outputs: O1_list, O2_list 
get_leg_angles = function(subfault) {
  
  # TODO: put in a coordinate_specification == 'triangular' check here
  # x = numpy.array(subfault$corners)
  # y = numpy.zeros(x.shape)
  x = subfault$corners
  y = matrix(0, nrow=nrow(x), ncol=ncol(x))
  
  # convert to meters
  DEG2RAD = 2*pi/360
  LAT2METER = 111133.84012073894 #/10^3
  # y[:,0] = LAT2METER * cos( DEG2RAD*self.latitude )*x[:,0]
  # y[:,1] = LAT2METER * x[:,1]
  # y[:,2] = - abs(x[:,2])    # force sign
  y[,1] = LAT2METER * cos( DEG2RAD*subfault$lat)*x[,1]
  y[,2] = LAT2METER * x[,2]
  y[,3] = - abs(x[,3])    # force sign
  
  # v_list = [y[0,:] - y[1,:], y[1,:] - y[2,:], y[2,:] - y[0,:]]
  v_list = list(y[1,] - y[2,], y[2,] - y[3,], y[3,] - y[1,])
  
  # e3 = numpy.array([0.,0.,-1.])
  e3 = c(0,0,-1)
  
  # O1_list = []
  # O2_list = []
  # alpha_list = []
  # beta_list = []
  # reverse_list = [FALSE,FALSE,FALSE]
  O1_list = list()
  O2_list = list()
  alpha_list = list()
  beta_list = list()
  reverse_list = c(FALSE,FALSE,FALSE)
  
  j = 0
  for(v in v_list) {
    # vn = v/numpy.linalg.norm(v)
    vn = v/norm(v)
    k = j
    # l = (j+1)%3
    l = (j+1)%%3
    # if(vn[2] > 0) {
    if(vn[3] > 0) {
      vn = -vn    # point vn in depth direction
      # k = (j+1)%3
      k = (j+1)%%3
      l = j
      reverse_list[j+1] = TRUE
    }
    
    # O1_list.append(x[k,:].copy())  # set origin for the vector v
    # O2_list.append(x[l,:].copy())  # set dest.  for the vector v
    O1_list = c(O1_list, list(x[k+1,]))  # set origin for the vector v
    O2_list = c(O2_list, list(x[l+1,]))  # set dest.  for the vector v
    
    # alpha = numpy.arctan2(vn[0],vn[1])
    alpha = atan2(vn[1],vn[2])
    alpha_list = c(alpha_list, list(alpha))
    
    # beta = numpy.pi/2 \
    # - numpy.arctan(\
    #                numpy.divide(abs(vn[2]),
    #                             abs(numpy.sqrt(vn[0]**2 + vn[1]**2))))
    beta = pi/2 - 
      atan(abs(vn[3]) / abs(sqrt(vn[1]^2 + vn[2]^2)))
    beta_list = c(beta_list, list(beta))
    
    j = j+1
  }
  
  list(reverse_list=reverse_list,O1_list=O1_list,O2_list=O2_list, 
       alpha_list=alpha_list,beta_list=beta_list)
}

# compute angular dislocations at the *free surface* of the half space, 
# according to the paper
# 
# M. Comninou and J. Dundurs 
# Journal of Elasticity, Vol. 5, Nos.3-4, Nov 1975
# 
# The specific equations used in the papers are (1-29). 
# 
# :Input:
#   -  Y1, Y2, Y3
# Z1, Z2, Z3
# Yb1,Yb2,Yb3
# Zb1,Zb2,Zb3   : coordinates in meters
# 
# :Output:
#   - v11,v12,v13
# v21,v22,v23
# v31,v32,v33   : dislocation vectors
# 
# 
# For a more recent reference with comprehensive review see:
#   
#   Brendan J. Meade
# Computers & Geosciences, Vol. 33, Issue 8, pp 1064-1075
get_angular_dislocations_surface = function(self,Y1,Y2,Y3,beta,Odepth) {
  
  a = abs(Odepth)   #lazy
  
  nu = 0.25        # .5 * lambda / (lambda + mu) poisson ratio
  
  C = (2*pi)
  
  Z1 = cos(beta)*Y1 + a*sin(beta)
  Z3 = sin(beta)*Y1 - a*cos(beta)
  R = sqrt(Y1^2 + Y2^2 + a^2)
  
  F =  - atan2(Y2,Y1) + atan2(Y2*R*sin(beta),Y1*Z1 + Y2^2*cos(beta)) + atan2(Y2,Z1) 
  
  # Burgers vector (1,0,0)
  
  v11 = 1/C*((1 - (1 - 2*nu)/(tan(beta)^2))*F + Y2/(R + a)*
               ((1 - 2*nu)*(1/tan(beta) + Y1/(2*(R + a))) - Y1/R) - 
               (Y2/R)*((R*sin(beta) - Y1)*cos(beta)/(R - Z3)))
  
  v21 = 1/C*((1 - 2*nu)*((.5 + 1/(tan(beta)^2))*log(R + a) - 
                           1/tan(beta)/sin(beta)*log(R - Z3)) - 
               1/(R + a)*((1 - 2*nu)*(Y1/tan(beta) - .5*a - (Y2/(2*(R + a)))*Y2)) - 
               (Y2/R)*(Y2/(R+a)) + (Y2/R)*cos(beta)*(Y2/(R - Z3)))
  
  v31 = 1/C*((1 - 2*nu)*F/tan(beta) + Y2/(R + a)*(2*nu + a/R) - 
               (Y2/(R - Z3))*cos(beta)*(cos(beta) + a/R))
  
  # Burgers vector (0,1,0)
  
  v12 = 1/C*(-(1 - 2*nu)*((.5 - 1/(tan(beta)^2))*log(R + a) + 
                            cos(beta)/(tan(beta)^2)*log(R - Z3)) - 
               1/(R + a)*((1 - 2*nu)*(Y1/tan(beta) + .5*a + (Y1/(2*(R+a)))*Y1) - 
                            (Y1/R)*Y1) + (Z1/R)*((R*sin(beta) - Y1)/(R - Z3)))
  
  v22 = 1/C*((1 + (1 - 2*nu)/(tan(beta)^2))*F - Y2/(R + a)*
               ((1 - 2*nu)*(1/tan(beta) + Y1/(2*(R+a))) - Y1/R) - 
               (Y2/R)*(Z1/(R - Z3)))
  
  v32 = 1/C*(-(1 - 2*nu)/tan(beta)*(log(R + a) - cos(beta)*log(R - Z3)) - 
               Y1/(R + a)*(2*nu + a/R) + Z1/(R - Z3)*(cos(beta) + a/R))
  
  # Burgers vectors (0,0,1)
  
  v13 = 1/C*(Y2*(R*sin(beta) - Y1)*sin(beta)/(R*(R-Z3)))
  v23 = 1/C*(-Y2^2*sin(beta))/(R*(R - Z3))
  v33 = 1/C*(F + Y2*(R*cos(beta) + a)*sin(beta)/(R*(R - Z3)))
  
  
  list(v11=v11,v12=v12,v13=v13,
       v21=v21,v22=v22,v23=v23,
       v31=v31,v32=v32,v33=v33)
}

# compute a unit vector in the slip-direction (rake-direction)
# for a triangular fault
get_unit_slip_vector = function(subfault) {
  DEG2RAD = 2*pi/360
  # strike = numpy.deg2rad(self.strike)
  # dip = numpy.deg2rad(self.dip)
  # rake = numpy.deg2rad(self.rake)
  strike = DEG2RAD * subfault$strike
  dip = DEG2RAD * subfault$dip
  rake = DEG2RAD * subfault$rake
  
  # e1 = numpy.array([1.,0.,0.])
  # e2 = numpy.array([0.,1.,0.])
  # e3 = numpy.array([0.,0.,-1.])
  e1 = c(1,0,0)
  e2 = c(0,1,0)
  e3 = c(0,0,-1)
  
  u = sin(strike)*e1 + cos(strike)*e2
  v = cos(strike)*e1 - sin(strike)*e2
  
  w = sin(dip)*e3 + cos(dip)*v
  z = sin(-rake)*w + cos(-rake)*u
  
  z
}

# Calculate geometry for triangular subfaults
# 
# - Uses *corners* to calculate *centers*, *longitude*, *latitude*,
# *depth*, *strike*, *dip*, *length*, *width*.
# 
# - sets *coordinate_specification* as "triangular"
# 
# - Note that calculate_geometry() computes 
# long/lat/strike/dip/length/width to calculate centers/corners
# JP's notes:
# corners is 3x3 matrix where each row is (x, y, depth) of a single corner, 
# where x and y are in Easting and Northing in meters.
calculate_geometry_triangles = function(corners) {
  
  x0 = corners
  x = x0
  x[,3] = -abs(x[,3]) # set depth to be always negative
  
  if(FALSE) {
    # old coordinate transform
    
    # compute strike and dip direction
    # e3: vertical 
    # v1,v2: tangents from x0
    e3 = c(0,0,1)
    v1 = x[2,] - x[1,]
    v2 = x[3,] - x[1,]
    
    DEG2RAD = 2*pi/360
    v1[1] = v1[1] * LAT2METER * cos( DEG2RAD*x[1,2] ) 
    v2[1] = v2[1] * LAT2METER * cos( DEG2RAD*x[1,2] ) 
    v1[2] = v1[2] * LAT2METER 
    v2[2] = v2[2] * LAT2METER
  } 
  
  # x[:,0],x[:,1] = self._llz2utm(x[:,0],x[:,1],\
  #                               projection_zone=self._projection_zone)
  # x[,1:2] = projCSZ(x[,1:2], units="m")
  
  v1 = x[2,] - x[1,]
  v2 = x[3,] - x[1,]
  
  e3 = c(0,0,1)
  normal = pracma::cross(v1,v2)
  if(normal[3] < 0) {
    normal = -normal
    fix_orientation = TRUE
  } else {
    fix_orientation = FALSE
  }
  # not used for some reason:
  # strikev = cross(normal,e3)   # vector in strike direction
  
  a = normal[1]
  b = normal[2]
  c = normal[3]
  
  #Compute strike
  # strike_deg = rad2deg(arctan(-b/a))
  RAD2DEG = 360 / (2*pi)
  strike_deg = RAD2DEG * atan(-b/a)
  
  #Compute dip
  # beta = deg2rad(strike_deg + 90)
  # beta = 2*pi / 360 * ((strike_deg + 90) %% 360)
  beta = 2*pi / 360 * (strike_deg + 90)
  # m = numpy.array([sin(beta),cos(beta),0]) #Points in dip direction
  # n = numpy.array([a,b,c]) #Normal to the plane
  m = c(sin(beta),cos(beta),0) #Points in dip direction
  n = c(a,b,c) #Normal to the plane
  
  if(abs(c) < 1e-8) {
    dip_deg = 90   # vertical fault
  } else {
    # dip_deg = rad2deg(arcsin(m.dot(n)/(norm(m)*norm(n))))
    # browser() # dip and strike are not exactly correct. I think this is due to 
    # projection issues. Issue seems minor, though
    dip_deg = RAD2DEG * asin((m %*% n)/(norm(as.matrix(m))*norm(as.matrix(n))))
  }
  
  # dip should be between 0 and 90. If negative, reverse strike:
  if(dip_deg < 0) {
    strike_deg = strike_deg - 180.
    dip_deg = -dip_deg
  }
  if((0 > dip_deg) || (dip_deg > 90)) {
    stop(paste0("dip_deg = ", dip_deg, ", but must be between 0 and 90"))
  }
  
  #if(dip_deg > 30) {
   # browser()
  #}
  
  # keep strike_deg positive
  if(strike_deg < 0) {
    strike_deg = 360 + strike_deg
  }
  if(strike_deg < 0) {
    stop(paste0("strike_deg = ", strike_deg, ", but must be positive"))
  }
  
  strike = strike_deg
  dip = c(dip_deg)
  
  # find the center line
  xx = matrix(0, nrow=3, ncol=3)
  xx[1,] = (x0[2,] + x0[3,]) / 2. # midpt opposite a
  xx[2,] = (x0[1,] + x0[3,]) / 2. # midpt opposite b
  xx[3,] = (x0[1,] + x0[2,]) / 2. # midpt opposite c
  
  # centers seems to be [longitude, easting (m)] of furthest south point?
  # Doesn't make sense, but centers doesn't seem to be used anywhere so it's 
  # fine
  # i = numpy.argmin(xx[2,:])
  i = which.min(xx[3,])
  
  # centers = [x[,i].tolist(), xx[,i].tolist()]
  centers = cbind(x[,i], xx[,i])
  
  if(x[3,i] <= xx[3,i]) {
    # self._centers.reverse()
    centers = rev(centers)
  }
  
  # xcenter = numpy.mean(xx, axis=0)
  # longitude = xcenter[0]
  # latitude = xcenter[1]
  # depth = xcenter[2]
  xcenter = colMeans(xx)
  longitude = xcenter[1]
  latitude = xcenter[2]
  depth = xcenter[3]
  
  # length and width are set to sqrt(area): 
  # this is set temporarily so that Fault.Mw() can be computed
  
  area = norm(normal) / 2.
  length = sqrt(area)
  width = sqrt(area)
  
  # a list summarizing the subfault
  list(corners=corners, centers=centers, 
       lon=longitude, lat=latitude, 
       depth=depth, strike=strike, dip=dip, 
       fix_orientation=fix_orientation)
}

# compute coordinates
# 
# :Input:
#   - X1,X2,X3: longitude,latitude,depth
# - alpha: angle of the vertical hyperplane
# (measured from north-direction, =strike)
# - beta: angle of the angular dislocation 
# (angle between two inf lines, 
#   measured from depth-direction)
# - Olong,Olat,Odepth: longitude,latitude,depth of 
# origin of y1-y2-y3 coordinates
# 
# :Output:
#   
#   - tuple (y,z,ybar,zbar)
# - each numpy array of same shape as x
# containing *angular* coordinates as well as its mirrored image
get_halfspace_coords = function(subfault,X1,X2,X3,alpha,beta,Olong,Olat,Odepth) {
  
  dims = dim(X1)
  Odepth = abs(Odepth)
  
  # convert lat/long to meters
  DEG2RAD = 2*pi/360
  LAT2METER = 111133.84012073894 #/10^3
  X1 = LAT2METER * cos( DEG2RAD*subfault$lat ) * (X1 - Olong)
  X2 = LAT2METER * (X2 - Olat)
  
  # Y1 = numpy.zeros(dims)       # yi-coordinates
  # Y2 = numpy.zeros(dims)       # yi-coordinates
  # Y3 = numpy.zeros(dims)       # yi-coordinates
  # 
  # Z1 = numpy.zeros(dims)       # yi coordinates rot. by beta
  # Z2 = numpy.zeros(dims)       # yi coordinates rot. by beta
  # Z3 = numpy.zeros(dims)       # yi coordinates rot. by beta
  # 
  # Yb1 = numpy.zeros(dims)    # mirrored yi-coordinates
  # Yb2 = numpy.zeros(dims)    # mirrored yi-coordinates
  # Yb3 = numpy.zeros(dims)    # mirrored yi-coordinates
  # 
  # Zb1 = numpy.zeros(dims)    # mirrored yi-coordinates rot. by beta
  # Zb2 = numpy.zeros(dims)    # mirrored yi-coordinates rot. by beta
  # Zb3 = numpy.zeros(dims)    # mirrored yi-coordinates rot. by beta
  Y1 =  matrix(0, nrow=dims[1], ncol=dims[2])       # yi-coordinates
  Y2 =  matrix(0, nrow=dims[1], ncol=dims[2])       # yi-coordinates
  Y3 =  matrix(0, nrow=dims[1], ncol=dims[2])       # yi-coordinates
  
  Z1 =  matrix(0, nrow=dims[1], ncol=dims[2])       # yi coordinates rot. by beta
  Z2 =  matrix(0, nrow=dims[1], ncol=dims[2])       # yi coordinates rot. by beta
  Z3 =  matrix(0, nrow=dims[1], ncol=dims[2])       # yi coordinates rot. by beta
  
  Yb1 = matrix(0, nrow=dims[1], ncol=dims[2])    # mirrored yi-coordinates
  Yb2 = matrix(0, nrow=dims[1], ncol=dims[2])    # mirrored yi-coordinates
  Yb3 = matrix(0, nrow=dims[1], ncol=dims[2])    # mirrored yi-coordinates
  
  Zb1 = matrix(0, nrow=dims[1], ncol=dims[2])    # mirrored yi-coordinates rot. by beta
  Zb2 = matrix(0, nrow=dims[1], ncol=dims[2])    # mirrored yi-coordinates rot. by beta
  Zb3 = matrix(0, nrow=dims[1], ncol=dims[2])    # mirrored yi-coordinates rot. by beta
  
  # rotate by -alpha in long/lat plane
  Y1 = sin(alpha)*X1 + cos(alpha)*X2
  Y2 = cos(alpha)*X1 - sin(alpha)*X2
  Y3 = X3 - abs(Odepth)
  
  Z1 = cos(beta)*Y1 - sin(beta)*Y3
  Z2 = Y2
  Z3 = sin(beta)*Y1 + cos(beta)*Y3
  
  Yb1 = Y1
  Yb2 = Y2
  Yb3 = X3 + abs(Odepth)
  
  Zb1 =  cos(beta)*Y1 + sin(beta)*Yb3
  Zb2 =  Y2
  Zb3 = -sin(beta)*Y1 + cos(beta)*Yb3
  
  list(Y1=Y1,Y2=Y2,Y3=Y3,
       Z1=Z1,Z2=Z2,Z3=Z3,
       Yb1=Yb1,Yb2=Yb2,Yb3=Yb3,
       Zb1=Zb1,Zb2=Zb2,Zb3=Zb3)
}

# compute coordinate transforms by computing 
# 
#  [sin  cos   0] [v11 v12 v13] [sin  cos   0]
#  |cos -sin   0| |v21 v22 v23| |cos -sin   0|
#  [  0    0   1] [v31 v32 v33] [  0    0   1]
# 
coord_transform = function(subfault,v11,v12,v13,v21,v22,v23,v31,v32,v33,alpha) {
  
  w11 = sin(alpha)*v11 + cos(alpha)*v12
  w12 = cos(alpha)*v11 - sin(alpha)*v12
  w13 = v13
  
  w21 = sin(alpha)*v21 + cos(alpha)*v22
  w22 = cos(alpha)*v21 - sin(alpha)*v22
  w23 = v23
  
  w31 = sin(alpha)*v31 + cos(alpha)*v32
  w32 = cos(alpha)*v31 - sin(alpha)*v32
  w33 = v33
  
  v11 = sin(alpha)*w11 + cos(alpha)*w21
  v12 = sin(alpha)*w12 + cos(alpha)*w22
  v13 = sin(alpha)*w13 + cos(alpha)*w23
  
  v21 = cos(alpha)*w11 - sin(alpha)*w21
  v22 = cos(alpha)*w12 - sin(alpha)*w22
  v23 = cos(alpha)*w13 - sin(alpha)*w23
  
  v31 = w31
  v32 = w32
  v33 = w33
  
  list(v11=v11,v12=v12,v13=v13,
       v21=v21,v22=v22,v23=v23,
       v31=v31,v32=v32,v33=v33)
}

# _llz2utm = function(self,lon,lat,projection_zone=NULL) {
#   # Convert lat,lon to UTM
#   # 
#   # originally written by Diego Melgar (Univ of Oregon)
#   
#   
#   array_dims = lon.shape
# 
#   lon = lon.flatten()
#   lat = lat.flatten()
# 
#   x=zeros(lon.shape)
#   y=zeros(lon.shape)
#   zone=zeros(lon.shape)
#   
#   #b=chararray(lon.shape) # gives byte error
#   # b=len(lon)*['A']  # list of characters, modified in loop below
#   b=rep('A', length(lon))  # list of characters, modified in loop below
#   if(is.null(projection_zone)) {
#     #Determine most suitable UTM zone
#     for(k in range(len(lon))) {
#       x,y,zone[k],b[k]=utm.from_latlon(lat[k],lon[k])
#     }
#       zone_mode=mode(zone)
#       i=where(zone==zone_mode)[0]
#       letter=b[i[0]]
#       z=str(int(zone[0]))+letter
#   } else {
#     z=projection_zone
#   }
#   p0 = Proj(proj='utm',zone=z,ellps='WGS84')
#   x,y=p0(lon,lat)
#   
#   x = x.reshape(array_dims)
#   y = y.reshape(array_dims)
#   list(x,y)
# }