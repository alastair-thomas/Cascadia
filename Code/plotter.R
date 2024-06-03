plotCSZOutline = function(maxDepth = 30){
  # first load in the Slab 2.0 geometry
  slab = loadSlab2()
  lonLat = cbind(slab$lon - 360, slab$lat) # maybe need a -360 here
  depths = slab$depth
  
  # get projected coordinates
  xy = projCSZ(lonLat, units="m")
  
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
  # otherwise some depths become NA when depth is interpolated from Slab2
  v = terra::vect(concaveHull, "polygons") # convert to format needed for buffer function
  v = terra::buffer(v, width=-0.1) # shrink the hull by 0.1km
  v = terra::as.points(v) # convert to points
  concaveHull = terra::crds(v) # get the hull in same format as before
  
  g = plotBase(proj="UTM")
  
  plotDF = data.frame(x = concaveHull[,1],
                      y = concaveHull[,2])
  
  g = g +
    geom_polygon(data=plotDF,
                 aes(x=x, y=y),
                 colour="red", fill="gray", linewidth=1, alpha=0.75)
  
  # sort the plotting limits
  limits = data.frame(x = -c(128, 119.5),
                      y = c(40, 50))
  limits = st_as_sf(x=limits,
                    coords = c("x", "y"),
                    crs=st_crs("EPSG:4326"))
  limits = toUTM(limits)
  g = g +
    coord_sf(xlim = c(st_coordinates(limits)[,1]),
             ylim = c(st_coordinates(limits)[,2]),
             crs=st_crs("EPSG:32610"))
  
  plot(g)
}
plotTaper = function(lambdaDraws,
                     taperPrior=c(93.757027137, 0.001316695)){
  depths = seq(0,30,0.01)
  lambdaQ = stats::quantile(lambdaDraws, probs=c(0.025, 0.5, 0.975))
  tlower = exp(-lambdaQ[1]*depths)
  tmedian = exp(-lambdaQ[2]*depths)
  tupper = exp(-lambdaQ[3]*depths)
  
  lambdaPriorLower = qgamma(0.025, shape=taperPrior[1], scale=taperPrior[2])
  lambdaPriorUpper = qgamma(0.975, shape=taperPrior[1], scale=taperPrior[2])
  tPriorLower = exp(-lambdaPriorLower*depths)
  tPriorUpper = exp(-lambdaPriorUpper*depths)
  
  plotDF = data.frame(depth = depths,
                      lower = tlower,
                      taper = tmedian,
                      upper = tupper,
                      Plower = tPriorLower,
                      Pupper = tPriorUpper)
  ggplot(plotDF) +
    geom_ribbon(aes(x=depth, ymin=lower, ymax=upper), fill="lightblue", alpha=0.75) +
    geom_line(aes(x=depth, y=taper), colour="black") +
    geom_line(aes(x=depth, y=Plower), colour="orange", linetype="dashed") +
    geom_line(aes(x=depth, y=Pupper), colour="orange", linetype="dashed") +
    labs(x="Depth (km)", y="Taper")
}

plotInnerOuterError = function(fullScoreDF, validScoreDF){
  # Plot the inside vs outside CRPS scores
  scores = rbind(fullScoreDF, validScoreDF)
  
  scores$model = as.factor(c(rep("Full Model", dim(fullScoreDF)[1]),
                             rep("Validation", dim(validScoreDF)[1])))
  
  scores$Event = as.factor(scores$Event)
  
  ggplot(scores, aes(x=Event,y=CRPS,fill=model)) +
    geom_boxplot(outlier.shape = NA) +
    stat_summary(fun=mean,
                 geom='point',
                 shape=20,
                 size=3,
                 aes(colour=model),
                 position = position_dodge(width = 0.75)) +
    scale_fill_scico_d(palette="berlin") +
    scale_colour_scico_d(palette="roma", direction = -1) +
    labs(fill="Model", colour="Model Mean")
}

plotCorrelation = function(psiDraws, anisoPrior=0.38){
  
  lowerPsi = quantile(psiDraws, 0.025)
  medianPsi = quantile(psiDraws, 0.5)
  upperPsi = quantile(psiDraws, 0.975)
  
  # Generate points for a standard ellipse centered at (0, 0)
  t = seq(0, 2*pi, length.out = 1000)
  baseY = cos(t)
  baseX = sin(t)
  
  # Rotation
  # measured anticlockwise from the y
  thetaMajor = deg2rad(-12.39776)
  rotMat = matrix(c(cos(thetaMajor), -sin(thetaMajor),
                    sin(thetaMajor), cos(thetaMajor)),
                  nrow = 2)
  
  
  # Apply the rotation matrix to the points
  # Median
  medianEllipsePoints = rotMat %*% rbind((1/medianPsi) * baseX,
                                         medianPsi * baseY)
  medianDF = data.frame(x = medianEllipsePoints[1,],
                        y = medianEllipsePoints[2,])
  
  # Lower
  lowerEllipsePoints = rotMat %*% rbind((1/lowerPsi) * baseX,
                                         lowerPsi * baseY)
  lowerDF = data.frame(x = lowerEllipsePoints[1,],
                        y = lowerEllipsePoints[2,])
  
  # Upper
  upperEllipsePoints = rotMat %*% rbind((1/upperPsi) * baseX,
                                         upperPsi * baseY)
  upperDF = data.frame(x = upperEllipsePoints[1,],
                        y = upperEllipsePoints[2,])
  
  # Now sort the priors
  psiPriorL = exp(qnorm(0.025, mean=0, sd=anisoPrior))
  psiPriorU = exp(qnorm(0.975, mean=0, sd=anisoPrior))
  # Lower Prior
  lowerAnisoEllipsePoints = rotMat %*% rbind((1/psiPriorL) * baseX,
                                             psiPriorL * baseY)
  lowerAnisoDF = data.frame(x = lowerAnisoEllipsePoints[1,],
                            y = lowerAnisoEllipsePoints[2,])
  # Upper Prior
  upperAnisoEllipsePoints = rotMat %*% rbind((1/psiPriorU) * baseX,
                                             psiPriorU * baseY)
  upperAnisoDF = data.frame(x = upperAnisoEllipsePoints[1,],
                            y = upperAnisoEllipsePoints[2,])
    
  # Plot the ellipse
  ggplot() +
    geom_polygon(data=lowerDF, aes(x=x, y=y), fill = "#88CCEE", alpha=0.75) +
    geom_polygon(data=upperDF, aes(x=x, y=y), fill = "#88CCEE", alpha=0.75) +
    geom_path(data=medianDF, aes(x=x, y=y), linewidth=1, colour="black") +
    geom_path(data=lowerAnisoDF, aes(x=x, y=y), linewidth=1, linetype="dashed", colour="orange") +
    geom_path(data=upperAnisoDF, aes(x=x, y=y), linewidth=1, linetype="dashed", colour="orange") +
    coord_fixed() +
    labs(x = "X",
         y = "Y")
}

plotFixedParameters = function(DF, histLabel){
  
  plotDF = data.frame()
  for (p in 1:dim(DF)[2]){
    lower = quantile(DF[,p], 0.025)
    upper = quantile(DF[,p], 0.975)
    mu    = mean(DF[,p])
    thisDF = data.frame(value = DF[,p],
                        var   = rep(paste(histLabel[p],
                                         "\nMean=", signif(mean(DF[,p]), 3),
                                         ", Median=", signif(median(DF[,p]), 3)), dim(DF)[1]),
                        lower = rep(lower, dim(DF)[1]),
                        upper = rep(upper, dim(DF)[1]),
                        mu    = mu)
    
    plotDF = rbind(plotDF, thisDF)
  }
  
  # Create a grid of histograms
  # Colours chosen from:
  # https://stackoverflow.com/questions/9563711/r-color-palettes-for-many-data-classes
  
  g = ggplot(plotDF) +
    geom_histogram(aes(x = value), bins=25, alpha = 0.75, fill="#88CCEE") +
    geom_vline(aes(xintercept=lower), linetype="dashed", linewidth=1, colour="#CC6677") +
    geom_vline(aes(xintercept=upper), linetype="dashed", linewidth=1, colour="#CC6677") +
    geom_vline(aes(xintercept=mu), linetype="dashed", linewidth=1, colour="#DDCC77") +
    facet_wrap(~var, scales="free") +
    labs(x    = "Parameter Value",
         y    = "Count")
  plot(g)
}

# data - a dataframe with:
#       - Latitude
#       - Subsidence
#       - Uncertainty
plotSubsidenceLat = function(DR){
  sites = c("Port Alberni",
            "Tofino area",
            "Quinault River",
            "Copalis River Estuary",
            "Grays Harbour",
            "Willapa Bay",
            "Columbia River",
            "Necanicum River Estuary",
            "Nehalem River",
            "Tillamook Bay",
            "Netarts Bay",
            "Nestucca Bay",
            "Salmon River Estuary",
            "Siletz Bay",
            "Yaquina Bay",
            "Alsea Bay",
            "Siuslaw River",
            "Umpqua River",
            "Coos Bay",
            "Coquille River",
            "Sixes River",
            "Humboldt Bay",
            "Eel River")
  DR$Site = factor(DR$Site,
                   levels=rev(sites))
  
  Sites = aggregate(Lat~Site, data=DR4, FUN=mean)
  yBreaks = Sites$Lat
  yLabels = Sites$Site
  
  g = ggplot()+
       geom_point(data=DR,
                  aes(x=subsidence, y=Lat, colour=Uncertainty),
                  size=1)+
    scale_colour_viridis(option="magma") +
    scale_y_continuous(breaks = yBreaks, labels = yLabels) +
    theme(panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          legend.key.height = unit(2, "cm")) +
    labs(x = "Subsidence (m)",
         y = "Latitude (°)",
         colour="Uncertainty (m)")
  return(g)
}

plotSites = function(data){
  # gets the base map of the CSZ
  g = plotBase(proj="UTM", scale=1.4)
  
  # create data to plot
  Sites = data.frame(Site=data$Site,
                     Lat=data$Lat,
                     Lon=data$Lon,
                     Sub=data$subsidence)
  Sites = aggregate(.~Site, data=Sites, mean)
  Sites$count = as.factor(data.frame(table(data$Site))$Freq)
  
  # convert coordinates to UTM
  xy = cbind(Sites$Lon, Sites$Lat)
  xy = projCSZ(xy, units="m")
  Sites$Lon = xy[,1]
  Sites$Lat = xy[,2]
  
  # plot locations of sites and colour by number of estimates
  g = g +
    geom_point(data=Sites,
               aes(x=Lon, y=Lat, colour=count),
               size=2.5) +
    scale_colour_viridis_d(option="magma", direction=-1)
  
  # add the site name as a label
  xlims = data.frame(x=c(-130, -127),
                     y=c(40, 50))
  xlims = st_as_sf(x=xlims,
                  coords = c("x", "y"),
                  crs=st_crs("EPSG:4326"))
  xlims = toUTM(xlims)
  xlims = c(st_coordinates(xlims)[,1])
  g = g +
    ggrepel::geom_text_repel(data=Sites,
                             aes(x=Lon, y=Lat, label=Site),
                             size=3, hjust = 0,
                             xlim = xlims, force = 5,
                             direction="y", box.padding=0.15)
  
  # sort the plotting limits
  limits = data.frame(x = -c(130, 120),
                      y = c(40, 50))
  limits = st_as_sf(x=limits,
                    coords = c("x", "y"),
                    crs=st_crs("EPSG:4326"))
  limits = toUTM(limits)
  g = g +
    coord_sf(xlim = c(st_coordinates(limits)[,1]),
             ylim = c(st_coordinates(limits)[,2]),
             crs=st_crs("EPSG:32610"))
  
  # change aesthetics
  g = g +
    labs(colour='Number of\nSubsidence Estimates') +
    theme(legend.position = "right",
          legend.key.height = unit(0.5, 'cm'))
  
  return(g)
}

# a function to plot the number of subsidence estimates for each earthquake
plotEventCounts = function(data){
  df = data.frame(data) # makes sure using a dataframe
  df = df[!grepl("a|b|R", df$event), ] # gets only full margin ruptures
  
  eventsOrder = c("T1", "T2", "T3", "T4", "T5", "T6", "T7", "T8", "T9", "T10", "T11", "T12")
  df$event = factor(df$event, levels=eventsOrder)
  
  ggplot(df, aes(x = event)) +
    geom_bar() +
    labs(title = "",
         x = "Seismic Event",
         y = "Number of Subsidence Estimates")
}

plotEventSites = function(data){
  # Calculate counts of each event factor from the original dataframe
  eventCounts = data %>%
                count(event) %>%
                rename(event = event, eventCount = n)
  
  # Summarize data by Site and event
  plotDF = data %>%
           group_by(Site, event) %>%
           summarise(siteLon = mean(Lon),
           siteLat = mean(Lat)) %>%
           as.data.frame()
  
  # Join event counts with the summarized dataframe
  plotDF = left_join(plotDF, eventCounts, by = "event")
  
  # make event a factor
  eventCount = c(197,24,26,42,52,30,27,11,8,7,3,2)
  plotDF$event = factor(plotDF$event,
                        levels=paste("T", 1:12, sep = ""))
  levels(plotDF$event) = paste("T", 1:12, "\n", eventCount, " Estimates", sep = "")
  
  plotDF$eventCount = as.factor(plotDF$eventCount)
  
  
  # plot base without labels
  g = plotBase(proj="UTM", labels=F, countryBoundary=F)
  
  # convert coordinates to UTM
  xy = cbind(plotDF$siteLon, plotDF$siteLat)
  xy = projCSZ(xy, units="m")
  plotDF$siteLon = xy[,1]
  plotDF$siteLat = xy[,2]
  
  # plot
  g = g +
    geom_point(data=plotDF,
               aes(x=siteLon, y=siteLat),
               colour="red") +
    facet_wrap(~event,
               nrow=2)
    
  
  # sort the x axis labels
  g = g +
    scale_x_continuous(breaks=-c(124,126))
  
  
  # sort the legend
  g = g +
    labs(colour="Number of\nSubsidence Estimates")
  
  
  # sort the plotting limits
  limits = data.frame(x = -c(127, 122),
                      y = c(40, 50))
  limits = st_as_sf(x=limits,
                    coords = c("x", "y"),
                    crs=st_crs("EPSG:4326"))
  limits = toUTM(limits)
  g = g +
    coord_sf(xlim = c(st_coordinates(limits)[,1]),
             ylim = c(st_coordinates(limits)[,2]),
             crs=st_crs("EPSG:32610"))
  return(g)
}

# plot the number of observations per site on the map
# size of spot is number of observations
# only for the full margin ruptures.

plotNumObersvations = function(scale=2){
  # gets the base map of the CSZ
  g = plotBase(scale=scale)
  
  Sites = data.frame(Site=dr1$Site, Lat=dr1$Lat, Lon=dr1$Lon)
  Sites = aggregate(.~Site, data=Sites, mean)
  Sites$count = as.numeric(data.frame(table(dr1$Site))$Freq)
  
  g = g +
    geom_point(data=Sites, aes(x=Lon, y=Lat, size=count, alpha=count, colour=count)) +
    scale_size_continuous(range = c(1, scale*8)) +
    scale_alpha_continuous(range = c(0.5, 1), trans="reverse") +
    scale_colour_gradient2(low="red", mid="yellow", high="green", midpoint=50) +
    coord_sf(xlim=-c(131, 120), ylim=c(40, 50)) + 
    guides(
      size = guide_legend(title = "Number of\nSubsidence Estimates"),
      alpha = guide_legend(title = "Number of\nSubsidence Estimates"),
      colour = guide_legend(title = "Number of\nSubsidence Estimates")
    )
  
  return(g)
}

# plots the Slab2 grids
plotGrid = function(scale=2, varName = "depth", limitDepth=T, proj="UTM", legendTitle="Depth (km)"){
  
  # read in the slab2 data for depth
  grid = readGrid(varName)
  
  # make strike centered around zero
  if (varName == "strike"){
    grid$Z[grid$Z > 100] = grid$Z[grid$Z > 100] - 360
  }
  
  # limits the region of depths up to 30km if wanted
  if (limitDepth){
    # limit to 30km depths
    if (varName == "depth"){
      depths = -grid$Z
      grid$Z = -grid$Z
    } else{
      tempGrid = readGrid("depth")
      depths = -tempGrid$Z
    }
    depthMask = which(depths <= 30)
    grid2 = grid[depthMask,]
    
    # create the limits
    limits = data.frame(x = -c(128.5, 122),
                        y = c(39.8, 50.2))
    
  } else{
    grid2 = grid
    grid2$Z = -grid2$Z
    
    # create the limits
    limits = data.frame(x = -c(128.5, 118),
                        y = c(38, 50))
  }
  
  # set the map in the correct projection
  # works because everything is an sf object
  if (proj == "UTM"){
    # gets the base map in given projection
    g = plotBase(scale=scale, proj="UTM", labels=F, countryBoundary = F)
    
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
    
    limits = st_as_sf(x=limits,
                      coords = c("x", "y"),
                      crs=st_crs("EPSG:4326"))
    limits = toUTM(limits)
    
    # plot
    g = g +
      geom_raster(data=grid4, aes(x=x, y=y, fill=z)) +
      scale_fill_viridis(alpha=0.75, option="mako") +
      labs(fill=legendTitle) +
      theme(legend.key.height = unit(2, "cm")) +
      coord_sf(xlim = c(st_coordinates(limits)[,1]),
               ylim = c(st_coordinates(limits)[,2]),
               crs=st_crs("EPSG:32610"))
  }
  else{
    # gets the base map in given projection
    g = plotBase(scale=scale, labels=FALSE)
    
    g = g +
      geom_raster(data=grid2, aes(x=X, y=Y, fill=Z)) +
      scale_fill_viridis(alpha=0.75, option="mako") +
      labs(fill=legendTitle) +
      theme(legend.key.height = unit(2, "cm")) + 
      coord_sf(xlim = -c(128.5, 122), ylim = c(39.8, 50.2), crs=st_crs("EPSG:4326"))
  }
  
  return(g)
}

plotErrors = function(error, signError, scale=2){
  # gets the base map of the CSZ
  g = plotBase(scale=scale)
  
  Sites = data.frame(Site=DR3$Site[DR3$event == earthquake],
                     Lat=DR3$Lat[DR3$event == earthquake],
                     Lon=DR3$Lon[DR3$event == earthquake])
  Sites$Error = error
  Sites$sign = signError
  
  
  g = g +
    geom_point(data=Sites, aes(x=Lon, y=Lat, colour=Error, shape=sign), size=scale*1.5) +
    scale_colour_gradient(low="red", high="green", name="Absolute Error (m)", trans = "reverse") +
    scale_shape_manual(name = "Same Sign?", values = c(15, 20), labels = c("No", "Yes")) + 
    coord_sf(xlim=-c(128, 122), ylim=c(40, 50))
  
  return(g)
}

plotOneSubsidencePrediction = function(subsidence, DR,
                                       event="T1"){
  
  Nj = dim(subsidence)[1]
  
  sites = c("Port Alberni",
            "Tofino area",
            "Quinault River",
            "Copalis River Estuary",
            "Grays Harbour",
            "Willapa Bay",
            "Columbia River",
            "Necanicum River Estuary",
            "Nehalem River",
            "Tillamook Bay",
            "Netarts Bay",
            "Nestucca Bay",
            "Salmon River Estuary",
            "Siletz Bay",
            "Yaquina Bay",
            "Alsea Bay",
            "Siuslaw River",
            "Umpqua River",
            "Coos Bay",
            "Coquille River",
            "Sixes River",
            "Humboldt Bay",
            "Eel River")
  DR$Site = factor(DR$Site,
                   levels=rev(sites))
  
  Sites = aggregate(Lat~Site, data=DR, FUN=mean)
  yBreaks = Sites$Lat
  yLabels = Sites$Site
  
  subReal = DR$subsidence[DR$event == event]
  subMean = apply(subsidence, 1, mean)
  subSD   = apply(subsidence, 1, sd)
  
  plotDF = data.frame(subsidence = c(subReal, subMean),
                      lat        = c(rep(DR$Lat[DR$event == event], 2)),
                      sd         = c(DR$Uncertainty[DR$event == event],
                                     subSD),
                      realPred   = c(rep("Estimate", Nj),
                                     rep("Prediction", Nj)))
  g = ggplot(plotDF) +
    geom_point(aes(x=subsidence,
                   y=lat,
                   colour=realPred)) +
    labs(x      = "Subsidence (m)",
         y      = "Latitude (°)",
         colour = "Data Type") +
    scale_y_continuous(breaks = yBreaks, labels = yLabels)
  
  
  error = abs(subReal - subMean)
  sign  = as.factor(ifelse(subMean*subMean >= 0, "Yes", "No"))
  
  plotDF2 = data.frame(subsidence = subMean,
                       lat        = DR$Lat[DR$event == event],
                       error      = error,
                       sign       = sign)
  g2 = ggplot(plotDF2) +
    geom_point(aes(x=subsidence, y=lat,
                   colour=error, shape=sign)) +
    scale_colour_viridis(option="plasma") +
    labs(x = "Subsidence (m)",
         y = "Latitude (°)",
         colour="Error (m)",
         shape="Same Sign?") +
    scale_y_continuous(breaks = yBreaks, labels = yLabels)
    
  plot(g)
  plot(g2)
}

plotAllSubsidencePredicition = function(subsidences, DR){
  E = length(subsidences)
  
  plotDF = data.frame()
  plotDF2 = data.frame()
  earthquakes = c("T1","T2","T3","T4","T5","T6","T7","T8","T9","T10","T11","T12")
  for (e in 1:E){
    Nj = dim(subsidences[[e]])[1]
    
    # mean subsidence
    meanSub = rowMeans(subsidences[[e]])
    
    thisDF = data.frame(subsidence = c(DR$subsidence[DR$event == earthquakes[e]],
                                       meanSub),
                        lat        = rep(DR$Lat[DR$event == earthquakes[e]], 2),
                        quake      = rep(earthquakes[e], 2*Nj),
                        realPred   = c(rep("Estimate", Nj), rep("Prediction", Nj)))
    
    thisDF2 = data.frame(subsidence = meanSub,
                         lat        = DR$Lat[DR$event == earthquakes[e]],
                         error      = abs(DR$subsidence[DR$event == earthquakes[e]] - meanSub),
                         residual   = DR$subsidence[DR$event == earthquakes[e]] - meanSub,
                         sign       = as.factor(ifelse(DR$subsidence[DR$event == earthquakes[e]]*meanSub >= 0, "Yes", "No")),
                         quake      = rep(earthquakes[e], Nj))
    
    plotDF = rbind(plotDF, thisDF)
    plotDF2 = rbind(plotDF2, thisDF2)
  }
  
  plotDF$quake = factor(plotDF$quake,
                        levels = earthquakes,
                        ordered = T)
  
  plotDF2$quake = factor(plotDF2$quake,
                         levels = earthquakes,
                         ordered = T)
  
  plotDF$realPred = factor(plotDF$realPred,
                           levels = c("Prediction", "Estimate"),
                           ordered = T)
  
  g = ggplot(plotDF) +
    geom_point(aes(x=subsidence, y=lat, colour=realPred, shape=realPred)) +
    scale_colour_manual(values = c("Prediction"="#661100", "Estimate"="#117733")) +
    scale_shape_manual(values = c("Prediction"=17, "Estimate"=1)) +
    facet_wrap(~quake, nrow=2) +
    labs(x = "Subsidence (m)",
         y = "Latitude (°)") +
    theme(legend.position = "none")
  
  ga = ggplot(plotDF) +
    geom_point(aes(x=subsidence, y=lat, colour=realPred, shape=realPred)) +
    scale_shape_manual(values = c("Prediction"=17, "Estimate"=1)) +
    scale_colour_manual(values = c("Prediction"="#661100", "Estimate"="#117733")) +
    labs(x = "Subsidence (m)",
         y = "Latitude (°)") +
    theme(legend.position = "none")
  
  # g2 = ggplot(plotDF2) +
  #   geom_point(aes(x=subsidence, y=lat, colour=error, shape=sign)) +
  #   scale_colour_viridis(option="plasma") +
  #   facet_wrap(~quake, nrow=2) +
  #   labs(x = "Subsidence (m)",
  #        y = "Latitude (°)",
  #        colour = "Absolute Error (m)",
  #        shape = "Same Sign?")
  # 
  # g2a = ggplot(plotDF2) +
  #   geom_point(aes(x=subsidence, y=lat, colour=error, shape=sign)) +
  #   scale_colour_viridis(option="plasma") +
  #   labs(x = "Subsidence (m)",
  #        y = "Latitude (°)",
  #        colour = "Absolute Error (m)",
  #        shape = "Same Sign?")
  
  g3 = ggplot(plotDF2) +
    geom_point(aes(x=residual, y=lat)) +
    labs(x = "Residual (m)",
         y = "Latitude (°)")
    
  g3a = ggplot(plotDF2) +
      geom_point(aes(x=residual, y=lat)) +
      facet_wrap(~quake, nrow=2) +
    labs(x = "Residual (m)",
         y = "Latitude (°)")
  
  plot(g)
  plot(ga)
  # plot(g2)
  # plot(g2a)
  plot(g3)
  plot(g3a)
}

plotSubsidenceGrid = function(subsidences, lonLat, DR){
  # gets the base map of the CSZ
  g = plotBase(proj="UTM", labels=F, countryBoundary=F)
  g2 = plotBase(proj="UTM", labels=F, countryBoundary=F)
  xy = projCSZ(lonLat, unit="m")
  
  meanSub = apply(subsidences, 1, mean)
  sdSub   = apply(subsidences, 1, sd)
  
  # Interpolate onto a regular grid
  xs = seq(min(xy[,1]), max(xy[,1]), length.out = 200)
  ys = seq(min(xy[,2]), max(xy[,2]), length.out = 400)
  meanMesh = akima::interp(xy[,1], xy[,2], meanSub, xo=xs, yo=ys)
  meanGrid = akima::interp2xyz(meanMesh, data.frame=T)
  meanGrid = na.omit(meanGrid)
  
  sdMesh = akima::interp(xy[,1], xy[,2], sdSub, xo=xs, yo=ys)
  sdGrid = akima::interp2xyz(sdMesh, data.frame=T)
  sdGrid = na.omit(sdGrid)
  
  g = g +
    geom_raster(data=meanGrid, aes(x=x, y=y, fill=z), alpha=0.75) +
    scale_fill_scico(midpoint=0, palette="vik", direction=1) +
    labs(fill="Subsidence\nMean (m)") +
    theme(legend.key.height = unit(2, 'cm'))
  
  g2 = g2 +
    geom_raster(data=sdGrid, aes(x=x, y=y, fill=z), alpha=0.75) +
    scale_fill_scico(palette="lajolla", direction=-1) +
    labs(fill="Subsidence\nSD (m)") +
    theme(legend.key.height = unit(2, 'cm'))
  
  # Add data points
  # add the points where subsidence estimates are
  quake="T1"
  xy = cbind(DR$Lon, DR$Lat)
  xy = projCSZ(xy, units="m")
  plotDF2 = data.frame(x     = xy[DR$event == quake, 1],
                       y     = xy[DR$event == quake, 2],
                       site  = DR$Site[DR$event == quake])
  
  plotDF2 = plotDF2 %>%
    group_by(site) %>%
    summarize(
      x     = mean(x, na.rm = TRUE),
      y     = mean(y, na.rm = TRUE),
      count = n()
    )
  
  g = g +
    geom_point(data=plotDF2,
               aes(x=x, y=y), colour="red")
  g2 = g2 +
    geom_point(data=plotDF2,
               aes(x=x, y=y), colour="red")
  
  # Fix the limits
  limits = data.frame(x = -c(128, 123),
                      y = c(40, 50))
  limits = st_as_sf(x=limits,
                    coords = c("x", "y"),
                    crs=st_crs("EPSG:4326"))
  limits = toUTM(limits)
  g = g +
    coord_sf(xlim = c(st_coordinates(limits)[,1]),
             ylim = c(st_coordinates(limits)[,2]),
             crs=st_crs("EPSG:32610"))
  
  g2 = g2 +
    coord_sf(xlim = c(st_coordinates(limits)[,1]),
             ylim = c(st_coordinates(limits)[,2]),
             crs=st_crs("EPSG:32610"))
  
  plot(g)
  plot(g2)
}


plotFutureSubsidenceGrid = function(subsidences, lonLat){
  # gets the base map of the CSZ
  g = plotBase(proj="UTM", labels=F, countryBoundary=F)
  g2 = plotBase(proj="UTM", labels=F, countryBoundary=F)
  xy = projCSZ(lonLat, unit="m")
  
  meanSub = apply(subsidences, 1, mean)
  sdSub   = apply(subsidences, 1, sd)
  
  # Interpolate onto a regular grid
  xs = seq(min(xy[,1]), max(xy[,1]), length.out = 200)
  ys = seq(min(xy[,2]), max(xy[,2]), length.out = 400)
  meanMesh = akima::interp(xy[,1], xy[,2], meanSub, xo=xs, yo=ys)
  meanGrid = akima::interp2xyz(meanMesh, data.frame=T)
  meanGrid = na.omit(meanGrid)
  
  sdMesh = akima::interp(xy[,1], xy[,2], sdSub, xo=xs, yo=ys)
  sdGrid = akima::interp2xyz(sdMesh, data.frame=T)
  sdGrid = na.omit(sdGrid)
  
  g = g +
    geom_raster(data=meanGrid, aes(x=x, y=y, fill=z), alpha=0.75) +
    scale_fill_scico(midpoint=0, palette="vik", direction=1, limits=c(-9.5, 2.5)) +
    labs(fill="Subsidence\nMean (m)") +
    theme(legend.key.height = unit(2, 'cm'))
  
  g2 = g2 +
    geom_raster(data=sdGrid, aes(x=x, y=y, fill=z)) +
    scale_fill_gradientn(
      colours = scale_alpha_viridis(option = "magma", range=c(0.25, 1)),
      limits = c(0, 5.5)) + 
    labs(fill="Subsidence\nSD (m)") +
    theme(legend.key.height = unit(2, 'cm'))
  
  # Fix the limits
  limits = data.frame(x = -c(128, 123),
                      y = c(40, 50))
  limits = st_as_sf(x=limits,
                    coords = c("x", "y"),
                    crs=st_crs("EPSG:4326"))
  limits = toUTM(limits)
  g = g +
    coord_sf(xlim = c(st_coordinates(limits)[,1]),
             ylim = c(st_coordinates(limits)[,2]),
             crs=st_crs("EPSG:32610"))
  
  g2 = g2 +
    coord_sf(xlim = c(st_coordinates(limits)[,1]),
             ylim = c(st_coordinates(limits)[,2]),
             crs=st_crs("EPSG:32610"))
  
  plot(g)
  plot(g2)
}

plotAllSubsidenceGrid = function(AllSubsidences, lonLat, DR){
  # gets the base map of the CSZ
  g = plotBase(proj="UTM", labels=F, countryBoundary=F)
  g2 = plotBase(proj="UTM", labels=F, countryBoundary=F)
  xy = projCSZ(lonLat, unit="m")
  
  E = length(AllSubsidences)
  allQuakes = c("T1", "T2", "T3", "T4", "T5", "T6",
                "T7", "T8", "T9", "T10", "T11", "T12") 
  plotDF = data.frame()
  for (j in 1:E){
    
    # Interpolate onto a regular grid
    xs = seq(min(xy[,1]), max(xy[,1]), length.out = 100)
    ys = seq(min(xy[,2]), max(xy[,2]), length.out = 200)
    
    meanMesh = akima::interp(xy[,1], xy[,2],
                             apply(AllSubsidences[[j]], 1, mean),
                             xo=xs, yo=ys)
    meanGrid = akima::interp2xyz(meanMesh, data.frame=T)
    meanGrid = na.omit(meanGrid)
    
    sdMesh = akima::interp(xy[,1], xy[,2],
                           apply(AllSubsidences[[j]], 1, sd),
                           xo=xs, yo=ys)
    sdGrid = akima::interp2xyz(sdMesh, data.frame=T)
    sdGrid = na.omit(sdGrid)
    
    thisDF = data.frame(x     = meanGrid$x,
                        y     = meanGrid$y,
                        u     = meanGrid$z,
                        o     = sdGrid$z,
                        quake = rep(allQuakes[j], length(meanGrid$x)))
    
    plotDF = rbind(plotDF, thisDF)
  }
  # make sure events are ordered
  plotDF$quake = factor(plotDF$quake,
                        levels=allQuakes,
                        ordered=T)
  
  g = g +
    geom_raster(data=plotDF, aes(x=x, y=y, fill=u), alpha=0.75) +
    scale_fill_scico(midpoint=0, palette="vik", direction=1, limits=c(-11.5, 3)) +
    scale_x_continuous(breaks=-c(124,126)) +
    facet_wrap(~quake, nrow=2) +
    labs(fill="Subsidence\nMean (m)") +
    theme(legend.key.height = unit(2, 'cm'))
  
  g2 = g2 +
    geom_raster(data=plotDF, aes(x=x, y=y, fill=o)) +
    scale_fill_gradientn(
      colours = scale_alpha_viridis(option = "magma", range=c(0.25, 1)),
      limits = c(0, 6)) + 
    scale_x_continuous(breaks=-c(124,126)) +
    facet_wrap(~quake, nrow=2) +
    labs(fill="Subsidence\nSD (m)",
         alpha="") +
    theme(legend.key.height = unit(2, 'cm'))
  
  # Add data points
  # add the points where subsidence estimates are
  
  # Add the subsidence loaction
  xy = cbind(DR$Lon, DR$Lat)
  xy = projCSZ(xy, units="m")
  plotDF2 = data.frame(x     = xy[,1],
                       y     = xy[,2],
                       quake = DR$event,
                       site  = DR$Site,
                       sub   = DR$subsidence,
                       unc   = DR$Uncertainty)
  
  plotDF2 = plotDF2 %>%
    group_by(site, quake) %>%
    summarize(
      x     = mean(x, na.rm = TRUE),
      y     = mean(y, na.rm = TRUE),
      sub   = mean(sub, na.rm=TRUE),
      unc   = mean(unc, na.rm=TRUE),
      count = n(),
      .groups= 'drop'
    )
  
  plotDF2$quake = factor(plotDF2$quake,
                         levels = allQuakes,
                         ordered=T)
  
  g = g +
    geom_point(data=plotDF2,
               aes(x=x, y=y), colour="red")
  g2 = g2 +
    geom_point(data=plotDF2,
               aes(x=x, y=y), colour="red")
  
  # Fix the limits
  limits = data.frame(x = -c(128, 123),
                      y = c(40, 50))
  limits = st_as_sf(x=limits,
                    coords = c("x", "y"),
                    crs=st_crs("EPSG:4326"))
  limits = toUTM(limits)
  g = g +
    coord_sf(xlim = c(st_coordinates(limits)[,1]),
             ylim = c(st_coordinates(limits)[,2]),
             crs=st_crs("EPSG:32610"))
  
  g2 = g2 +
    coord_sf(xlim = c(st_coordinates(limits)[,1]),
             ylim = c(st_coordinates(limits)[,2]),
             crs=st_crs("EPSG:32610"))
  
  plot(g)
  plot(g2)
}

# Custom colour palette which scales viridis with alpha
scale_alpha_viridis = function(option = "magma", range = c(0.1, 1)){
  n = 256
  colours = viridis_pal(option = option)(n)
  alphas  = seq(range[1], range[2], length.out = n)
  
  alphaColours = mapply(function(color, alpha) scales::alpha(color, alpha),
                        colours, alphas,
                        SIMPLIFY = TRUE)
  return(alphaColours)
}

plotBase = function(scale=1, proj="lonLat", labels=TRUE, countryBoundary=TRUE){
  
  # read in the border data
  borders = readCountries()
  
  # extract each country and simplfy the geometry so it plots nicely
  canadaBorders = ms_simplify(borders$Canada, weighting=0.7, keep_shapes=TRUE)
  usBorders = ms_simplify(borders$US, weighting=0.7, keep_shapes=TRUE)
  canUSBorder = ms_simplify(borders$canUSBorder, weighting=0.7, keep_shapes=TRUE)
  
  placeNames = data.frame(Lon = c(-126),
                          Lat = c(49.9),
                          Place = c("Vancouver Island"))
  placeNames = st_as_sf(x=placeNames, coords = c("Lon", "Lat"), crs=st_crs("EPSG:4326"))
  
  stateNames = data.frame(Lon = c(-121.75, -121.5, -121, -121.5),
                          Lat = c(40.9, 44, 47.5, 50),
                          State = c("California", "Oregon", "Washington", "British Columbia"))
  stateNames = st_as_sf(x=stateNames, coords = c("Lon", "Lat"), crs=st_crs("EPSG:4326"))
  
  countryNames = data.frame(Lon=c(-119.5, -119.5),
                            Lat=c(49.2, 48.8),
                            Country=c("Canada", "USA"))
  countryNames = st_as_sf(x=countryNames, coords = c("Lon", "Lat"), crs=st_crs("EPSG:4326"))
  
  
  # Converts to UTM if needed
  if (proj == "UTM"){
    canadaBorders = toUTM(canadaBorders)
    usBorders = toUTM(usBorders)
    canUSBorder = toUTM(canUSBorder)
    placeNames = toUTM(placeNames)
    stateNames = toUTM(stateNames)
    countryNames = toUTM(countryNames)
  }
  
  
  # add land
  g = ggplot() +
    geom_sf(data = usBorders, fill = "white", colour="black") +
    geom_sf(data = canadaBorders, fill = "white", colour="black")
  
  # if the boundary boundary line is wanted
  if (countryBoundary==TRUE){
    g = g+
      geom_sf(data=canUSBorder, linewidth=scale*0.6, colour="black", linetype=11)
  }
  
  # if the labels are wanted
  if (labels){
    g = g +
      geom_sf_text(data=stateNames, aes(label=State), size=scale*2.5) +
      geom_sf_text(data=countryNames, aes(label=Country), size=scale*3, hjust=1) +
      geom_sf_text(data=placeNames, aes(label=Place), size=scale*2, angle=-40)
  }
  
  # control the appearance
  g = g +
    theme_bw() +
    theme(panel.background = element_rect('#f0feff')) +
    labs(x = "",
         y = "",
         title = "")
  
  return(g)
}


# a fucntion to get the border data for the country and required states / provinces
readCountries = function(){
  
  # reads in the Canada country data from GADM
  if (file.exists("~/Uni/NTNU/Masters Project/CSZ/R/Data/Canada/gadm/gadm41_CAN_1_pk.rds")){
    canada = readRDS("~/Uni/NTNU/Masters Project/CSZ/R/Data/Canada/gadm/gadm41_CAN_1_pk.rds")
  } else {
    canada = gadm("CAN", level=1, path="~/Uni/NTNU/Masters Project/CSZ/R/Data/Canada", resolution=2) # provinces
  }
  
  # reads in the US country data from GADM
  if (file.exists("~/Uni/NTNU/Masters Project/CSZ/R/Data/US/gadm/gadm41_USA_1_pk.rds")){
    us = readRDS("~/Uni/NTNU/Masters Project/CSZ/R/Data/US/gadm/gadm41_USA_1_pk.rds")
  } else {
    us = gadm("USA", level=1, path="~/Uni/NTNU/Masters Project/CSZ/R/Data/US", resolution=2) # states
  }
  
  usSF = st_as_sf(us, crs=st_crs("EPSG:4326"))
  canadaSF = st_as_sf(canada, crs=st_crs("EPSG:4326"))
  
  usStates = usSF[is.element(usSF$NAME_1,
                             c("California", "Oregon", "Washington",
                               "Idaho", "Nevada", "Arizona")),]
  canadaProv = canadaSF[is.element(canadaSF$NAME_1,
                                   c("British Columbia", "Alberta")),]
  
  border = read_sf('~/Uni/NTNU/Masters Project/CSZ/R/Data/Border/Canada_and_US_Border.shp')
  border2 = border[is.element(border$SectionEng, c('Straits of Georgia and Juan de Fuca',
                                                   'The 49th Parallel Boundary, Columbia Valley to Pacific Ocean',
                                                   'The 49th Parallel Boundary, Similkameen River to Columbia Valley',
                                                   'The 49th Parallel Boundary, West Kootenay to the Similkameen River')), ]
  
  return(list(Canada=canadaProv, US=usStates, canUSBorder=border2))
}

# A function that plots a given inla.mesh, with or without colours
#
# mesh        - an inla.mesh object
#             - Assumes is given in easting/ northing coordinates
# z           - colours that relate to the verticies in the inla.mesh
#               If an NA vector then just the mesh is plotted
# proj        - c("UTM", "lonLat") which projection to plot in
# legendTitle - The name of the z variable to put on colourbar legend
# colourScale - To differentiate between means and uncertainties
plotX = function(mesh, z, proj="UTM", legendTitle="Spatial Effect", colourScale="viridis"){
  
  # transform mesh object to correct projection
  if (proj=="lonLat"){
    thiscrs = "EPSG:4326"
    xy = mesh$loc[,1:2]
    latLon = projCSZ(xy, inverse=TRUE, units="km")
    mesh$loc[,1:2] = latLon
    mesh$crs = st_crs(thiscrs)
  } else{
    # coord_sf expects meters
    mesh$loc[,1:2] = 1000*mesh$loc[,1:2]
    
    thiscrs = "EPSG:32610"
    mesh$crs = st_crs(thiscrs)
  }
  
  g = plotBase(proj=proj, labels=FALSE, countryBoundary=FALSE)
  
  if (all(is.na(z))){
    g = g +
      gg(mesh, edge.colour="black", interior=FALSE, exterior=FALSE)
    
  } else{
    
    boundaryXY = as.matrix(inla.mesh.interior(mesh)[[1]]$loc[,1:2])
    boundaryXY2 = concaveman(boundaryXY)
    # Create the mask object
    points = SpatialPoints(boundaryXY2, proj4string=CRS(thiscrs))
    hullInt = SpatialPolygons(list(Polygons(list(Polygon(points)), "ExteriorBoundary")), proj4string=CRS(thiscrs))
    
    g = g +
      gg(mesh,
         colour=z,
         mask=hullInt, nx=300, ny=300) +
      scale_fill_viridis_c(alpha=0.75, name = legendTitle, option = colourScale) +
      theme(legend.position="right", legend.key.height = unit(3, 'cm'))
  }
  
  # create limits of plotting
  limits = data.frame(x = -c(128, 123),
                      y = c(40, 50))
  limits = st_as_sf(x=limits,
                    coords = c("x", "y"),
                    crs=st_crs("EPSG:4326"))
  if (proj == "UTM"){
    limits = toUTM(limits)
  }

  g = g  +
    coord_sf(xlim = c(st_coordinates(limits)[,1]),
             ylim = c(st_coordinates(limits)[,2]),
             crs=st_crs(thiscrs))
  
  return(g)
}

# A function that plots a given inla.mesh, with or without colours
#
# mesh        - an inla.mesh object
#             - Assumes is given in easting/ northing coordinates
# z           - colours that relate to the verticies in the inla.mesh
#               If an NA vector then just the mesh is plotted
# proj        - c("UTM", "lonLat") which projection to plot in
# legendTitle - The name of the z variable to put on colourbar legend
# colourScale - To differentiate between means and uncertainties
plotOneX = function(mesh, X, DR,
                    limitMean = c(NA), limitSD = c(NA),
                    event = "T1", proj="UTM", xw = "X"){
  
  
  # transform mesh object to correct projection
  if (proj=="lonLat"){
    thiscrs = "EPSG:4326"
    xy = mesh$loc[,1:2]
    latLon = projCSZ(xy, inverse=TRUE, units="km")
    mesh$loc[,1:2] = latLon
    mesh$crs = st_crs(thiscrs)
  } else{
    # coord_sf expects meters
    mesh$loc[,1:2] = 1000*mesh$loc[,1:2]
    thiscrs = "EPSG:32610"
    mesh$crs = st_crs(thiscrs)
  }
  
  g = plotBase(proj=proj, labels=FALSE, countryBoundary=FALSE)
  g2 = plotBase(proj=proj, labels=FALSE, countryBoundary=FALSE)
  
  # make the plotting boundary
  boundaryXY = as.matrix(inla.mesh.interior(mesh)[[1]]$loc[,1:2])
  hull = concaveman(boundaryXY)
  
  # Points to interpolate onto
  xs = seq(min(mesh$loc[,1]), max(mesh$loc[,1]), length.out = 1000)
  ys = seq(min(mesh$loc[,2]), max(mesh$loc[,2]), length.out = 1000)
  XY = expand.grid(xs, ys)
  
  # find interpolation points inside the hull
  inside = point.in.polygon(point.x = XY[,1],
                            point.y = XY[,2],
                            pol.x   = hull[,1],
                            pol.y   = hull[,2])
  
  # interpolate
  meanGrid = akima::interp(mesh$loc[,1], mesh$loc[,2],
                           apply(X, 1, mean),
                           xo=xs, yo=ys)
  sdGrid   = akima::interp(mesh$loc[,1], mesh$loc[,2],
                           apply(X
                                 
                                 , 1, sd),
                           xo=xs, yo=ys)
  meanGrid = akima::interp2xyz(meanGrid, data.frame=T)
  sdGrid   = akima::interp2xyz(sdGrid, data.frame=T)
  
  # just take inside points
  meanGrid = meanGrid[inside == 1,]
  sdGrid   = sdGrid[inside == 1,]
  
  plotDF = data.frame(x = meanGrid$x,
                      y = meanGrid$y,
                      mean = meanGrid$z,
                      sd = sdGrid$z)
  # plot means
  g = g +
    geom_raster(data=plotDF, aes(x=x, y=y, fill=mean)) +
    labs(fill=paste(xw,"\nMean")) +
    theme(legend.key.height = unit(2, "cm"))
  
  # plot sds
  g2 = g2 +
    geom_raster(data=plotDF, aes(x=x, y=y, fill=sd)) +
    labs(fill=paste(xw,"\nSD")) +
    theme(legend.key.height = unit(2, "cm"))
  
  # sort out colour limits
  if (all(is.na(limitMean))){
    g = g +
      scale_fill_scico(alpha=0.75, palette="cork", midpoint=0)
  } else{
    g = g +
      scale_fill_scico(alpha=0.75, palette="cork", midpoint=0, limits=limitMean)
  }
  
  if (all(is.na(limitMean))){
    g2 = g2 +
      scale_fill_viridis(alpha=0.75, option="magma")
  } else{
    g = g +
      scale_fill_viridis(alpha=0.75, option="lajolla", limits=limitSD)
  }
  
  # add the points where subsidence estimates are
  xy = cbind(DR$Lon, DR$Lat)
  xy = projCSZ(xy, units="m")
  plotDF2 = data.frame(x     = xy[DR$event == event, 1],
                       y     = xy[DR$event == event, 2],
                       site  = DR$Site[DR$event == event],
                       sub   = DR$subsidence[DR$event == event],
                       unc   = DR$Uncertainty[DR$event == event])
  
  plotDF2 = plotDF2 %>%
    group_by(site) %>%
    summarize(
      x     = mean(x, na.rm = TRUE),
      y     = mean(y, na.rm = TRUE),
      sub   = mean(sub, na.rm=TRUE),
      unc   = mean(unc, na.rm=TRUE),
      count = n()
    )
  
  g = g +
    geom_point(data=plotDF2,
               aes(x=x, y=y), colour="red")
  g2 = g2 +
    geom_point(data=plotDF2,
               aes(x=x, y=y), colour="red")
  
  # create limits of plotting
  limits = data.frame(x = -c(128, 123),
                      y = c(40, 50))
  limits = st_as_sf(x=limits,
                    coords = c("x", "y"),
                    crs=st_crs("EPSG:4326"))
  if (proj == "UTM"){
    limits = toUTM(limits)
  }
  
  g = g  +
    coord_sf(xlim = c(st_coordinates(limits)[,1]),
             ylim = c(st_coordinates(limits)[,2]),
             crs=st_crs(thiscrs))
  g2 = g2  +
    coord_sf(xlim = c(st_coordinates(limits)[,1]),
             ylim = c(st_coordinates(limits)[,2]),
             crs=st_crs(thiscrs))
  
  plot(g)
  plot(g2)
  grid.arrange(g, g2, ncol=2)
  
  
}

plotXRealisation = function(X, mesh, DR,
                            event="T1", proj="UTM", xw="X"){
  # transform mesh object to correct projection
  if (proj=="lonLat"){
    thiscrs = "EPSG:4326"
    xy = mesh$loc[,1:2]
    latLon = projCSZ(xy, inverse=TRUE, units="km")
    mesh$loc[,1:2] = latLon
    mesh$crs = st_crs(thiscrs)
  } else{
    # coord_sf expects meters
    mesh$loc[,1:2] = 1000*mesh$loc[,1:2]
    thiscrs = "EPSG:32610"
    mesh$crs = st_crs(thiscrs)
  }
  
  g = plotBase(proj=proj, labels=FALSE, countryBoundary=FALSE)
  
  # make the plotting boundary
  boundaryXY = as.matrix(inla.mesh.interior(mesh)[[1]]$loc[,1:2])
  hull = concaveman(boundaryXY)
  
  # Points to interpolate onto
  xs = seq(min(mesh$loc[,1]), max(mesh$loc[,1]), length.out = 1000)
  ys = seq(min(mesh$loc[,2]), max(mesh$loc[,2]), length.out = 1000)
  XY = expand.grid(xs, ys)
  
  # find interpolation points inside the hull
  inside = point.in.polygon(point.x = XY[,1],
                            point.y = XY[,2],
                            pol.x   = hull[,1],
                            pol.y   = hull[,2])
  
  # interpolate
  meanGrid = akima::interp(mesh$loc[,1], mesh$loc[,2],
                           X,
                           xo=xs, yo=ys)
  meanGrid = akima::interp2xyz(meanGrid, data.frame=T)
  
  # just take inside points
  meanGrid = meanGrid[inside == 1,]
  
  plotDF = data.frame(x = meanGrid$x,
                      y = meanGrid$y,
                      z = meanGrid$z)
  # plot means
  g = g +
    geom_raster(data=plotDF, aes(x=x, y=y, fill=z)) +
    scale_fill_scico(alpha=0.75, palette="cork", midpoint=0) +
    labs(fill=paste(xw,"Realisation"))
  
  
  # add the points where subsidence estimates are
  xy = cbind(DR$Lon, DR$Lat)
  xy = projCSZ(xy, units="m")
  plotDF2 = data.frame(x     = xy[DR$event == event, 1],
                       y     = xy[DR$event == event, 2],
                       site  = DR$Site[DR$event == event],
                       sub   = DR$subsidence[DR$event == event],
                       unc   = DR$Uncertainty[DR$event == event])
  
  plotDF2 = plotDF2 %>%
    group_by(site) %>%
    summarize(
      x     = mean(x, na.rm = TRUE),
      y     = mean(y, na.rm = TRUE),
      sub   = mean(sub, na.rm=TRUE),
      unc   = mean(unc, na.rm=TRUE),
      count = n()
    )
  
  g = g +
    geom_point(data=plotDF2,
               aes(x=x, y=y), colour="purple")
  
  # create limits of plotting
  limits = data.frame(x = -c(128, 123),
                      y = c(40, 50))
  limits = st_as_sf(x=limits,
                    coords = c("x", "y"),
                    crs=st_crs("EPSG:4326"))
  if (proj == "UTM"){
    limits = toUTM(limits)
  }
  
  g = g  +
    coord_sf(xlim = c(st_coordinates(limits)[,1]),
             ylim = c(st_coordinates(limits)[,2]),
             crs=st_crs(thiscrs))
  
  plot(g)
}

# A function that plots a given inla.mesh, with or without colours
#
# mesh        - an inla.mesh object
#             - Assumes is given in easting/ northing coordinates
# z           - colours that relate to the verticies in the inla.mesh
#               If an NA vector then just the mesh is plotted
# proj        - c("UTM", "lonLat") which projection to plot in
# legendTitle - The name of the z variable to put on colourbar legend
# colourScale - To differentiate between means and uncertainties
plotAllX = function(mesh, X, DR,
                    limitMean = c(NA),
                    limitSD = c(NA),
                    proj="UTM"){
  
  # transform mesh object to correct projection
  if (proj=="lonLat"){
    thiscrs = "EPSG:4326"
    xy = mesh$loc[,1:2]
    latLon = projCSZ(xy, inverse=TRUE, units="km")
    mesh$loc[,1:2] = latLon
    mesh$crs = st_crs(thiscrs)
  } else{
    # coord_sf expects meters
    mesh$loc[,1:2] = 1000*mesh$loc[,1:2]
    thiscrs = "EPSG:32610"
    mesh$crs = st_crs(thiscrs)
  }
  
  g = plotBase(proj=proj, labels=FALSE, countryBoundary=FALSE)
  g2 = plotBase(proj=proj, labels=FALSE, countryBoundary=FALSE)
  
  # make the plotting boundary
  boundaryXY = as.matrix(inla.mesh.interior(mesh)[[1]]$loc[,1:2])
  hull = concaveman(boundaryXY)
  
  # Points to interpolate onto
  xs = seq(min(mesh$loc[,1]), max(mesh$loc[,1]), length.out = 1000)
  ys = seq(min(mesh$loc[,2]), max(mesh$loc[,2]), length.out = 1000)
  XY = expand.grid(xs, ys)
  
  # find interpolation points inside the hull
  inside = point.in.polygon(point.x = XY[,1],
                            point.y = XY[,2],
                            pol.x   = hull[,1],
                            pol.y   = hull[,2])
  
  # create plotting data frame
  plotDF = data.frame()
  allQuakes = c("T1", "T2", "T3", "T4", "T5", "T6", "T7", "T8", "T9", "T10", "T11", "T12") 
  for (j in 1:length(X)){
    thisX = X[[j]] # extract simulations
    # interpolate
    meanGrid = akima::interp(mesh$loc[,1], mesh$loc[,2],
                             apply(thisX, 1, mean),
                             xo=xs, yo=ys)
    sdGrid   = akima::interp(mesh$loc[,1], mesh$loc[,2],
                             apply(thisX, 1, sd),
                             xo=xs, yo=ys)
    meanGrid = akima::interp2xyz(meanGrid, data.frame=T)
    sdGrid   = akima::interp2xyz(sdGrid, data.frame=T)
    
    # just take inside points
    meanGrid = meanGrid[inside == 1,]
    sdGrid   = sdGrid[inside == 1,]
    
    thisDF = data.frame(x = meanGrid$x,
                        y = meanGrid$y,
                        mean = meanGrid$z,
                        sd = sdGrid$z,
                        quake = rep(allQuakes[j], length(meanGrid$x)))
    
    plotDF = rbind(plotDF, thisDF)
  }

  # make sure events are ordered
  plotDF$quake = factor(plotDF$quake,
                        levels=allQuakes,
                        ordered=T)
  # plot means
  g = g +
    geom_raster(data=plotDF, aes(x=x, y=y, fill=mean)) +
    facet_wrap(~quake,
               nrow=2) +
    scale_x_continuous(breaks=-c(124,126)) +
    labs(fill="X\nMean") +
    theme(legend.key.height = unit(2, "cm"))
  
  # plot sds
  g2 = g2 +
    geom_raster(data=plotDF, aes(x=x, y=y, fill=sd)) +
    facet_wrap(~quake,
               nrow=2) +
    scale_x_continuous(breaks=-c(124,126)) +
    labs(fill="X\nSD") +
    theme(legend.key.height = unit(2, "cm"))

  # sort colour scaling
  if (all(is.na(limitMean))){
    g = g +
      scale_fill_scico(alpha=0.75, palette="cork", midpoint=0)
  } else{
    g = g +
      scale_fill_scico(alpha=0.75, palette="cork", midpoint=0, limits=limitMean)
  }
  
  if (all(is.na(limitMean))){
    g2 = g2 +
      scale_fill_viridis(alpha=0.75, option="magma")
  } else{
    g = g +
      scale_fill_viridis(alpha=0.75, palette="magma", limits=limitSD)
  }
  
  # Add the subsidence loaction
  xy = cbind(DR$Lon, DR$Lat)
  xy = projCSZ(xy, units="m")
  plotDF2 = data.frame(x     = xy[,1],
                       y     = xy[,2],
                       quake = DR$event,
                       site  = DR$Site,
                       sub   = DR$subsidence,
                       unc   = DR$Uncertainty)
  
  plotDF2 = plotDF2 %>%
    group_by(site, quake) %>%
    summarize(
      x     = mean(x, na.rm = TRUE),
      y     = mean(y, na.rm = TRUE),
      sub   = mean(sub, na.rm=TRUE),
      unc   = mean(unc, na.rm=TRUE),
      count = n(),
      .groups= 'drop'
    )
  
  plotDF2$quake = factor(plotDF2$quake,
                         levels = allQuakes,
                         ordered=T)
  g = g +
    geom_point(data=plotDF2,
               aes(x=x, y=y), colour="red")
  
  g2 = g2 +
    geom_point(data=plotDF2,
               aes(x=x, y=y), colour="red")
  
  # create limits of plotting
  limits = data.frame(x = -c(128, 123),
                      y = c(40, 50))
  limits = st_as_sf(x=limits,
                    coords = c("x", "y"),
                    crs=st_crs("EPSG:4326"))
  if (proj == "UTM"){
    limits = toUTM(limits)
  }

  g = g  +
    coord_sf(xlim = c(st_coordinates(limits)[,1]),
             ylim = c(st_coordinates(limits)[,2]),
             crs=st_crs(thiscrs))
  g2 = g2  +
    coord_sf(xlim = c(st_coordinates(limits)[,1]),
             ylim = c(st_coordinates(limits)[,2]),
             crs=st_crs(thiscrs))

  plot(g)
  plot(g2)
}

plotAllX2 = function(mesh, X, DR,
                     proj="UTM", legendTitle="Spatial Effect", colourScale="viridis"){
  
  # transform mesh object to correct projection
  if (proj=="lonLat"){
    thiscrs = "EPSG:4326"
    xy = mesh$loc[,1:2]
    latLon = projCSZ(xy, inverse=TRUE, units="km")
    mesh$loc[,1:2] = latLon
    mesh$crs = st_crs(thiscrs)
  } else{
    # coord_sf expects meters
    mesh$loc[,1:2] = 1000*mesh$loc[,1:2]
    thiscrs = "EPSG:32610"
    mesh$crs = st_crs(thiscrs)
  }
  
  g = plotBase(proj=proj, labels=FALSE, countryBoundary=FALSE)
  
  # make the plotting boundary
  boundaryXY = as.matrix(inla.mesh.interior(mesh)[[1]]$loc[,1:2])
  hull = concaveman(boundaryXY)
  
  # Points to interpolate onto
  xs = seq(min(mesh$loc[,1]), max(mesh$loc[,1]), length.out = 1000)
  ys = seq(min(mesh$loc[,2]), max(mesh$loc[,2]), length.out = 1000)
  XY = expand.grid(xs, ys)
  
  # find interpolation points inside the hull
  inside = point.in.polygon(point.x = XY[,1],
                            point.y = XY[,2],
                            pol.x   = hull[,1],
                            pol.y   = hull[,2])
  
  # create plotting data frame
  plotDF = data.frame()
  allQuakes = c("T1","T2","T3","T4","T5","T6","T7","T8","T9","T10","T11","T12") 
  for (j in 1:dim(X)[2]){
    thisX = X[,j] # extract simulations
    # interpolate
    meanGrid = akima::interp(mesh$loc[,1], mesh$loc[,2],
                             thisX,
                             xo=xs, yo=ys)
    meanGrid = akima::interp2xyz(meanGrid, data.frame=T)
    
    # just take inside points
    meanGrid = meanGrid[inside == 1,]
    
    thisDF = data.frame(x = meanGrid$x,
                        y = meanGrid$y,
                        mean = meanGrid$z,
                        quake = rep(allQuakes[j], length(meanGrid$x)))
    
    plotDF = rbind(plotDF, thisDF)
  }
  
  # make sure events are ordered
  plotDF$quake = factor(plotDF$quake,
                        levels=allQuakes,
                        ordered=T)
  # plot means
  g = g +
    geom_raster(data=plotDF, aes(x=x, y=y, fill=mean)) +
    scale_fill_scico(alpha=0.75, palette="cork", midpoint=0) +
    facet_wrap(~quake,
               nrow=2) +
    scale_x_continuous(breaks=-c(124,126)) +
    labs(fill="X\nMean")
  

  # Add the subsidence loaction
  xy = cbind(DR$Lon, DR$Lat)
  xy = projCSZ(xy, units="m")
  plotDF2 = data.frame(x     = xy[,1],
                       y     = xy[,2],
                       quake = DR$event,
                       site  = DR$Site,
                       sub   = DR$subsidence,
                       unc   = DR$Uncertainty)
  
  plotDF2 = plotDF2 %>%
    group_by(site, quake) %>%
    summarize(
      x     = mean(x, na.rm = TRUE),
      y     = mean(y, na.rm = TRUE),
      sub   = mean(sub, na.rm=TRUE),
      unc   = mean(unc, na.rm=TRUE),
      count = n(),
      .groups= 'drop'
    )
  
  plotDF2$quake = factor(plotDF2$quake,
                         levels = allQuakes,
                         ordered=T)
  g = g +
    geom_point(data=plotDF2,
               aes(x=x, y=y), colour="purple")
  
  
  # create limits of plotting
  limits = data.frame(x = -c(128, 123),
                      y = c(40, 50))
  limits = st_as_sf(x=limits,
                    coords = c("x", "y"),
                    crs=st_crs("EPSG:4326"))
  if (proj == "UTM"){
    limits = toUTM(limits)
  }
  
  g = g  +
    coord_sf(xlim = c(st_coordinates(limits)[,1]),
             ylim = c(st_coordinates(limits)[,2]),
             crs=st_crs(thiscrs))
  
  plot(g)
}

plotAllX3 = function(mesh, X, DR,
                     proj="UTM", legendTitle="Spatial Effect", colourScale="viridis"){
  
  # transform mesh object to correct projection
  if (proj=="lonLat"){
    thiscrs = "EPSG:4326"
    xy = mesh$loc[,1:2]
    latLon = projCSZ(xy, inverse=TRUE, units="km")
    mesh$loc[,1:2] = latLon
    mesh$crs = st_crs(thiscrs)
  } else{
    # coord_sf expects meters
    mesh$loc[,1:2] = 1000*mesh$loc[,1:2]
    thiscrs = "EPSG:32610"
    mesh$crs = st_crs(thiscrs)
  }
  
  g = plotBase(proj=proj, labels=FALSE, countryBoundary=FALSE)
  
  # make the plotting boundary
  boundaryXY = as.matrix(inla.mesh.interior(mesh)[[1]]$loc[,1:2])
  hull = concaveman(boundaryXY)
  
  # Points to interpolate onto
  xs = seq(min(mesh$loc[,1]), max(mesh$loc[,1]), length.out = 1000)
  ys = seq(min(mesh$loc[,2]), max(mesh$loc[,2]), length.out = 1000)
  XY = expand.grid(xs, ys)
  
  # find interpolation points inside the hull
  inside = point.in.polygon(point.x = XY[,1],
                            point.y = XY[,2],
                            pol.x   = hull[,1],
                            pol.y   = hull[,2])
  
  # create plotting data frame
  plotDF = data.frame()
  allQuakes = c("T1","T2","T3","T4","T5","T6","T7","T8","T9","T10","T11","T12") 
  for (j in 1:dim(X)[2]){
    thisX = X[,j] # extract simulations
    # interpolate
    meanGrid = akima::interp(mesh$loc[,1], mesh$loc[,2],
                             thisX,
                             xo=xs, yo=ys)
    meanGrid = akima::interp2xyz(meanGrid, data.frame=T)
    
    # just take inside points
    meanGrid = meanGrid[inside == 1,]
    
    thisDF = data.frame(x = meanGrid$x,
                        y = meanGrid$y,
                        mean = meanGrid$z,
                        quake = rep(allQuakes[j], length(meanGrid$x)))
    
    plotDF = rbind(plotDF, thisDF)
  }
  
  # make sure events are ordered
  plotDF$quake = factor(plotDF$quake,
                        levels=allQuakes,
                        ordered=T)
  # plot means
  g = g +
    geom_raster(data=plotDF, aes(x=x, y=y, fill=mean)) +
    scale_fill_scico(alpha=0.75, palette="lajolla", direction=-1) +
    facet_wrap(~quake,
               nrow=2) +
    scale_x_continuous(breaks=-c(124,126)) +
    labs(fill="X\nSD")
  
  
  # Add the subsidence loaction
  xy = cbind(DR$Lon, DR$Lat)
  xy = projCSZ(xy, units="m")
  plotDF2 = data.frame(x     = xy[,1],
                       y     = xy[,2],
                       quake = DR$event,
                       site  = DR$Site,
                       sub   = DR$subsidence,
                       unc   = DR$Uncertainty)
  
  plotDF2 = plotDF2 %>%
    group_by(site, quake) %>%
    summarize(
      x     = mean(x, na.rm = TRUE),
      y     = mean(y, na.rm = TRUE),
      sub   = mean(sub, na.rm=TRUE),
      unc   = mean(unc, na.rm=TRUE),
      count = n(),
      .groups= 'drop'
    )
  
  plotDF2$quake = factor(plotDF2$quake,
                         levels = allQuakes,
                         ordered=T)
  g = g +
    geom_point(data=plotDF2,
               aes(x=x, y=y), colour="purple")
  
  
  # create limits of plotting
  limits = data.frame(x = -c(128, 123),
                      y = c(40, 50))
  limits = st_as_sf(x=limits,
                    coords = c("x", "y"),
                    crs=st_crs("EPSG:4326"))
  if (proj == "UTM"){
    limits = toUTM(limits)
  }
  
  g = g  +
    coord_sf(xlim = c(st_coordinates(limits)[,1]),
             ylim = c(st_coordinates(limits)[,2]),
             crs=st_crs(thiscrs))
  
  plot(g)
}

# Function to plot the slab geometry
#
# fault       - like the output of getFullFaultGeom
# z           - Colour of subfaults changes with z
#             - Set to a zero vector for no colour fill
# proj        - "lonLat" or "UTM", the projection to plot in
#             - Default "UTM"
# legendTitle - The title of colour bar legend for z
# colourScale - To differentiate between means and uncertainties
#             - Options are viridis package options
plotFault = function(fault, z,
                     proj="UTM", legendTitle="Slip (m)", colourScale="viridis"){
  
  # check inputs are correct
  if (length(z) != length(fault)){
    stop("The z vector must correspond to the number of subfaults")
  }
  if (!(proj %in% c("UTM", "lonLat"))){
    stop("Not a recognised projection system")
  }
  
  # number sub faults
  nsf = length(fault)
  
  # create the dataframe to plot
  ids = c()
  lons = c()
  lats = c()
  zs = c()
  for (i in 1:nsf){
    ids = c(ids, rep(i, 4))
    lons = c(lons, fault[[i]]$corners[,1], fault[[i]]$corners[1,1])
    lats = c(lats, fault[[i]]$corners[,2], fault[[i]]$corners[1,2])
    zs = c(zs, rep(z[i], 4))
  }
  datapoly = data.frame(id=ids, lon=lons, lat=lats, z=zs)
  
  if (proj == "UTM"){
    # plot the base map
    g = plotBase(proj="UTM", labels=FALSE, countryBoundary=FALSE)
    
    # plot the subfault
    xy = cbind(datapoly$lon, datapoly$lat)
    xy = projCSZ(xy, units="m")
    datapoly$x = xy[,1]
    datapoly$y = xy[,2]
    # Add a colourscale if required
    if (!(all(is.na(datapoly$z)))){
      g = g +
        geom_polygon(data=datapoly, aes(x=x, y=y, group=id, fill=z), alpha=0.8, linewidth=NA) +
        scale_fill_viridis_c(name = legendTitle, option = colourScale)
    } else{
      g = g +
        geom_polygon(data=datapoly, aes(x=x, y=y, group=id), fill=NA, color="brown", linewidth=0.2)
    }
    
    # g = g +
    #   theme(axis.text = element_text(size = 12))
    # # sort the x axis labels
    # g = g +
    #   scale_x_continuous(breaks=-c(123,127))
    
    # Fix the limits
    limits = data.frame(x = -c(128, 123),
                        y = c(40, 50))
    limits = st_as_sf(x=limits,
                      coords = c("x", "y"),
                      crs=st_crs("EPSG:4326"))
    limits = toUTM(limits)
    
    g = g +
      coord_sf(xlim = c(st_coordinates(limits)[,1]),
               ylim = c(st_coordinates(limits)[,2]),
               crs=st_crs("EPSG:32610"))
  } else{
    # plot the base map
    g = plotBase(proj="lonLat", labels=FALSE, countryBoundary=FALSE)
    
    # Add a colourscale if required
    if (!(all(is.na(datapoly$z)))){
      g = g +
        geom_polygon(data=datapoly, aes(x=lon, y=lat, group=id, fill=z), alpha=0.8, color="black", linewidth=NA) +
        scale_fill_viridis_c(name = legendTitle, option = colourScale)
    } else{
      g = g +
        geom_polygon(data=datapoly, aes(x=lon, y=lat, group=id), fill=NA, color="brown", linewidth=0.1)
    }
    
    g = g +
      coord_sf(xlim = -c(128.5, 123),
               ylim = c(40, 50),
               crs=st_crs("EPSG:4326"))
  }
  
  g = g +
      theme(legend.position="right", legend.key.height = unit(2, 'cm'))
  
  return(g)
}

plotFutureSlip = function(fault, slips,
                          proj="UTM"){
  
  # number sub faults
  nsf = length(fault)
  
  # create the geometry dataframe to plot
  ids = c()
  lons = c()
  lats = c()
  for (i in 1:nsf){
    ids = c(ids, rep(i, 4))
    lons = c(lons, fault[[i]]$corners[,1], fault[[i]]$corners[1,1])
    lats = c(lats, fault[[i]]$corners[,2], fault[[i]]$corners[1,2])
  }
  datapoly = data.frame(id=ids, lon=lons, lat=lats)
  
  g  = plotBase(proj="UTM", labels=FALSE, countryBoundary=FALSE)
  g2 = plotBase(proj="UTM", labels=FALSE, countryBoundary=FALSE)
  g3 = plotBase(proj="UTM", labels=FALSE, countryBoundary=FALSE)
  
  # plot the subfault
  xy = cbind(datapoly$lon, datapoly$lat)
  xy = projCSZ(xy, units="m")
  datapoly$x = xy[,1]
  datapoly$y = xy[,2]
  
  # add the slip data
  datapoly$slipMean  = rep(apply(slips, 1, mean),
                           each=4)
  datapoly$slipSD    = rep(apply(slips, 1, sd),
                           each=4)
  datapoly$slipWidth = rep(apply(slips, 1, quantile, 0.975) -
                           apply(slips, 1, quantile, 0.025),
                           each=4)
  
  plotDF = datapoly
  
  # plot means
  g = g +
    geom_polygon(data=plotDF,
                 aes(x=x, y=y, group=id, fill=slipMean),
                 alpha=0.75, linewidth=NA) +
    scale_fill_viridis(name = "Future Slip\nMean (m)",
                       option = "viridis",
                       limits=c(1, 28.5)) +
  theme(legend.key.height = unit(2, "cm"))
  
  # plot sds
  g2 = g2 +
    geom_polygon(data=plotDF,
                 aes(x=x, y=y, group=id, fill=slipSD),
                 alpha=0.75, linewidth=NA) +
    scale_fill_viridis(name = "Future Slip\nSD (m)",
                       option = "magma",
                       limits=c(0.5, 19.0)) +
    theme(legend.key.height = unit(2, "cm"))
  
  # plot width
  g3 = g3 +
    geom_polygon(data=plotDF,
                 aes(x=x, y=y, group=id, fill=slipWidth),
                 alpha=0.75, linewidth=NA) +
    scale_fill_viridis(name = "Future Slip\n95% PI Width (m)",
                       option = "mako",
                       limits=c(2, 70)) +
    theme(legend.key.height = unit(2, "cm"))
  
  # Fix the limits
  limits = data.frame(x = -c(128, 123),
                      y = c(40, 50))
  limits = st_as_sf(x=limits,
                    coords = c("x", "y"),
                    crs=st_crs("EPSG:4326"))
  limits = toUTM(limits)
  g = g +
    coord_sf(xlim = c(st_coordinates(limits)[,1]),
             ylim = c(st_coordinates(limits)[,2]),
             crs=st_crs("EPSG:32610"))
  g2 = g2 +
    coord_sf(xlim = c(st_coordinates(limits)[,1]),
             ylim = c(st_coordinates(limits)[,2]),
             crs=st_crs("EPSG:32610"))
  g3 = g3 +
    coord_sf(xlim = c(st_coordinates(limits)[,1]),
             ylim = c(st_coordinates(limits)[,2]),
             crs=st_crs("EPSG:32610"))
  
  plot(g)
  plot(g2)
  grid.arrange(g, g2, ncol=2)
}


plotOneSlip = function(fault, slips, DR,
                       proj="UTM"){
  
  # number sub faults
  nsf = length(fault)
  
  # create the geometry dataframe to plot
  ids = c()
  lons = c()
  lats = c()
  for (i in 1:nsf){
    ids = c(ids, rep(i, 4))
    lons = c(lons, fault[[i]]$corners[,1], fault[[i]]$corners[1,1])
    lats = c(lats, fault[[i]]$corners[,2], fault[[i]]$corners[1,2])
  }
  datapoly = data.frame(id=ids, lon=lons, lat=lats)
  
  g  = plotBase(proj="UTM", labels=FALSE, countryBoundary=FALSE)
  g2 = plotBase(proj="UTM", labels=FALSE, countryBoundary=FALSE)
  g3 = plotBase(proj="UTM", labels=FALSE, countryBoundary=FALSE)
  
  # plot the subfault
  xy = cbind(datapoly$lon, datapoly$lat)
  xy = projCSZ(xy, units="m")
  datapoly$x = xy[,1]
  datapoly$y = xy[,2]
  
  # add the slip data
  datapoly$slipMean  = rep(apply(slips, 1, mean),
                           each=4)
  datapoly$slipSD    = rep(apply(slips, 1, sd),
                           each=4)
  datapoly$slipWidth = rep(apply(slips, 1, quantile, 0.975) -
                           apply(slips, 1, quantile, 0.025),
                           each=4)
  
  plotDF = datapoly
  
  # plot means
  g = g +
    geom_polygon(data=plotDF,
                 aes(x=x, y=y, group=id, fill=slipMean),
                 alpha=0.75, linewidth=NA) +
    scale_fill_viridis(name = "Posterior Slip\nMean (m)",
                       option = "viridis")
    theme(legend.key.height = unit(2, "cm"))
  
  # plot sds
  g2 = g2 +
    geom_polygon(data=plotDF,
                 aes(x=x, y=y, group=id, fill=slipSD),
                 alpha=0.75, linewidth=NA) +
    scale_fill_viridis(name = "Posterior Slip\nStandard Deviation (m)",
                       option = "plasma") +
    theme(legend.key.height = unit(2, "cm"))
  
  # plot width
  g3 = g3 +
    geom_polygon(data=plotDF,
                 aes(x=x, y=y, group=id, fill=slipWidth),
                 alpha=0.75, linewidth=NA) +
    scale_fill_viridis(name = "Posterior Slip\n95% PI Width (m)",
                       option = "mako") +
    theme(legend.key.height = unit(2, "cm"))
  
  # Add the subsidence location
  quake = "T1"
  allQuakes = c("T1","T2","T3","T4","T5","T6","T7","T8","T9","T10","T11","T12")
  xy = cbind(DR$Lon, DR$Lat)
  xy = projCSZ(xy, units="m")
  plotDF2 = data.frame(x     = xy[,1],
                       y     = xy[,2],
                       quake = DR$event,
                       site  = DR$Site,
                       sub   = DR$subsidence,
                       unc   = DR$Uncertainty)
  
  plotDF2 = plotDF2 %>%
    group_by(site, quake) %>%
    summarize(
      x     = mean(x, na.rm = TRUE),
      y     = mean(y, na.rm = TRUE),
      sub   = mean(sub, na.rm=TRUE),
      unc   = mean(unc, na.rm=TRUE),
      count = n(),
      .groups= 'drop'
    )
  
  plotDF2$quake = factor(plotDF2$quake,
                         levels = allQuakes,
                         ordered=T)
  g = g +
    geom_point(data=plotDF2,
               aes(x=x, y=y), colour="red")
  
  g2 = g2 +
    geom_point(data=plotDF2,
               aes(x=x, y=y), colour="red")
  
  g3 = g3 +
    geom_point(data=plotDF2,
               aes(x=x, y=y), colour="red")
  
  # Fix the limits
  limits = data.frame(x = -c(128, 123),
                      y = c(40, 50))
  limits = st_as_sf(x=limits,
                    coords = c("x", "y"),
                    crs=st_crs("EPSG:4326"))
  limits = toUTM(limits)
  g = g +
    coord_sf(xlim = c(st_coordinates(limits)[,1]),
             ylim = c(st_coordinates(limits)[,2]),
             crs=st_crs("EPSG:32610"))
  g2 = g2 +
    coord_sf(xlim = c(st_coordinates(limits)[,1]),
             ylim = c(st_coordinates(limits)[,2]),
             crs=st_crs("EPSG:32610"))
  g3 = g3 +
    coord_sf(xlim = c(st_coordinates(limits)[,1]),
             ylim = c(st_coordinates(limits)[,2]),
             crs=st_crs("EPSG:32610"))
  
  plot(g)
  plot(g2)
  plot(g3)
}

plotAllSlips = function(fault, slips, DR,
                        proj="UTM"){
  
  # number sub faults
  nsf = length(fault)
  
  # create the geometry dataframe to plot
  ids = c()
  lons = c()
  lats = c()
  for (i in 1:nsf){
    ids = c(ids, rep(i, 4))
    lons = c(lons, fault[[i]]$corners[,1], fault[[i]]$corners[1,1])
    lats = c(lats, fault[[i]]$corners[,2], fault[[i]]$corners[1,2])
  }
  datapoly = data.frame(id=ids, lon=lons, lat=lats)
  
  g  = plotBase(proj="UTM", labels=FALSE, countryBoundary=FALSE)
  g2 = plotBase(proj="UTM", labels=FALSE, countryBoundary=FALSE)
  g3 = plotBase(proj="UTM", labels=FALSE, countryBoundary=FALSE)
    
  # plot the subfault
  xy = cbind(datapoly$lon, datapoly$lat)
  xy = projCSZ(xy, units="m")
  datapoly$x = xy[,1]
  datapoly$y = xy[,2]

  # add the slip data
  earthquakes = as.factor(c("T1","T2","T3","T4","T5","T6","T7","T8","T9","T10","T11","T12"))
  plotDF = data.frame()
  for (j in 1:length(slips)){
    meanSlip  = apply(slips[[j]], 1, mean)
    sdSlip    = apply(slips[[j]], 1, sd)
    lower     = apply(slips[[j]], 1, quantile, 0.025)
    upper     = apply(slips[[j]], 1, quantile, 0.975)
    widthMean = upper - lower
    
    thisDF = datapoly
    thisDF$slipMean  = rep(meanSlip, each=4)
    thisDF$slipSD    = rep(sdSlip, each=4)
    thisDF$slipWidth = rep(widthMean, each=4)
    thisDF$quake     = rep(earthquakes[j], 4*nsf)
    
    plotDF = rbind(plotDF, thisDF)
  }
  # make sure events are ordered
  plotDF$quake = factor(plotDF$quake,
                        levels=earthquakes,
                        ordered=T)
  
  # plot means
  g = g +
    geom_polygon(data=plotDF,
                 aes(x=x, y=y, group=id, fill=slipMean),
                 alpha=0.75, linewidth=NA) +
    scale_fill_viridis(name = "Posterior Slip\nMean (m)",
                       option = "viridis",
                       limits=c(1, 34.5)) + # limits found from slips of all models
    facet_wrap(~quake, nrow=2) +
    scale_x_continuous(breaks=-c(124,126)) +
    theme(legend.key.height = unit(2, "cm"))
  
  # plot sds
  g2 = g2 +
    geom_polygon(data=plotDF,
                 aes(x=x, y=y, group=id, fill=slipSD),
                 alpha=0.75, linewidth=NA) +
    scale_fill_viridis(name = "Posterior Slip\nSD (m)",
                       option = "magma",
                       limits=c(0.5, 19)) + # limits found from slips of all models
    facet_wrap(~quake, nrow=2) +
    scale_x_continuous(breaks=-c(124,126)) +
    theme(legend.key.height = unit(2, "cm"))
  
  # plot width
  g3 = g3 +
    geom_polygon(data=plotDF,
                 aes(x=x, y=y, group=id, fill=slipWidth),
                 alpha=0.75, linewidth=NA) +
    scale_fill_viridis(name = "Posterior Slip\n95% PI Width (m)",
                       option = "mako",
                       limits=c(2.0, 71.5)) + # limits found from slips of all models
    facet_wrap(~quake, nrow=2) +
    scale_x_continuous(breaks=-c(124,126)) +
    theme(legend.key.height = unit(2, "cm"))
  
  # Add the subsidence loaction
  allQuakes = c("T1","T2","T3","T4","T5","T6","T7","T8","T9","T10","T11","T12")
  xy = cbind(DR$Lon, DR$Lat)
  xy = projCSZ(xy, units="m")
  plotDF2 = data.frame(x     = xy[,1],
                       y     = xy[,2],
                       quake = DR$event,
                       site  = DR$Site,
                       sub   = DR$subsidence,
                       unc   = DR$Uncertainty)
  
  plotDF2 = plotDF2 %>%
    group_by(site, quake) %>%
    summarize(
      x     = mean(x, na.rm = TRUE),
      y     = mean(y, na.rm = TRUE),
      sub   = mean(sub, na.rm=TRUE),
      unc   = mean(unc, na.rm=TRUE),
      count = n(),
      .groups= 'drop'
    )
  
  plotDF2$quake = factor(plotDF2$quake,
                         levels = allQuakes,
                         ordered=T)
  g = g +
    geom_point(data=plotDF2,
               aes(x=x, y=y), colour="red")
  
  g2 = g2 +
    geom_point(data=plotDF2,
               aes(x=x, y=y), colour="red")
  
  g3 = g3 +
    geom_point(data=plotDF2,
               aes(x=x, y=y), colour="red")
  
  # Fix the limits
  limits = data.frame(x = -c(128, 123),
                        y = c(40, 50))
  limits = st_as_sf(x=limits,
                      coords = c("x", "y"),
                      crs=st_crs("EPSG:4326"))
  limits = toUTM(limits)
  g = g +
      coord_sf(xlim = c(st_coordinates(limits)[,1]),
               ylim = c(st_coordinates(limits)[,2]),
               crs=st_crs("EPSG:32610"))
  g2 = g2 +
    coord_sf(xlim = c(st_coordinates(limits)[,1]),
             ylim = c(st_coordinates(limits)[,2]),
             crs=st_crs("EPSG:32610"))
  g3 = g3 +
    coord_sf(xlim = c(st_coordinates(limits)[,1]),
             ylim = c(st_coordinates(limits)[,2]),
             crs=st_crs("EPSG:32610"))
  
  plot(g)
  plot(g2)
  plot(g3)
}

plotOneMagnitude = function(magnitudes){
  data = data.frame(mag=magnitudes)
  data$lower = rep(quantile(magnitudes, 0.025), length(magnitudes))
  data$upper = rep(quantile(magnitudes, 0.975), length(magnitudes))
  
  u = mean(magnitudes)
  med = median(magnitudes)
  
  data$mean = rep(u, length(magnitudes))
  
  g = ggplot(data=data) +
    geom_histogram(aes(x=mag), bins = 25, fill="#88CCEE", colour="#88CCEE") +
    geom_vline(aes(xintercept=mean), linetype="dashed", linewidth=1, colour="#DDCC77") +
    geom_vline(aes(xintercept=lower), linetype="dashed", linewidth=1, colour="#CC6677") +
    geom_vline(aes(xintercept=upper), linetype="dashed", linewidth=1, colour="#CC6677") +
    labs(x    = "Earthquake Magnitude",
         y    = "Count")
  
  # Extract the x and y axis ranges
  plot_build = ggplot_build(g)
  x_range = range(plot_build$layout$panel_scales_x[[1]]$range$range)
  y_range = range(plot_build$layout$panel_scales_y[[1]]$range$range)

  g = g +
    annotate(geom="rect",
             xmin=-Inf, xmax=Inf,
             ymin=y_range[2]+3, ymax=Inf, fill="lightgray") +
    annotate(geom="text",
             x = (x_range[2] + x_range[1]) / 2,
             y = y_range[2]+10,
             label=paste("Mean=", round(u, 4), ", Median=", round(med, 4)))
  
  plot(g)
}

plotAllMagnitudes = function(magnitudes){
  magDF = data.frame(magnitudes)
  colnames(magDF) = paste0("T", 1:12)
  
  # convert to long format and merge
  # quakes
  allQuakes = c("T1","T2","T3","T4","T5","T6","T7","T8","T9","T10","T11","T12")
  plotDF = data.frame()
  for (j in 1:dim(magDF)[2]){
    thisMean   = round(mean(magDF[,j]), 3)
    thisMedian = round(median(magDF[,j]), 3)
    thisLower = quantile(magDF[,j], 0.025)
    thisUpper = quantile(magDF[,j], 0.975)
    
    tempDF = cbind(magDF[,j],
                   rep(paste0(allQuakes[j],
                              "\n Mean=", thisMean,
                              ", Median=", thisMedian), dim(magDF)[1]),
                   rep(thisMean, dim(magDF)[1]),
                   rep(thisLower, dim(magDF)[1]),
                   rep(thisUpper, dim(magDF)[1]))
    
    plotDF = rbind(plotDF,
                   tempDF)
  }
  colnames(plotDF) = c("mag", "quake", "mean", "lower", "upper")
  
  plotDF$quake = factor(plotDF$quake,
                        levels = unique(plotDF$quake),
                        ordered = T)
  plotDF$mag = as.numeric(plotDF$mag)
  plotDF$mean = as.numeric(plotDF$mean)
  plotDF$lower = as.numeric(plotDF$lower)
  plotDF$upper = as.numeric(plotDF$upper)
  rownames(plotDF) = NULL

  g = ggplot(data=plotDF) +
    geom_histogram(aes(x=mag), bins=25, alpha=0.75, fill="#88CCEE") +
    geom_vline(aes(xintercept=mean), linetype="dashed", linewidth=1, colour="#DDCC77") +
    geom_vline(aes(xintercept=lower), linetype="dashed", linewidth=1, colour="#CC6677") +
    geom_vline(aes(xintercept=upper), linetype="dashed", linewidth=1, colour="#CC6677") +
    facet_wrap(~quake, nrow=6) +
    labs(x    = "Earthquake Magnitude",
         y    = "Count")
  
  plot(g)
}

# A function to plot the mesh onto the base map
#
# mesh - inla.mesh object
# pts - options co-ordinates to plot
# scale - doesn't do anything atm
# proj - which projection the given data is in (should all be in same format)
plotMesh = function(mesh, pts=NULL, scale=1.5, proj="UTM"){
  
  if (proj == "lonLat"){
    g = plotBase(scale=scale, proj="lonLat", labels=FALSE, countryBoundary=FALSE)
    
    xy = mesh$loc[,1:2]
    latLon = projCSZ(xy, inverse=TRUE, units="km")
    mesh$loc[,1:2] = latLon
    
    g = g +
      gg(mesh,
         edge.color=CB_color_cycle[2], edge.linewidth=0.15,
         int.linewidth=0.75, int.color=CB_color_cycle[6],
         ext.linewidth = 0.75, ext.color = CB_color_cycle[7]) +
      coord_sf(xlim=-c(131, 120), ylim=c(36, 54))
  }
  else{
    
    g = plotBase(scale=scale, proj="UTM", labels=FALSE, countryBoundary=FALSE)
    
    # create the limits
    limits = data.frame(x = -c(131, 120),
                        y = c(36, 54))
    limits = st_as_sf(x=limits,
                      coords = c("x", "y"),
                      crs=st_crs("EPSG:4326"))
    limits = toUTM(limits)
    
    # coord_sf expects meters
    mesh$loc[,1:2] = 1000*mesh$loc[,1:2]
    
    g = g +
      gg(mesh,
         edge.color="blue", edge.linewidth=0.15,
         int.linewidth=0.75, int.color="red",
         ext.linewidth = 0.75, ext.color = "red") +
      coord_sf(xlim = c(st_coordinates(limits)[,1]),
               ylim = c(st_coordinates(limits)[,2]),
               crs=st_crs("EPSG:32610"))
  }
  
  # if want to plot the centers of the subfaults
  if (is.null(pts) == FALSE){
    if (proj == "UTM"){
      pts = projCSZ(pts, inverse=TRUE, units="km")
    }
    g = g +
      geom_point(data=data.frame(pts), aes(x=X, y=Y), size=0.1, color="red")
  }
  
  return(g)
}

# A function to plot the inla mesh and the subfault ontop of each other
#
# fault - like the output of getFullFaultGeom, describes the subfaults
#       - assumes is given in Lon/Lat coordinates
# mesh - inla.mesh object
#      - assumes is made in UTM
# scale - controls the size of text on basemap
plotBothMesh = function(fault, mesh, scale=1.5){
  
  g = plotBase(scale=scale, proj="UTM", labels=F, countryBoundary=F)
  #g = ggplot()
  # Convert inlaMesh to a polygon dataframe and plot
  
  mesh$loc[,1:2] = 1000*mesh$loc[,1:2] # coord_sf expects things in meters
  corners = mesh$loc[,1:2] # corners of the triangles, i.e. vertices
  # t: triangles, v: vertices
  tv = mesh$graph$tv
  nV = nrow(corners) # number of points in mesh
  nT = nrow(tv) # number of triangles in mesh
  
  # for each triangle, calculate its center
  allCorners = c()
  id = c()
  for(ti in 1:nT) { # loop over each subfault
    vInds = tv[ti,] # get the indicies for corners of this subfault
    thisCoords = corners[vInds,] # get the corners of the subfault
    
    allCorners = rbind(allCorners,
                       thisCoords,
                       thisCoords[1,])
    id = c(id, rep(ti, 4))
  }
  
  meshDF = data.frame(id=id,
                      x=allCorners[,1],
                      y=allCorners[,2])
  
  g = g +
    geom_polygon(data=meshDF,
                 aes(x=x, y=y, group=id, colour="blue"),
                 fill=NA, linewidth=1)
  
  # create polygon dataframe from trigeom object
  nsf = length(fault)
  ids = c()
  lons = c()
  lats = c()
  for (i in 1:nsf){
    ids = c(ids, rep(i, 4))
    lons = c(lons, fault[[i]]$corners[,1], fault[[i]]$corners[1,1])
    lats = c(lats, fault[[i]]$corners[,2], fault[[i]]$corners[1,2])
  }
  xy = cbind(lons, lats)
  xy = projCSZ(xy, units="m")
  
  subfaultDF = data.frame(
    id = ids,
    x = xy[,1],
    y = xy[,2]
  )
  
  # add subfault mesh
  g = g +
    geom_polygon(data=subfaultDF,
                 aes(x=x, y=y, group=id, colour="brown"),
                 fill=NA, linewidth=1)
  
  
  # Finally add the centroids of subfaults
  x = rep(0, nsf)
  y = rep(0, nsf)
  for (i in 1:nsf){
    x[i] = fault[[i]]$lon
    y[i] = fault[[i]]$lat
  }
  xy = cbind(x, y)
  xy = projCSZ(xy, units="m")
  
  centers = data.frame(x=xy[,1], y=xy[,2])
  
  g = g +
    geom_point(data=centers, aes(x=x, y=y, colour="#CD7F32"), size=4)
  
  # finally sort out the projection and limits
  
  limits = data.frame(x = -c(124.2, 125.2),
                      y = c(44, 45))
  limits = st_as_sf(x=limits,
                    coords = c("x", "y"),
                    crs=st_crs("EPSG:4326"))
  limits = toUTM(limits)
  
  colourLegend = c('blue'    = "SPDE Mesh",
                   'brown'  = "Subfaults",
                   '#CD7F32' = "Subfault Centroids")
  g = g +
    coord_sf(xlim = c(st_coordinates(limits)[,1]),
             ylim = c(st_coordinates(limits)[,2]),
             crs=st_crs("EPSG:32610")) +
    scale_colour_identity(labels=colourLegend, guide = "legend") +
    labs(colour="")
  
  return(g)
}






# Plots values over a triangulated fault geometry
# faultGeom: output of getFullFaultGeom
# plotVar: If null, plot the areal boundaries only. 
#          If numeric, plot plotVar values for each area.
#          If character, plot variable from faultGeom with the given name
# zlim: range of the response
# project: if FALSE, plot with lon/lat coordinates.  Otherwise, plot with projected coords 
#          using myProjection function.  This can be used when plotting the projected `easting' 
#          and `northing' variables for instance.
# cols: color vector representing the color scale
# legend.mar, legend.args, n.ticks: see ?image.plot
# plotArgs: arguments to the plot function
# scaleFun, scaleFunInverse: how to scale the color scale and its inverse. For example, log and exp
# asp: aspect ratio
# addColorBar: whether to add the color bar/legend
# forceColorsInRange: whether or not to force colors in the plotted range. Useful if 
#   you have a value outside of the range or that is NA after being transformed via the scale 
#   that you still want to plot at the edge of the color scale
# crosshatchNADensity: Adds crosshatching for areas with NA values. See ?polygon density argument.
# # min.n: approximate number of ticks in color scale. See ?pretty
# myProjection: a map projection function taking a 2 column matrix of coordinates 
#   and projects them.
# ...: arguments to polygon function
plotFaultDatTri = function(faultGeom, plotVar=NULL, zlim=NULL, project=FALSE, cols=tim.colors(), 
                           legend.mar=7, new=TRUE, plotArgs=NULL, main=NULL, xlim=NULL, xlab=NULL, scaleFun = function(x) {x}, scaleFunInverse = function(x) {x}, 
                           ylim=NULL, ylab=NULL, n.ticks=5, min.n=5, ticks=NULL, tickLabels=NULL, asp=1, legend.width=1.2, addColorBar=TRUE, 
                           legendArgs=list(), leaveRoomForLegend=TRUE, forceColorsInRange=FALSE, 
                           crosshatchNADensity=10, myProjection=NULL, ...) {
  
  polyList = lapply(faultGeom, function(x) {x$corners})
  if(!is.null(plotVar) && is.character(plotVar)) {
    # if plotVar is a variable name, extract the variable from the fault geometry
    plotVar = sapply(faultGeom, function(x) {x[[plotVar]]})
  }
  plotPolyDat(polyList=polyList, plotVar=plotVar, zlim=zlim, project=project, cols=cols, 
              legend.mar=legend.mar, new=new, plotArgs=plotArgs, main=main, xlim=xlim, xlab=xlab, scaleFun=scaleFun, scaleFunInverse=scaleFunInverse, 
              ylim=ylim, ylab=ylab, n.ticks=n.ticks, min.n=min.n, ticks=ticks, tickLabels=tickLabels, asp=asp, legend.width=legend.width, addColorBar=addColorBar, 
              legendArgs=legendArgs, leaveRoomForLegend=leaveRoomForLegend, forceColorsInRange=forceColorsInRange, 
              crosshatchNADensity=crosshatchNADensity, myProjection=myProjection, ...)
}

# Same as plotMapDat in genericSpatialPlottingFunctions, but plots values over list 
# of polygons instead of shapefile data.
# polyList: a list of matrices defining points forming polygons. If points are 3D or above, 
#           removes last coordinates to make them 2D
# plotVar: if null, plot the areal boundaries only. Else plot plotVar values for each area
# zlim: range of the response
# project: if FALSE, plot with lon/lat coordinates.  Otherwise, plot with projected coords 
#          using myProjection function.  This can be used when plotting the projected `easting' 
#          and `northing' variables for instance.
# cols: color vector representing the color scale
# legend.mar, legend.args, n.ticks: see ?image.plot
# plotArgs: arguments to the plot function
# scaleFun, scaleFunInverse: how to scale the color scale and its inverse. For example, log and exp
# asp: aspect ratio
# addColorBar: whether to add the color bar/legend
# forceColorsInRange: whether or not to force colors in the plotted range. Useful if 
#   you have a value outside of the range or that is NA after being transformed via the scale 
#   that you still want to plot at the edge of the color scale
# crosshatchNADensity: Adds crosshatching for areas with NA values. See ?polygon density argument.
# # min.n: approximate number of ticks in color scale. See ?pretty
# myProjection: a map projection function taking a 2 column matrix of coordinates 
#   and projects them.
# ...: arguments to polygon function
plotPolyDat = function(polyList, plotVar=NULL, zlim=NULL, project=FALSE, cols=tim.colors(), 
                       legend.mar=7, new=TRUE, plotArgs=NULL, main=NULL, xlim=NULL, xlab=NULL, scaleFun = function(x) {x}, scaleFunInverse = function(x) {x}, 
                       ylim=NULL, ylab=NULL, n.ticks=5, min.n=5, ticks=NULL, tickLabels=NULL, asp=1, legend.width=1.2, addColorBar=TRUE, 
                       legendArgs=list(), leaveRoomForLegend=TRUE, forceColorsInRange=FALSE, 
                       crosshatchNADensity=10, myProjection=NULL, ...) {
  require(sp)
  
  # construct mapDat from polyList:
  #  1. Make Polygon object from each polygon
  #  2. Make list of Polygon objects
  #  3. Make Polygons object from the list and ID vector
  #  4. Make SpatialPolygons object from Polygons object and a plotting order
  #  5. Make SpatialPolygonsDataFrame object from SpatialPolygons object and dataframe
  polys = list()
  for(i in 1:length(polyList)) {
    thisPolyCoords = polyList[[i]]
    thisPoly = Polygon(thisPolyCoords[,1:2]) # force points to be 2D
    thisPPoly = Polygons(list(thisPoly), ID=as.character(i))
    polys = c(polys, list(thisPPoly))
  }
  # Polys = Polygons(polys, ID=as.character(1:length(polys)))
  SpPolys = SpatialPolygons(polys, pO=as.integer(1:length(polys)))
  SpPolysDFrame = SpatialPolygonsDataFrame(SpPolys, 
                                           data.frame(list(ID=1:length(polys))), 
                                           match.ID=FALSE)
  
  # call plotMapDat
  varAreas = as.character(1:length(polys))
  regionNames = varAreas
  plotMapDat(mapDat=SpPolysDFrame, plotVar=plotVar, varAreas=varAreas, regionNames=regionNames, zlim=zlim, project=project, cols=cols, 
             legend.mar=legend.mar, new=new, plotArgs=plotArgs, main=main, xlim=xlim, xlab=xlab, scaleFun=scaleFun, scaleFunInverse=scaleFunInverse, 
             ylim=ylim, ylab=ylab, n.ticks=n.ticks, min.n=min.n, ticks=ticks, tickLabels=tickLabels, asp=asp, legend.width=legend.width, addColorBar=addColorBar, 
             legendArgs=legendArgs, leaveRoomForLegend=leaveRoomForLegend, forceColorsInRange=forceColorsInRange, 
             crosshatchNADensity=crosshatchNADensity, myProjection=myProjection, ...)
}

# old plotting functions no longer used ----

# makes plots comparing the predicted subsidence levels to the subsidence data.
# if prediction inputs are left NULL, they are computed using params
comparePredsToSubs = function(params, slipPreds=NULL, slipPredsGPS=NULL, subPreds=NULL, 
                              subPredsGPS=NULL, nsim=100, plotNameRoot="full", 
                              savePlots=TRUE, G=NULL, fileNameRoot=plotNameRoot, 
                              muVec=NULL, useGPS=FALSE, tvec=NULL, subDat=dr1, 
                              logScale=FALSE, fault=csz, latRange=c(40, 50), 
                              posNormalModel=FALSE, normalModel=posNormalModel, doGPSPred=FALSE, 
                              useMVNApprox=FALSE, taperedGPSDat=FALSE, dStar=28000, normalizeTaper=FALSE, 
                              noTitle=FALSE, lwd=.5, magLimits=c(8.6, 9.4), anisotropic=FALSE) {
  # get parameters
  if(is.null(muVec)) {
    lambdaMLE = params[1]
    muZetaMLE = params[2]
    sigmaZetaMLE = params[3]
    muXi = params[5]
    muZetaGPS = rep(muZetaMLE, nrow(slipDatCSZ))
    muZetaCSZ = rep(muZetaMLE, nrow(fault))
  }
  else {
    lambdaMLE = params[1]
    sigmaZetaMLE = params[3]
    muXi = params[5]
    muZetaMLE = muVec
    muZetaGPS = muVec[1:nrow(slipDatCSZ)]
    muZetaCSZ = muVec[(nrow(slipDatCSZ)+1):length(muVec)]
  }
  
  # get taper if necessary
  if(is.null(tvec))
    tvec = taper(getFaultCenters(fault)[,3], lambda=lambdaMLE, dStar=dStar, normalize=normalizeTaper)
  
  #
  if(taperedGPSDat)
    phiZeta = params[length(params) - anisotropic]
  else
    phiZeta = NULL
  
  # generate predictions if they are left NULL by the user
  if(is.null(slipPreds))
    slipPreds = preds(params, nsim=nsim, fault=fault, muVec=c(muZetaGPS, muZetaCSZ), tvec=tvec, 
                      posNormalModel=posNormalModel, normalModel=normalModel, phiZeta=phiZeta, 
                      taperedGPSDat=taperedGPSDat, anisotropic=anisotropic)
  if(is.null(slipPredsGPS) && doGPSPred)
    slipPredsGPS = predsGivenGPS(params, nsim=nsim, muVec=c(muZetaGPS, muZetaCSZ), fault=fault, tvec=tvec, 
                                 posNormalModel=posNormalModel, normalModel=normalModel)
  if(is.null(subPreds)) {
    if(is.null(G))
      subPreds = predsToSubsidence(params, slipPreds, useMVNApprox = useMVNApprox, subDat=subDat, 
                                   posNormalModel=posNormalModel, normalModel=normalModel, tvec=tvec)
    else
      subPreds = predsToSubsidence(params, slipPreds, G=G, useMVNApprox = useMVNApprox, subDat=subDat, 
                                   posNormalModel=posNormalModel, normalModel=normalModel, tvec=tvec)
  }
  if(is.null(subPredsGPS) && doGPSPred) {
    if(is.null(G))
      subPredsGPS = predsToSubsidence(params, slipPredsGPS, useMVNApprox = useMVNApprox, subDat=subDat, 
                                      posNormalModel=posNormalModel, normalModel=normalModel, tvec=tvec)
    else
      subPredsGPS = predsToSubsidence(params, slipPredsGPS, G=G, useMVNApprox = useMVNApprox, subDat=subDat, 
                                      posNormalModel=posNormalModel, normalModel=normalModel, tvec=tvec)
  }
  meanSlip = slipPreds$meanSlip
  meanSub = subPreds$meanSub
  u95 = subPreds$u95
  l95 = subPreds$l95
  u95Noise = subPreds$u95Noise
  l95Noise = subPreds$l95Noise
  slipSD = apply(slipPreds$slipSims, 1, sd)
  if(doGPSPred) {
    meanSlipGPS = slipPredsGPS$meanSlip
    meanSubGPS = subPredsGPS$meanSub
    u95GPS = subPredsGPS$u95
    l95GPS = subPredsGPS$l95
    u95NoiseGPS = subPredsGPS$u95Noise
    l95NoiseGPS = subPredsGPS$l95Noise
  }
  
  ##### generate mean seaDef field from Okada model
  # set Okada subsidence grid
  lonRange=c(-128, -122.5)
  nx = 300
  ny=  900
  lonGrid = seq(lonRange[1], lonRange[2], l=nx)
  latGrid = seq(latRange[1], latRange[2], l=ny)
  coordGrid = make.surface.grid(list(lon=lonGrid, lat=latGrid))
  
  ##### plot subsidence predictions
  if(savePlots)
    pdf(file=paste0(fileNameRoot, "subsidencePredictions.pdf"), width=8, height=10)
  
  # slip mean
  if(!logScale) {
    pl1 = plotFaultDat(fault, plotVar=meanSlip, main=paste0(plotNameRoot, "Mean Slip (m)"), 
                       xlim=lonRange, ylim=latRange, clab="", lwd=lwd) + geom_point(aes(x=Lon, y=Lat, col="red"), shape=3, data=subDat)
  }
  else {
    pl1 = plotFaultDat(fault, plotVar=meanSlip, main=paste0(plotNameRoot, "Mean Slip (m)"), 
                       xlim=lonRange, ylim=latRange, logScale=TRUE, clab="", lwd=lwd) + geom_point(aes(x=Lon, y=Lat, col="red"), shape=3, data=subDat)
  }
  
  # slip SE
  pl2 = plotFaultDat(fault, plotVar = slipSD, main=paste0(plotNameRoot, "Slip SD (m)"), 
                     xlim=lonRange, ylim=latRange, clab="", lwd=lwd) + geom_point(aes(x=Lon, y=Lat, col="red"), shape=3, data=subDat)
  
  # simulated subsidence data from Okada model using marginal distribution
  # subRange = range(c(-meanSub, -l95Noise, -u95Noise, subDat$subsidence))
  ord = order(subDat$Lat)
  ordDat = subDat[ord,]
  ordL95 = l95Noise[ord]
  ordU95 = u95Noise[ord]
  ordLat = subDat$Lat[ord]
  ordMeanSub = subPreds$meanSub[ord]
  # tmp = cbind(ordDat, ordL95, ordU95)
  # pl3 = ggplot() + geom_point(aes(x=subsidence, y=Lat, col="red"), shape=3, data=tmp) + 
  #   coord_fixed(xlim=subRange, ylim=latRange, expand=FALSE) + 
  #   labs(x="Subsidence (m)", y="Latitude") + geom_path(aes(x=-ordL95, y=Lat), col="blue", data=tmp) +
  #   geom_path(aes(x=-ordU95, y=Lat), col="blue", data=tmp) + guides(col=FALSE) + 
  #   ggtitle(paste0(plotNameRoot, "95% Prediction Band")) + 
  #   theme(panel.border = element_blank(), panel.grid.major = element_blank(),
  #         panel.grid.minor = element_blank(), 
  #         panel.background = element_rect(fill='white'))
  tmp = data.frame(meanSub=ordMeanSub, Lat=ordLat)
  pl3 = ggplot() + geom_point(aes(x=subsidence, y=Lat), col="red", shape=3, size=.3, data=ordDat) + 
    geom_point(aes(x=-meanSub, y=Lat), col="blue", shape=19, size=.3, data=tmp) + 
    ggtitle(paste0(plotNameRoot, "Subsidence Predictions")) + 
    scale_y_continuous("Latitude", limits=latRange) + 
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          panel.background = element_rect(fill='white')) + labs(x="Subsidence (m)", y="Latitude") + guides(shape=FALSE, fill=FALSE)
  
  # plot magnitude distribution
  mags = apply(slipPreds$slipSims, 2, getMomentFromSlip, fault=fault, dStar=dStar, normalizeTaper=normalizeTaper)
  cleanMags = mags[is.finite(mags)]
  tmp = data.frame(cleanMags=cleanMags)
  pl4 = ggplot(tmp) + geom_histogram(aes(cleanMags, y=..density..)) + labs(x="Magnitudes", y="Density") + 
    geom_vline(col="purple", xintercept=quantile(cleanMags, probs=.975), linetype=2) + 
    geom_vline(col="purple", xintercept=quantile(cleanMags, probs=.025), linetype=2) + 
    geom_vline(col="purple", xintercept=mean(cleanMags)) + 
    ggtitle(paste0(plotNameRoot, "Histogram of earthquake magnitudes")) + 
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          panel.background = element_rect(fill='white')) + 
    scale_x_continuous(limits=magLimits)
  
  if(noTitle) {
    pl1 = pl1 + ggtitle(NULL)
    pl2 = pl2 + ggtitle(NULL)
    pl3 = pl3 + ggtitle(NULL)
    pl4 = pl4 + ggtitle(NULL)
  }
  
  # combine plots into one
  multiplot(pl1, pl3, pl2, pl4, layout=matrix(1:4, ncol=2))
  
  dev.off()
  # tmp = arrangeGrob(pl1, pl3, pl2, pl4, layout_matrix=matrix(1:4, ncol=2), heights=rep(41,4), widths=rep(4,4))
  # 
  # grid.arrange(pl1, pl3, pl2, pl4, layout_matrix=matrix(1:4, ncol=2), heights=1:4, widths=1:4)
  
  invisible(NULL)
}

# plots 2x2 grid of fault plots
plotFixedSlip = function(meanSlip, medSlip=NULL, l95, u95, slipSD=NULL, plotNameRoot="full", 
                         savePlots=TRUE, fileNameRoot=plotNameRoot, logScale=FALSE, 
                         event="All", subDat=dr1) {
  if(event == "All")
    inds = 1:nrow(subDat)
  else
    inds = as.character(subDat$event) == event
  sortDat = subDat[inds,]
  
  # mean
  pl1 = plotFaultDat(csz, meanSlip, logScale=logScale, xlim=lonRange, ylim=latRange, 
                     main=paste0(plotNameRoot, " Mean Slip (m)"), clab="") + 
    geom_point(aes(x=Lon, y=Lat, col="red"), shape=3, data=sortDat)
  
  # median/sd
  if(is.null(slipSD)) {
    pl2 = plotFaultDat(csz, medSlip, logScale=logScale, xlim=lonRange, ylim=latRange, 
                       main=paste0(plotNameRoot, " Median Slip (m)"), clab="") + 
      geom_point(aes(x=Lon, y=Lat, col="red"), shape=3, data=sortDat)
  }
  else if(is.null(medSlip)) {
    pl2 = plotFaultDat(csz, slipSD, logScale=logScale, xlim=lonRange, ylim=latRange, 
                       main=paste0(plotNameRoot, " Slip SD (m)"), clab="") + 
      geom_point(aes(x=Lon, y=Lat, col="red"), shape=3, data=sortDat)
  }
  # 2.5th percentile
  pl3 = plotFaultDat(csz, l95, logScale=logScale, xlim=lonRange, ylim=latRange, 
                     main=paste0(plotNameRoot, " 2.5th Percentile Slip (m)"), clab="") + 
    geom_point(aes(x=Lon, y=Lat, col="red"), shape=3, data=sortDat)
  
  ## 97.5th percentile
  pl4 = plotFaultDat(csz, u95, logScale=logScale, xlim=lonRange, ylim=latRange, 
                     main=paste0(plotNameRoot, " 97.5th Percentile Slip (m)"), clab="") + 
    geom_point(aes(x=Lon, y=Lat, col="red"), shape=3, data=sortDat)
  
  if(savePlots)
    pdf(file=paste0(fileNameRoot, "slipDistn.pdf"), width=8, height=10)
  
  ## now put all plots together onto grid
  multiplot(pl1, pl3, pl2, pl4, cols=2)
  
  if(savePlots)
    dev.off()
}

# plot subsidence predictions against each other
compareSubs = function(params, 
                       subPreds1, subPreds2, subPreds3, subPreds4, 
                       subDat1=dr1, subDat2=subDat1, subDat3=subDat1, subDat4=subDat3, 
                       tvec=NULL, 
                       plotNameRoot1="full", plotNameRoot2="full", plotNameRoot3="full", plotNameRoot4="full", 
                       savePlots=TRUE, fileNameRoot="", 
                       logScale=FALSE, fault=csz, latRange=c(40, 50), 
                       posNormalModel=FALSE, normalModel=posNormalModel, 
                       useMVNApprox=FALSE, taperedGPSDat=FALSE, dStar=25000, 
                       normalizeTaper=FALSE, noTitle=FALSE) {
  
  # get parameters
  lambdaMLE = params[1]
  muZetaMLE = params[2]
  sigmaZetaMLE = params[3]
  muXi = params[5]
  muZetaGPS = rep(muZetaMLE, nrow(slipDatCSZ))
  muZetaCSZ = rep(muZetaMLE, nrow(fault))
  
  # get taper if necessary
  if(is.null(tvec))
    tvec = taper(getFaultCenters(fault)[,3], lambda=lambdaMLE, dStar=dStar, normalize=normalizeTaper)
  
  #
  if(taperedGPSDat)
    phiZeta = params[length(params)]
  else
    phiZeta = NULL
  
  meanSub1 = subPreds1$meanSub
  u95Noise1 = subPreds1$u95Noise
  l95Noise1 = subPreds1$l95Noise
  meanSub2 = subPreds2$meanSub
  u95Noise2 = subPreds2$u95Noise
  l95Noise2 = subPreds2$l95Noise
  meanSub3 = subPreds3$meanSub
  u95Noise3 = subPreds3$u95Noise
  l95Noise3 = subPreds3$l95Noise
  meanSub4 = subPreds4$meanSub
  u95Noise4 = subPreds4$u95Noise
  l95Noise4 = subPreds4$l95Noise
  
  ## Make plots
  
  #plot 1
  subRange = range(c(-meanSub1, -l95Noise1, -u95Noise1, subDat1$subsidence, 
                     -meanSub2, -l95Noise2, -u95Noise2, subDat2$subsidence, 
                     -meanSub3, -l95Noise3, -u95Noise3, subDat3$subsidence, 
                     -meanSub4, -l95Noise4, -u95Noise4, subDat4$subsidence))
  ord = order(subDat1$Lat)
  ordDat = subDat1[ord,]
  ordL95 = l95Noise1[ord]
  ordU95 = u95Noise1[ord]
  tmp = cbind(ordDat, ordL95, ordU95)
  pl1 = ggplot() + geom_point(aes(x=subsidence, y=Lat, col="red"), shape=3, data=tmp) + 
    coord_fixed(xlim=subRange, ylim=latRange, expand=FALSE) + 
    labs(x="Subsidence (m)", y="Latitude") + geom_path(aes(x=-ordL95, y=Lat), col="blue", data=tmp) +
    geom_path(aes(x=-ordU95, y=Lat), col="blue", data=tmp) + guides(col=FALSE) + 
    ggtitle(paste0(plotNameRoot1, "95% Prediction Band")) + 
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          panel.background = element_rect(fill='white'))
  
  #plot 2
  ord = order(subDat2$Lat)
  ordDat = subDat2[ord,]
  ordL95 = l95Noise2[ord]
  ordU95 = u95Noise2[ord]
  tmp = cbind(ordDat, ordL95, ordU95)
  pl2 = ggplot() + geom_point(aes(x=subsidence, y=Lat, col="red"), shape=3, data=tmp) + 
    coord_fixed(xlim=subRange, ylim=latRange, expand=FALSE) + 
    labs(x="Subsidence (m)", y="") + geom_path(aes(x=-ordL95, y=Lat), col="blue", data=tmp) +
    geom_path(aes(x=-ordU95, y=Lat), col="blue", data=tmp) + guides(col=FALSE) + 
    ggtitle(paste0(plotNameRoot2, "95% Prediction Band")) + 
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          panel.background = element_rect(fill='white'))
  
  #plot 3
  ord = order(subDat3$Lat)
  ordDat = subDat3[ord,]
  ordL95 = l95Noise3[ord]
  ordU95 = u95Noise3[ord]
  tmp = cbind(ordDat, ordL95, ordU95)
  pl3 = ggplot() + geom_point(aes(x=subsidence, y=Lat, col="red"), shape=3, data=tmp) + 
    coord_fixed(xlim=subRange, ylim=latRange, expand=FALSE) + 
    labs(x="Subsidence (m)", y="Latitude") + geom_path(aes(x=-ordL95, y=Lat), col="blue", data=tmp) +
    geom_path(aes(x=-ordU95, y=Lat), col="blue", data=tmp) + guides(col=FALSE) + 
    ggtitle(paste0(plotNameRoot3, "95% Prediction Band")) + 
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          panel.background = element_rect(fill='white'))
  
  #plot 4
  ord = order(subDat4$Lat)
  ordDat = subDat4[ord,]
  ordL95 = l95Noise4[ord]
  ordU95 = u95Noise4[ord]
  tmp = cbind(ordDat, ordL95, ordU95)
  pl4 = ggplot() + geom_point(aes(x=subsidence, y=Lat, col="red"), shape=3, data=tmp) + 
    coord_fixed(xlim=subRange, ylim=latRange, expand=FALSE) + 
    labs(x="Subsidence (m)", y="") + geom_path(aes(x=-ordL95, y=Lat), col="blue", data=tmp) +
    geom_path(aes(x=-ordU95, y=Lat), col="blue", data=tmp) + guides(col=FALSE) + 
    ggtitle(paste0(plotNameRoot4, "95% Prediction Band")) + 
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          panel.background = element_rect(fill='white'))
  
  # remove titles if necessary
  if(noTitle) {
    pl1 = pl1 + ggtitle(NULL)
    pl2 = pl2 + ggtitle(NULL)
    pl3 = pl3 + ggtitle(NULL)
    pl4 = pl4 + ggtitle(NULL)
  }
  
  ## save plots
  if(savePlots)
    pdf(file=paste0(fileNameRoot, "compareSubs.pdf"), width=8, height=10)
  
  # put all plots together onto grid
  multiplot(pl1, pl3, pl2, pl4, cols=2)
  
  if(savePlots)
    dev.off()
}

# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, cols=1, layout=NULL, byrow=FALSE) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols), 
                     byrow=byrow)
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


# plots nr by nc grid of fault plots
# allSlips is a list of slips to put on csz fault
# nr and nc is number of rows and columns of grid
# byrow is whether to put plots in row-major or column-major order
plotSlipGrid = function(allSlips, allPlotNames=NULL, savePlots=TRUE, 
                        fileNameRoot="", logScale=FALSE, 
                        event="All", subDat=dr1, nr=NULL, nc=2, byrow=TRUE, 
                        lwd=.5) {
  if(is.null(allPlotNames)) {
    for(i in 1:ncol(allSlips))
      allPlotNames[i] = list(NULL)
  }
  
  if(event == "All")
    inds = 1:nrow(subDat)
  else
    inds = as.character(subDat$event) == event
  sortDat = subDat[inds,]
  
  # mean
  plots = list()
  slipRange = range(allSlips)
  for(i in 1:ncol(allSlips)) {
    pl = plotFaultDat(csz, allSlips[,i], varRange=slipRange, logScale=logScale, 
                      xlim=lonRange, ylim=latRange, 
                      main=allPlotNames[[i]], clab="", lwd=lwd) + 
      geom_point(aes(x=Lon, y=Lat, col="red"), shape=3, data=sortDat)
    plots = c(plots, list(pl))
  }
  
  if(savePlots)
    pdf(file=paste0(fileNameRoot, "SlipGrid.pdf"), width=8, height=10)
  
  ## now put all plots together onto grid
  multiplot(plotlist=plots, cols=nc, byrow=byrow)
  
  if(savePlots)
    dev.off()
  else
    return(plots)
}


# function for plotting model locking normalized residuals versus latitude
plotSubsidenceResiduals = function(modelFit, tvec, subDat, G, latRange=c(40,50), 
                                   main=NULL, fault=csz) {
  # get model parameters
  params = modelFit$optPar
  muZ = params[1]
  sigmaZ = params[2]
  phiZ = params[length(params)]
  
  # compute covariance matrix of T %*% Z
  coordsZ = cbind(fault$longitude, fault$latitude)
  distMatZ = rdist.earth(coordsZ, miles=FALSE)
  corMatZ = stationary.cov(coordsZ, Covariance="Matern", theta=phiZ,
                           onlyUpper=FALSE, distMat=distMatZ, smoothness=3/2)
  covMatZ = sigmaZ^2 * corMatZ
  covMatTZ = sweep(sweep(covMatZ, 1, tvec, "*"), 2, tvec, "*")
  
  # compute marginal variances of subsidences
  covMatSubs = G %*% covMatTZ %*% t(G) + diag(subDat$Uncertainty^2)
  sigmaSubs = sqrt(diag(covMatSubs))
  
  # compute predicted subsidences
  slipPreds = muZ * tvec
  subPreds = -(G %*% slipPreds)
  resids = subDat$subsidence - subPreds
  normalizedResids = resids/sigmaSubs
  lats = subDat$Lat
  ggplot() + geom_point(aes(x=normalizedResids, y=lats), col="blue", size=.5) + 
    labs(x="Normalized Residuals", y="Latitude") + 
    geom_vline(col="black", xintercept=0) + 
    ggtitle(main) + coord_fixed(xlim=range(normalizedResids), ylim=latRange, expand=FALSE) + 
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          panel.background = element_rect(fill='white'))
}

plotFaultDat = function(rows, plotVar="depth", varRange=NULL, plotData=TRUE, 
                        logScale=FALSE, xlim=c(-128, -122), ylim=c(39.5, 50.5), 
                        xlab=NULL, ylab=NULL, main="Cascadia Subduction Zone", 
                        clab="Depth (m)", addLegend=plotData, lwd=1, 
                        xName="longitude", yName="latitude", 
                        projection=NULL, parameters=NULL, orientation=NULL, 
                        scale=1, roundDigits=2, coordsAlreadyScaled=FALSE) {
  
  if(is.null(orientation) && !is.null(projection)) {
    warning("no orientation specified, so projection being set to the last used projection")
    projection = ""
    parameters = NULL
  }
  
  # get relevant map data
  if(!is.null(projection)) {
    states <- map_data("state", projection=projection, parameters=parameters, orientation=orientation)
    west_coast <- subset(states, region %in% c("california", "oregon", "washington"))
    canada = map_data("world", "Canada", projection=projection, parameters=parameters, orientation=orientation)
  } else {
    states <- map_data("state")
    west_coast <- subset(states, region %in% c("california", "oregon", "washington"))
    canada = map_data("world", "Canada")
  }
  
  # set some reasonable plotting parameters
  if(is.null(projection)) {
    xlab = "Longitude"
    ylab = "Latitude"
  } else {
    xlab = "Easting (km)"
    ylab = "Northing (km)"
  }
  
  if(!is.data.frame(rows))
    rows = data.frame(rows)
  
  # rename row$Fault so that each fault gets unique name
  rows$Fault=1:nrow(rows)
  
  # if the user supplies a variable to plot in plotVar:
  if(!is.character(plotVar)) {
    rows$tmp = plotVar
    plotVar = "tmp"
  }
  
  # 
  if("topRightX" %in% names(rows))
    faultCornerTable = rows[,9:ncol(rows)]
  else
    faultCornerTable = NULL
  
  # make fault polygons
  faultDat = getFaultPolygons(rows, faultCornerTable=faultCornerTable)
  faultDat = merge(faultDat, rows[,c("Fault", plotVar)], by=c("Fault"))
  faultDat$plotVar=faultDat[,plotVar]
  
  # rescale coordinates of the faulty geometry if necessary
  if(coordsAlreadyScaled) {
    faultDat$longitude = faultDat$longitude / scale
    faultDat$latitude = faultDat$latitude / scale
  }
  
  # set x and y limits if necessary
  if(is.null(xlim) && is.null(ylim)) {
    xlim=c(-128, -122)
    ylim=c(39.5, 50.5)
    proj<- mapproject(xlim, ylim) # if projection unspecified, last projection is used
    xlim = proj$x
    ylim = proj$y
    xScale = scale_x_continuous(xlab, c(-.02, 0, .02, .04), labels=as.character(round(scale*c(-.02, 0, .02, .04), digits=roundDigits)), limits=c(-360, 360))
    yScale = scale_y_continuous(ylab, seq(-.74, -.60, by=.02), labels=as.character(round(scale*seq(-.74, -.60, by=.02), digits=roundDigits)), limits=c(-360, 360))
  } else if(is.null(xlim)) {
    xlim = range(c(rows[[xName]], faultDat$longitude))
  } else if(is.null(ylim)) {
    ylim = range(c(rows[[yName]], faultDat$latitude))
  } else {
    xScale = scale_x_continuous(xlab, c(-127, -125, -123), labels=c("-127", "", "-123"), limits=c(-360, 360))
  }
  
  # grey maps plot:
  if(identical(projection, "")) {
    bg = ggplot(faultDat, aes(x=longitude, y=latitude)) + 
      geom_polygon(aes(x = long, y = lat, group = group), fill = "grey", color = "black", data=canada) +
      geom_polygon(aes(x = long, y = lat, group = group), fill = "grey", color = "black", data=west_coast) + 
      coord_fixed(xlim = xlim,  ylim = ylim, expand=FALSE) + 
      labs(x=xlab, y=ylab) + xScale + yScale + 
      theme(panel.border = element_blank(), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), 
            panel.background = element_rect(fill='lightblue1'))
  } else {
    bg = ggplot(faultDat, aes(x=longitude, y=latitude)) + 
      geom_polygon(aes(x = long, y = lat, group = group), fill = "grey", color = "black", data=canada) +
      geom_polygon(aes(x = long, y = lat, group = group), fill = "grey", color = "black", data=west_coast) + 
      coord_fixed(xlim = xlim,  ylim = ylim, ratio = 1.3, expand=FALSE) + 
      labs(x=xlab, y=ylab) + xScale + 
      theme(panel.border = element_blank(), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), 
            panel.background = element_rect(fill='lightblue1'))
  }
  
  # generate fault polygons portion of plot
  if(!plotData) {
    # faultPoly = plotFault(rows, color=rgb(.2,.2,.2))
    pl = bg + geom_polygon(aes(group=factor(Fault)), color="black", size=lwd, fill=NA)
  }
  else {
    faultPoly = geom_polygon(aes(fill=plotVar, group=factor(Fault)), color="black", size=lwd)
    
    if(is.null(varRange)) {
      if(logScale)
        pl = bg + faultPoly + scale_fill_distiller(clab, palette = "Spectral", direction=-1, trans="log") + 
          ggtitle(main)
      else
        pl = bg + faultPoly + scale_fill_distiller(clab, palette = "Spectral", direction=-1) +
          ggtitle(main)
    }
    else {
      if(logScale)
        pl = bg + faultPoly + scale_fill_distiller(clab, palette = "Spectral", direction=-1, trans="log", 
                                                   limits=varRange) + 
          ggtitle(main)
      else
        pl = bg + faultPoly + scale_fill_distiller(clab, palette = "Spectral", direction=-1, 
                                                   limits=varRange) +
          ggtitle(main)
    }
    # ii <- cut(values, breaks = seq(min(values), max(values), len = 100), 
    #           include.lowest = TRUE)
    # ## Use bin indices, ii, to select color from vector of n-1 equally spaced colors
    # colors <- colorRampPalette(c("lightblue", "blue"))(99)[ii]
    # if(logScale)
    #   faultPoly = faultPoly + scale_fill_manual("", palette = "Spectral", direction=-1, trans="log")
    # else
    #   faultPoly = faultPoly + scale_fill_distiller("", palette = "Spectral", direction=-1)
  }
  
  if(addLegend)
    pl + guides(color=FALSE)
  else
    pl + guides(color=FALSE, fill=FALSE)
}

##### Plot a smooth spatial surface over google maps
plotSurfGoogleMaps = function(gpsDat, plotVar=gpsDat$Depth, varRange=NULL, plotData=TRUE, 
                              logScale=FALSE, xlim=NULL, ylim=NULL, 
                              xlab="Longitude", ylab="Latitude", main="Cascadia Subduction Zone", 
                              clab="", addLegend=TRUE, lwd=1, includeFault=TRUE, 
                              gpsHull=NULL, xName="lon", yName="lat", 
                              projection=NULL, parameters=NULL, orientation=NULL) {
  library(ggmap)
  library(ggplot2)
  library(scales)
  library(RColorBrewer)
  library(latex2exp)
  
  ### Get a map
  # https://mapstyle.withgoogle.com/
  # https://console.cloud.google.com/home/dashboard?consoleReturnUrl=https:%2F%2Fcloud.google.com%2Fmaps-platform%2Fterms%2F%3Fapis%3Dmaps%26project%3Dspatialstatisticalmapping&consoleUI=CLOUD&mods=metropolis_maps&project=spatialstatisticalmapping&organizationId=657476903663
  # map <- get_map(location = c(lon[1], lat[1], lon[2], lat[2]), zoom = 6,
  #                maptype = "terrain", source = "google")
  # style1=c(feature="administrative", element="labels", visibility="off")
  # style2=c("&style=", feature="road", element="geometry", visibility="off")
  # style3=c("&style=", feature="poi", element="labels", visibility="off")
  # style4=c("&style=", feature="landscape", element="labels", visibility="off")
  # style5=c("&style=", feature="administrative", element="geometry.stroke", color="black")
  # style6=c("&style=", feature="administrative", element="geometry.stroke", weight=.75)
  # api_key =  'AIzaSyAA5-1cW2j6Q0_-xpuWF0YB6alAXkCBsmM' # key to my google account
  # map <- get_googlemap(center=c(lon = mean(xlim), lat = mean(ylim)), zoom=5,
  #                      style=c(style1, style2, style3, style4, style5, style6), 
  #                      key=api_key)
  
  ##### Let's try with greyed out land:
  library(maps)
  library(mapdata)
  
  # choose color if necessary (lightblue1 or white)
  # fields.color.picker()
  
  # get relevant map data
  if(is.null(orientation) && !is.null(projection)) {
    warning("no orientation specified, so projection being set to the last used projection")
    projection = ""
    parameters = NULL
  }
  
  # get relevant map data
  if(!is.null(projection)) {
    states <- map_data("state", projection=projection, parameters=parameters, orientation=orientation)
    west_coast <- subset(states, region %in% c("california", "oregon", "washington"))
    canada = map_data("world", "Canada", projection=projection, parameters=parameters, orientation=orientation)
  } else {
    states <- map_data("state")
    west_coast <- subset(states, region %in% c("california", "oregon", "washington"))
    canada = map_data("world", "Canada")
  }
  
  # try plotting on an interpolated grid
  library(gstat)
  library(sp)
  library(maptools)
  library(RColorBrewer)
  
  # set x and y limits if necessary
  if(is.null(xlim) && is.null(ylim)) {
    xlim=c(-128, -122)
    ylim=c(39.5, 50.5)
    proj<- mapproject(xlim, ylim) # if projection unspecified, last projection is used
    xlim = proj$x
    ylim = proj$y
    xScale = scale_x_continuous("Easting", c(-.02, 0, .02, .04), labels=c("-.02", "0", ".02", ".04"), limits=c(-360, 360))
    ylab = "Northing"
    xlab = "Easting"
  } else if(is.null(xlim)) {
    xlim = range(gpsDat[[xName]])
  } else if(is.null(ylim)) {
    ylim = range(gpsDat[[yName]])
  } else {
    xScale = scale_x_continuous("Longitude", c(-127, -125, -123), labels=c("-127", "", "-123"), limits=c(-360, 360))
  }
  
  # plot it (choose background color with fields.color.picker())
  # bg = ggplot() + 
  #   geom_polygon(aes(x = long, y = lat, group = group), fill = "grey", color = "black", data=canada) +
  #   geom_polygon(aes(x = long, y = lat, group = group), fill = "grey", color = "black", data=west_coast) + 
  #   coord_fixed(xlim = lon,  ylim = lat, ratio = 1.3, expand=FALSE) + 
  #   theme(panel.border = element_blank(), panel.grid.major = element_blank(),
  #         panel.grid.minor = element_blank(), 
  #         panel.background = element_rect(fill='lightblue1')) + 
  #   labs(x="Longitude", y="Latitude")
  if(identical(projection, "")) {
    bg = ggplot() + 
      geom_polygon(aes(x = long, y = lat, group = group), fill = "grey", color = "black", data=canada) +
      geom_polygon(aes(x = long, y = lat, group = group), fill = "grey", color = "black", data=west_coast) + 
      coord_fixed(xlim = xlim,  ylim = ylim, ratio = 1.3, expand=FALSE) + 
      labs(x=xlab, y=ylab) + xScale + 
      theme(panel.border = element_blank(), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), 
            panel.background = element_rect(fill='lightblue1'))
  } else {
    bg = ggplot() + 
      geom_polygon(aes(x = long, y = lat, group = group), fill = "grey", color = "black", data=canada) +
      geom_polygon(aes(x = long, y = lat, group = group), fill = "grey", color = "black", data=west_coast) + 
      coord_fixed(xlim = xlim,  ylim = ylim, expand=FALSE) + 
      labs(x=xlab, y=ylab) + xScale + 
      theme(panel.border = element_blank(), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), 
            panel.background = element_rect(fill='lightblue1'))
  }
  
  # construct bounding polygon around data:
  #make concave hull around prediction mask to get final prediction points
  locking30km = gpsDat
  if(is.null(gpsHull)) {
    source("model1/seaDefAnis.r")
    if(identical(projection, "")) {
      gpsHull = ahull(locking30km[[xName]], locking30km[[yName]], alpha=.05)
    } else {
      gpsHull = ahull(locking30km[[xName]], locking30km[[yName]], alpha=2)
    }
    
  }
  indx=gpsHull$arcs[,"end1"]  
  hullPts <- cbind(locking30km[[xName]], locking30km[[yName]])[indx,] # extract the boundary points from maskXY
  
  #plot hull and data to make sure it works
  # plotSubPoly(rbind(hullPts, hullPts[1,]), cbind(locking30km$lon, locking30km$lat))
  
  # now subset prediction data frame to only be predictions within the polygon from alphahull
  library(akima)
  
  # interpolate between observations for maximum purdyness
  predGrid = make.surface.grid(list(x=seq(xlim[1], xlim[2], l=500), lat=seq(ylim[1], ylim[2], l=500)))
  preds = interpp(locking30km[[xName]], locking30km[[yName]], plotVar, predGrid[,1], predGrid[,2])
  maskFinalPreds = in.poly(cbind(preds$x, preds$y), hullPts, convex.hull=FALSE)
  preds = data.frame(preds)
  finalPreds = preds[maskFinalPreds,]
  
  if(identical(projection, "")) {
    p1 = bg + geom_tile(data = finalPreds, aes(x = x, y = y, fill = z)) + 
      scale_fill_distiller(clab, palette = "Spectral", direction=-1) + 
      geom_point(aes(x=x, y=y), pch=20, col="black", size=.1, data=locking30km) + 
      ggtitle(main) + 
      geom_polygon(aes(x = long, y = lat, group = group), fill = rgb(0,0,0,0), color = "black", data=canada) +
      geom_polygon(aes(x = long, y = lat, group = group), fill = rgb(0,0,0,0), color = "black", data=west_coast)
  } else {
    p1 = bg + geom_tile(data = finalPreds, aes(x = x, y = y, fill = z)) + 
      scale_fill_distiller(clab, palette = "Spectral", direction=-1) + 
      geom_point(aes(x=lon, y=lat), pch=20, col="black", size=.1, data=locking30km) + 
      ggtitle(main) + 
      geom_polygon(aes(x = long, y = lat, group = group), fill = rgb(0,0,0,0), color = "black", data=canada) +
      geom_polygon(aes(x = long, y = lat, group = group), fill = rgb(0,0,0,0), color = "black", data=west_coast)
  }
  
  p1
}

getFaultPolygons = function(fault=csz, faultCornerTable=NULL) {
  if(length(unique(fault$Fault)) != length(fault$Fault)) {
    fault$Fault = 1:nrow(fault)
  }
  
  calcSubfaultPolygon = function(subfault) {
    tmp = subfault
    if(!is.list(subfault)) {
      tmp = matrix(subfault, ncol=length(subfault))
      colnames(tmp) = names(subfault)
      tmp = data.frame(tmp)
    }
    subfault = tmp
    
    if(is.null(faultCornerTable))
      subFaultPoly = calcGeom(subfault)$corners[,1:2]
    else
      subFaultPoly = rbind(c(subfault$topLeftX, subfault$topLeftY), 
                           c(subfault$topRightX, subfault$topRightY), 
                           c(subfault$bottomRightX, subfault$bottomRightY), 
                           c(subfault$bottomLeftX, subfault$bottomLeftY))
    return(cbind(subfault$Fault, subFaultPoly))
  }
  
  # faultPolys = apply(fault, 1, calcSubfaultPolygon)
  subfaultList = list()
  for(i in 1:nrow(fault)) {
    subfaultList = c(subfaultList, list(fault[i,]))
  }
  faultPolys = lapply(subfaultList, calcSubfaultPolygon)
  faultPolys = do.call("rbind", faultPolys)
  faultPolys = data.frame(list(Fault=faultPolys[,1], longitude=faultPolys[,2], latitude=faultPolys[,3]))
  return(faultPolys)
}
