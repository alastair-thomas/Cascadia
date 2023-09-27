
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

##### functions for plotting multiple fields on a grid

plotSubsidenceGrid = function(allSubs, allPlotNames=NULL, savePlots=TRUE, fileNameRoot="", 
                              event="All", subDat=dr1, nr=NULL, nc=2, byrow=TRUE, 
                              latRange=c(40,50)) {
  if(is.null(allPlotNames)) {
    for(i in 1:ncol(allSubs))
      allPlotNames[i] = list(NULL)
  }
  
  if(event == "All")
    inds = 1:nrow(subDat)
  else
    inds = as.character(subDat$event) == event
  sortDat = subDat[inds,]
  
  # mean
  plots = list()
  subRange = range(c(allSubs, sortDat$subsidence))
  for(i in 1:ncol(allSubs)) {
    sortDat$theseSubs = allSubs[,i]
    pl = ggplot() + 
      geom_point(aes(x=subsidence, y=Lat), col="red", shape=3, data=sortDat) +
      geom_point(aes(x=theseSubs, y=Lat), col="blue", size=.3, data=sortDat) + 
      coord_fixed(xlim=subRange, ylim=latRange, expand=FALSE) + 
      theme(panel.border = element_blank(), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), 
            panel.background = element_rect(fill='white')) + 
      guides(color=FALSE) + 
      labs(x="Subsidence (m)", y="Latitude")
    plots = c(plots, list(pl))
  }
  
  if(savePlots)
    pdf(file=paste0(fileNameRoot, "SubGrid.pdf"), width=8, height=10)
  
  ## now put all plots together onto grid
  multiplot(plotlist=plots, cols=nc, byrow=byrow)
  
  if(savePlots)
    dev.off()
  
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

# plot all subfaults in the fault using ggplot.  Don't add data, returns the geom_polygon object
plotFault = function(fault=csz, color="black") {
  faultPolys = getFaultPolygons(fault)
  geom_polygon(aes(x=longitude, y=latitude, group=factor(Fault)), data=faultPolys, 
               fill=rgb(1,1,1, 0), color=color)
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
