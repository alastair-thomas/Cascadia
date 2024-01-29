checkOptimisation = function(obj, opt0, mesh, fault, DR, earthquake){
  
  print(tail(obj$env$last.par))
  print(tail(obj$env$last.par.best))
  print(opt0)
  
  if (all(tail(obj$env$last.par, 4) == opt0$par)){
    thispar = obj$env$last.par
    print("Last parameters")
  } else if (all(tail(obj$env$last.par.best, 4) == opt0$par)){
    thispar = obj$env$last.par.best
    print("Not last parameters")
  } else{
    stop("Parameters aren't last params")
  }
  
  # get the report from the TMB object
  r = obj$report(thispar)
  
  
  # Find the effective range
  effRange = sqrt(8*1)/exp(opt0$par[3])
  print(paste("Effective Range: ", round(effRange, 2), "Km"))
  
  
  # check predicted subsidences
  # gets the base map of the CSZ
  g = plotBase(scale=2, labels=FALSE)
  
  subDF = data.frame(Lat  = DR$Lat[DR$event == earthquake],
                     Lon  = DR$Lon[DR$event == earthquake],
                     Error  = abs(DR$subsidence[DR$event == earthquake] - r$okadaSubsidence))
  
  
  g = g +
    geom_point(data=subDF, aes(x=Lon, y=Lat, colour=Error), size=2) +
    scale_colour_viridis_c(name="Subsidence Predictions\nAbsolute Error (m)", option="plasma") +
    coord_sf(xlim=-c(128, 122), ylim=c(40, 50))
  plot(g)
  
  # Check the plot of the spatial field
  g1 = plotX(mesh, z=r$x)
  plot(g1)
  
  # check average untapered slip
  g2 = plotFault(fault, z=r$untaperedSlips, legendTitle="Untapered Slips (m)")
  plot(g2)
  
  # plot tapered slips
  g3 = plotFault(fault, z=r$taperedSlips, legendTitle="Tapered Slips (m)")
  plot(g3)
}

simulateCSZ = function(SD, nSims=1000){
  # take samples from fitted model
  mu = c(SD$par.fixed, SD0$par.random)
  L = Cholesky(SD[['jointPrecision']], super = T)
  draws = rmvnorm_prec(mu = mu ,
                       chol_prec = L,
                       nSims = nSims)
  
  return(draws)
}

# simulate draws
rmvnorm_prec = function(mu, chol_prec, nSims) {
  z = matrix(rnorm(length(mu) * nSims), ncol=nSims)
  L = chol_prec #Cholesky(prec, super=TRUE)
  z = Matrix::solve(L, z, system = "Lt") ## z = Lt^-1 %*% z
  z = Matrix::solve(L, z, system = "Pt") ## z = Pt    %*% z
  z = as.matrix(z)
  mu + z
}

summariseFixedParameters = function(DF, histLabel){
  
  # Now get the histograms for the fixed parameters
  fixedDrawsLong = as.data.frame(pivot_longer(DF,
                                              cols = everything(),
                                              names_to = "Parameter",
                                              values_to = "Value"))
  
  # Create a 2x2 grid of histograms
  # cairo_pdf used so that the greek symbols are saved nicely
  g4 = ggplot(fixedDrawsLong, aes(x = Value, fill=Parameter)) +
    geom_histogram(bins=30, position = "identity", alpha = 0.8) +
    facet_wrap(~Parameter,
               scales = "free",
               labeller = as_labeller(histLabel)) +
    guides(fill = "none") +
    xlab("Parameter Value") +
    ylab("Count")
  
  plot(g4)
  # With Cairo
  #ggsave(g, filename = "Parameter Hists.pdf", 
  #       device = cairo_pdf, width=9, height=9)
  
  
  # Summarise the fixed parameters
  parSum = cbind(mean   = (apply(DF, 2, mean)),
                 median = (apply(DF, 2, median)),
                 sd     = (apply(DF, 2, sd)),
                 lower  = (apply(DF, 2, quantile, .05)),
                 upper  = (apply(DF, 2, quantile, .95)))
  print(parSum)
}

summariseRandomEffects = function(xDraws, mesh){
  # Visualise the x draws and their standard deviations
  xSum = cbind(mean   = (apply(xDraws, 1, mean)),
               median = (apply(xDraws, 1, median)),
               sd     = (apply(xDraws, 1, sd)))
  
  # plot the mean x distribution
  g5 = plotX(mesh, z=xSum[,1], legendTitle="Spatial Random Effect\nMean")
  plot(g5)
  
  # plot the standard deviation of x distribution
  g6 = plotX(mesh, z=xSum[,3], legendTitle="Spatial Random Effect\nStandard Deviation", colourScale="plasma")
  plot(g6)
}
