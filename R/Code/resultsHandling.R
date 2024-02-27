summariseMagnitudes = function(magnitudes){
  data = data.frame(M=magnitudes)
  q1 = quantile(magnitudes, 0.025)
  q2 = quantile(magnitudes, 0.975)
  
  g = ggplot(data = data, aes(x = M)) +
        geom_histogram(binwidth = 0.05, fill="skyblue") +
        geom_vline(xintercept = q1, color = "red", linetype = "dashed", linewidth=1) +
        geom_vline(xintercept = q2, color = "red", linetype = "dashed", linewidth=1) +
        labs(x = "Earthquake Magnitude",
             y = "Count")
  
  return(g)
}

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

# Check the output of the multivariate optimisation
checkOptimisationMulti = function(obj, opt0, mesh, fault, DR, plotAll=F){
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
  
  lambda = exp(opt0$par[1])
  mu = opt0$par[2]
  kappa = exp( opt0$par[3])
  tau =  exp(opt0$par[4])
  
  X = r$X
  
  # Find the effective range
  effRange = sqrt(8*1)/exp(kappa)
  print(paste("Effective Range: ", round(effRange, 2), "Km"))
  
  marVar = 1.0 / (4.0*3.14159265359*(kappa^2)*(tau^2))
  print(paste("Marginal Variance: ", round(marVar, 2)))
  
  M = dim(X)[2] # number of earthquakes
  B = dim(X)[1] # number of basis functions
  
  # Check that X is different between earthquakes
  plotDF = data.frame(quake=as.factor(rep(1:M, each=B)),
                      x = c(X),
                      triangle = rep(1:B, M))
  
  g = ggplot(plotDF, aes(x=triangle, y=x, colour=quake))+
    geom_line()
  plot(g)
  rm(g)
  
  ## extract the depths of the centers of each subfault
  K = length(fault)
  depths = rep(0, K) # empty vector for depths
  for (i in 1:K){
    # Use positive depths
    # as depth increases, taper goes from one to zero
    depths[i] = -fault[[i]]$depth / (10^3) # change to kms
  }
  
  taper = exp(-depths*lambda)
  
  g = plotFault(fault, z=taper, legendTitle="Taper")
  plot(g)
  rm(g)
  
  if (plotAll){
    # Check the underlying distributions
    earthquakes = c("T1","T2","T3","T4","T5","T6","T7","T8","T9","T10","T11","T12")
    for (i in 1:M){
      thisX = X[,i]
      g = plotX(mesh, z=thisX)
      plotDF = data.frame(x=DR$Lon[DR$event == earthquakes[i]],
                          y=DR$Lat[DR$event == earthquakes[i]])
      g = g +
        geom_point(data=plotDF, aes(x=x, y=y), colour="red") +
        ggtitle(paste("Earthquake ", i))
      
      plot(g)
      rm(g)
      
      thisSlip = taper * exp(mu + thisX[inla_mesh$idx$loc])
      # plot tapered slips
      g = plotFault(fault, z=thisSlip, legendTitle="Tapered Slips (m)")
      g = g +
        ggtitle(paste("Earthquake ", i))
      plot(g)
      rm(g)
    }
  } else{
    # Just plot the 1700 quake
    
    thisX = X[,1]
    g = plotX(mesh, z=thisX)
    plotDF = data.frame(x=DR$Lon[DR$event == earthquakes[1]],
                        y=DR$Lat[DR$event == earthquakes[1]])
    g = g +
      geom_point(data=plotDF, aes(x=x, y=y), colour="red") +
      ggtitle("1700 Earthquake")
    
    plot(g)
    rm(g)
    
    thisSlip = taper * exp(mu + thisX[inla_mesh$idx$loc])
    # plot tapered slips
    g = plotFault(fault, z=thisSlip, legendTitle="Tapered Slips (m)")
    g = g +
      ggtitle("1700 Earthquake")
    plot(g)
    rm(g)
  }
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

summariseRandomEffects = function(xDraws, mesh, quake="1"){
  # Visualise the x draws and their standard deviations
  xSum = cbind(mean   = (apply(xDraws, 1, mean)),
               median = (apply(xDraws, 1, median)),
               sd     = (apply(xDraws, 1, sd)))
  
  # plot the mean x distribution
  g5 = plotX(mesh, z=xSum[,1], legendTitle="Spatial Random Effect\nMean")
  g5 = g5 + ggtitle(paste("Earthquake " + quake))
  plot(g5)
  
  # plot the standard deviation of x distribution
  g6 = plotX(mesh, z=xSum[,3], legendTitle="Spatial Random Effect\nStandard Deviation", colourScale="plasma")
  g6 = g6 + ggtitle(paste("Earthquake " + quake))
  plot(g6)
}
