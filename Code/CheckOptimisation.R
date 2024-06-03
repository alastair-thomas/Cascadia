plotMagnitudes = function(magnitudes){
  data = data.frame(M=magnitudes)
  q1 = quantile(magnitudes, 0.025)
  q2 = quantile(magnitudes, 0.975)
  u = mean(magnitudes)
  
  g = ggplot(data = data, aes(x = M)) +
        geom_histogram(binwidth = 0.01, fill="skyblue") +
        geom_vline(xintercept = q1, color = "red", linetype = "dashed", linewidth=1) +
        geom_vline(xintercept = q2, color = "red", linetype = "dashed", linewidth=1) +
        geom_vline(xintercept = u, color = "black", linetype = "dashed", linewidth=1) +
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

checkOptimisationAniso = function(obj, opt0, allPar, mesh, fault, DR, earthquake){
  
  print(tail(obj$env$last.par), 5)
  print(tail(obj$env$last.par.best), 5)
  print(opt0)
  
  if (all(tail(obj$env$last.par, 5) == opt0$par)){
    thispar = obj$env$last.par
    print("Last parameters")
  } else if (all(tail(obj$env$last.par.best, 5) == opt0$par)){
    thispar = obj$env$last.par.best
    print("Not last parameters")
  } else{
    stop("Parameters aren't last params")
  }
  
  # get the report from the TMB object
  r = obj$report(thispar)
  
  # extract fitted best parameters
  lambda = exp(opt0$par[1])
  mu     = opt0$par[2]
  kappa  = exp(opt0$par[3])
  tau    = exp(opt0$par[4])
  psi    = exp(opt0$par[5])
  
  # plot the iterations of parameters
  I = length(allPar[,1])
  P = length(allPar[1,])
  for (j in 1:P){
    plot(x = c(1:I), y= allPar[,j],
         type='o',
         col = "blue",
         xlab = "Iteration",
         ylab = colnames(allPar)[j])
  }
  
  # Checking the anisotropy
  print(r$H)
  
  E = eigen(r$H)
  evals = E$values
  print(paste("Eigenvalue 1", evals[1]))
  print(paste("Eigenvalue 2", evals[2]))
  
  print(paste("H Determinant: ", det(r$H)))
  
  unitCircle = cbind(cos(seq(0, 2*pi, length.out = 100)),
                     sin(seq(0, 2*pi, length.out = 100)))
  plot(unitCircle, type="l", xlim=c(-2,2), ylim=c(-2,2))
  
  unitCircleH = t(r$H %*% t(unitCircle))
  lines(unitCircleH, type="l", col="blue")
  
  H2 = E$vectors
  unitCircleH2 = t(H2 %*% t(unitCircle))
  lines(unitCircleH, type="l", col="green")
  
  effRangeEasting = evals[1] * (sqrt(8*1)/kappa)
  print(paste("Effective Range New Easting: ", round(effRangeEasting, 2), "Km"))
  
  effRangeNorthing = evals[2] * (sqrt(8*1)/kappa)
  print(paste("Effective Range New Northing: ", round(effRangeNorthing, 2), "Km"))
  
  # check predicted subsidences
  # gets the base map of the CSZ
  g = plotBase(scale=2, labels=FALSE)
  
  subDF = data.frame(Lat  = DR$Lat[DR$event == earthquake],
                     Lon  = DR$Lon[DR$event == earthquake],
                     Error  = DR$subsidence[DR$event == earthquake] - r$okadaSubsidence)
  
  
  g = g +
    geom_point(data=subDF, aes(x=Lon, y=Lat, colour=Error), size=2) +
    scale_colour_viridis_c(name="Subsidence Predictions\n Error (m)", option="plasma") +
    coord_sf(xlim=-c(128, 122), ylim=c(40, 50))
  plot(g)
  
  # Check the plot of the spatial field
  g1 = plotX(mesh, proj="northing", z=r$x)
  plot(g1)
  
  # check average untapered slip
  g2 = plotFault(fault, z=r$untaperedSlips, legendTitle="Untapered Slips (m)")
  plot(g2)
  
  # plot tapered slips
  g3 = plotFault(fault, z=r$taperedSlips, legendTitle="Tapered Slips (m)")
  plot(g3)
  
  # Check the contribution of each part
  plotDF = data.frame(x=c("nll1", "nll2", "nll3", "nll4", "nll5"),
                      y=c(r$nll1, r$nll2, r$nll3, r$nll4, r$nll5))
  g4 = ggplot(plotDF, aes(x=x,y=y))+
        geom_bar(stat="identity")
  plot(g4)
}

# Check the output of the multivariate optimisation
checkOptimisationModel1 = function(obj, opt0, allPar, mesh, fault, DR, plotAll=F){
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
  
  # some nice colour blind friendly colours for plotting
  CB_color_cycle = c('#377eb8', '#ff7f00', '#4daf4a',
                     '#f781bf', '#a65628', '#984ea3',
                     '#999999', '#e41a1c', '#dede00')
  
  # plot the iterations of parameters
  I = length(allPar[,1])
  for (j in 1:4){
    plot(x = c(1:I), y= allPar[,j],
         type='o',
         col = CB_color_cycle[j],
         xlab = "Iteration",
         ylab = colnames(allPar)[j])
  }
  
  # get the report from the TMB object
  r = obj$report(thispar)
  
  # extract useful parameters
  lambda = exp(opt0$par[1])
  mu = opt0$par[2]
  kappa = exp(opt0$par[3])
  tau =  exp(opt0$par[4])
  
  X = r$X
  
  # Find the effective range
  effRange = sqrt(8*1)/kappa
  print(paste("Effective Range: ", round(effRange, 2), "Km"))
  
  marVar = 1.0 / (4.0*3.14159265359*(kappa^2)*(tau^2))
  print(paste("Marginal Variance: ", round(marVar, 2)))
  
  E = dim(X)[2] # number of earthquakes
  B = dim(X)[1] # number of basis functions
  
  # Check that X is different between earthquakes
  plotDF = data.frame(quake=as.factor(rep(1:E, each=B)),
                      x = c(X),
                      triangle = rep(1:B, E))
  
  g = ggplot(plotDF, aes(x=triangle, y=x, colour=quake))+
    geom_line()
  plot(g)
  rm(g)
  
  ## extract the depths of the centers of each subfault
  S = length(fault)
  depths = rep(0, S) # empty vector for depths
  for (i in 1:S){
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
    for (i in 1:E){
      thisX = X[,i]
      g = plotX(mesh, z=thisX)
      plotDF = data.frame(x=DR$Lon[DR$event == earthquakes[i]],
                          y=DR$Lat[DR$event == earthquakes[i]])
      g = g +
        geom_point(data=plotDF, aes(x=x, y=y), colour="red") +
        ggtitle(paste("Earthquake ", i))
      
      plot(g)
      rm(g)
      
      thisSlip = taper * exp(mu + thisX[mesh$idx$loc])
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
    
    thisSlip = taper * exp(mu + thisX[mesh$idx$loc])
    # plot tapered slips
    g = plotFault(fault, z=thisSlip, legendTitle="Tapered Slips (m)")
    g = g +
      ggtitle("1700 Earthquake")
    plot(g)
    rm(g)
  }
}

# Check the output of the multiQuake earthquake
# The different tau and kappa models
checkOptimisationModel2 = function(obj, opt0, allPar, mesh, fault, DR, plotAll=F){
  print("Last Parameters:")
  print(tail(obj$env$last.par, 6))
  print("Best Parameters")
  print(tail(obj$env$last.par.best, 6))
  print(opt0)
  
  if (all(tail(obj$env$last.par, 6) == opt0$par)){
    thispar = obj$env$last.par
    print("Last parameters")
  } else if (all(tail(obj$env$last.par.best, 6) == opt0$par)){
    thispar = obj$env$last.par.best
    print("Not last parameters")
  } else{
    stop("Parameters aren't last params")
  }
  
  # some nice colours for plotting
  CB_color_cycle = c('#377eb8', '#ff7f00', '#4daf4a',
                     '#f781bf', '#a65628', '#984ea3',
                     '#999999', '#e41a1c', '#dede00')
  
  # plot the iterations of parameters
  I = length(allPar[,1])
  for (j in 1:6){
    plot(x = c(1:I), y= allPar[,j],
         type='o',
         col = CB_color_cycle[j],
         xlab = "Iteration",
         ylab = colnames(allPar)[j])
  }
  
  # get the report from the TMB object
  r = obj$report(thispar)
  
  print(paste("Total nll: ", r$nll))
  print(paste("X prior nll1: ", r$nll1))
  print(paste("W prior nll2: ", r$nll2))
  print(paste("Taper prior nll3: ", r$nll3))
  print(paste("Spde nll4: ", r$nll4))
  print(paste("Data nll5: ", r$nll5))
  
  lambda = exp(opt0$par[1])
  mu = opt0$par[2]
  kappaX = exp(opt0$par[3])
  tauX =  exp(opt0$par[4])
  kappaW = exp(opt0$par[5])
  tauW =  exp(opt0$par[6])
  
  X = r$X
  w = r$w
  
  # Find the effective range
  effRangeX = sqrt(8*1)/kappaX
  print(paste("Effective Range X: ", round(effRangeX, 2), "Km"))
  
  marVarX = 1.0 / (4.0*3.14159265359*(kappaX^2)*(tauX^2))
  print(paste("Marginal Variance X: ", round(marVarX, 2)))
  
  effRangeW = sqrt(8*1)/kappaW
  print(paste("Effective Range W: ", round(effRangeW, 2), "Km"))
  
  marVarW = 1.0 / (4.0*3.14159265359*(kappaW^2)*(tauW^2))
  print(paste("Marginal Variance W: ", round(marVarW, 2)))
  
  E = dim(X)[2] # number of earthquakes
  B = dim(X)[1] # number of basis functions
  
  # Check that X is different between earthquakes
  plotDF = data.frame(Quake=as.factor(rep(1:E, each=B)),
                      x = c(X),
                      triangle = rep(1:B, E))
  g = ggplot(plotDF, aes(x=triangle, y=x, colour=Quake))+
    geom_line()
  plot(g)
  rm(g)
  
  ## extract the depths of the centers of each subfault
  S = length(fault)
  depths = rep(0, S) # empty vector for depths
  for (i in 1:S){
    # Use positive depths
    # as depth increases, taper goes from one to zero
    depths[i] = -fault[[i]]$depth / (10^3) # change to kms
  }
  
  # Calculate the taper
  taper = exp(-depths*lambda)
  
  # Plot the taper
  g = plotFault(fault, z=taper, legendTitle="Taper")
  plot(g)
  rm(g)
  
  # Plot the shared effect W
  g = plotX(mesh, z=w, legendTitle="Spatial Effect W")
  plot(g)
  rm(g)
  
  if (plotAll){
    # Check the underlying distributions
    earthquakes = c("T1","T2","T3","T4","T5","T6","T7","T8","T9","T10","T11","T12")
    
    for (i in 1:E){
      thisX = X[,i]
      g = plotX(mesh, z=thisX, legendTitle="Spatial Effect X")
      plotDF = data.frame(x=DR$Lon[DR$event == earthquakes[i]],
                          y=DR$Lat[DR$event == earthquakes[i]])
      g = g +
        geom_point(data=plotDF, aes(x=x, y=y), colour="red") +
        ggtitle(paste("Earthquake ", i))
      
      plot(g)
      rm(g)
      
      thisSlip = taper * exp(mu + w[mesh$idx$loc] + thisX[mesh$idx$loc])
      # plot tapered slips
      g = plotFault(fault, z=thisSlip, legendTitle="Tapered Slips (m)")
      g = g +
        ggtitle(paste("Earthquake ", i))
      plot(g)
      rm(g)
    }
  } else{
    # Just plot the 1700 quake
    
    # extract 1700 quake
    thisX = X[,1]
    g = plotX(mesh, z=thisX, legendTitle="Spatial Effect X")
    plotDF = data.frame(x=DR$Lon[DR$event == earthquakes[1]],
                        y=DR$Lat[DR$event == earthquakes[1]])
    g = g +
      geom_point(data=plotDF, aes(x=x, y=y), colour="red") +
      ggtitle("1700 Earthquake")
    plot(g)
    rm(g)
    
    thisSlip = taper * exp(mu + w[mesh$idx$loc] + thisX[mesh$idx$loc])
    # plot tapered slips
    g = plotFault(fault, z=thisSlip, legendTitle="Tapered Slips (m)")
    g = g +
      ggtitle("1700 Earthquake")
    plot(g)
    rm(g)
  }
}

# Check the output of the multivariate optimisation
checkOptimisationMultiAniso = function(obj, opt0, mesh, fault, DR, allPar, plotAll=F){
  print(tail(obj$env$last.par))
  print(tail(obj$env$last.par.best))
  print(opt0)
  
  if (all(tail(obj$env$last.par, 6) == opt0$par)){
    thispar = obj$env$last.par
    print("Last parameters")
  } else if (all(tail(obj$env$last.par.best, 6) == opt0$par)){
    thispar = obj$env$last.par.best
    print("Not last parameters")
  } else{
    stop("Parameters aren't last params")
  }
  
  # some nice colour blind friendly colours for plotting
  CB_color_cycle = c('#377eb8', '#ff7f00', '#4daf4a',
                     '#f781bf', '#a65628', '#984ea3',
                     '#999999', '#e41a1c', '#dede00')
  
  # plot the iterations of parameters
  I = length(allPar[,1]) # number of iterations
  P = length(allPar[1,]) # number of parameters
  for (j in 1:P){
    plot(x = c(1:I), y= allPar[,j],
         type='o',
         col = CB_color_cycle[j],
         xlab = "Iteration",
         ylab = colnames(allPar)[j])
  }
  
  # get the report from the TMB object
  r = obj$report(thispar)
  
  # extract useful parameters
  lambda = exp(opt0$par[1])
  mu = opt0$par[2]
  kappa = exp(opt0$par[3])
  tau =  exp(opt0$par[4])
  v1 = opt0$par[5]
  v2 = opt0$par[6]
  
  X = r$X
  
  # Find the effective range
  effRange = sqrt(8*1)/kappa
  print(paste("Effective Range: ", round(effRange, 2), "Km"))
  
  marVar = 1.0 / (4.0*3.14159265359*(kappa^2)*(tau^2))
  print(paste("Marginal Variance: ", round(marVar, 2)))
  
  # Check the Anisotropy
  plot(x=c(0,v1), y=c(0,v2))
  
  
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

# Check the output of the multiQuake earthquake
# The shared Phi model
checkOptimisationShared = function(obj, opt0, allPar, mesh, fault, DR, plotAll=F){
  print("Last Parameters:")
  print(tail(obj$env$last.par, 5))
  print("Best Parameters")
  print(tail(obj$env$last.par.best, 5))
  print(opt0)
  
  if (all(tail(obj$env$last.par, 5) == opt0$par)){
    thispar = obj$env$last.par
    print("Last parameters")
  } else if (all(tail(obj$env$last.par.best, 5) == opt0$par)){
    thispar = obj$env$last.par.best
    print("Not last parameters")
  } else{
    stop("Parameters aren't last params")
  }
  
  # some nice colours for plotting
  CB_color_cycle = c('#377eb8', '#ff7f00', '#4daf4a',
                     '#f781bf', '#a65628', '#984ea3',
                     '#999999', '#e41a1c', '#dede00')
  
  # plot the iterations of parameters
  I = length(allPar[,1])
  for (j in 1:5){
    plot(x = c(1:I), y= allPar[,j],
         type='o',
         col = CB_color_cycle[j],
         xlab = "Iteration",
         ylab = colnames(allPar)[j])
  }
  
  # get the report from the TMB object
  r = obj$report(thispar)
  
  print(paste("Total nll: ", r$nll))
  print(paste("X prior nll1: ", r$nll1))
  print(paste("Taper prior nll2: ", r$nll2))
  print(paste("Phi prior nll3: ", r$nll3))
  print(paste("Spde nll4: ", r$nll4))
  print(paste("Data nll5: ", r$nll5))
  
  phi = exp(opt0$par[1]) / (exp(opt0$par[1]) + 1)
  lambda = exp(opt0$par[2])
  mu = opt0$par[3]
  kappa = exp(opt0$par[4])
  tau =  exp(opt0$par[5])
  
  X = r$X
  w = r$w
  
  # Find the effective range
  effRange = sqrt(8*1)/kappa
  print(paste("Effective Range: ", round(effRange, 2), "Km"))
  
  marVar = 1.0 / (4.0*3.14159265359*(kappa^2)*(tau^2))
  print(paste("Marginal Variance: ", round(marVar, 2)))
  
  M = dim(X)[2] # number of earthquakes
  B = dim(X)[1] # number of basis functions
  
  # Check that X is different between earthquakes
  plotDF = data.frame(Quake=as.factor(rep(1:M, each=B)),
                      x = c(X),
                      triangle = rep(1:B, M))
  g = ggplot(plotDF, aes(x=triangle, y=x, colour=Quake))+
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
  
  # Calculate the taper
  taper = exp(-depths*lambda)
  
  # Plot the taper
  g = plotFault(fault, z=taper, legendTitle="Taper")
  plot(g)
  rm(g)
  
  # Plot the shared effect W
  g = plotX(mesh, z=w, legendTitle="Spatial Effect W")
  plot(g)
  rm(g)
  
  if (plotAll){
    # Check the underlying distributions
    earthquakes = c("T1","T2","T3","T4","T5","T6","T7","T8","T9","T10","T11","T12")
    
    for (i in 1:M){
      thisX = X[,i]
      g = plotX(mesh, z=thisX, legendTitle="Spatial Effect X")
      plotDF = data.frame(x=DR$Lon[DR$event == earthquakes[i]],
                          y=DR$Lat[DR$event == earthquakes[i]])
      g = g +
        geom_point(data=plotDF, aes(x=x, y=y), colour="red") +
        ggtitle(paste("Earthquake ", i))
      
      plot(g)
      rm(g)
      
      thisSlip = taper * exp(mu + sqrt(phi)*w[inla_mesh$idx$loc] + sqrt(1-phi)*thisX[inla_mesh$idx$loc])
      # plot tapered slips
      g = plotFault(fault, z=thisSlip, legendTitle="Tapered Slips (m)")
      g = g +
        ggtitle(paste("Earthquake ", i))
      plot(g)
      rm(g)
    }
  } else{
    # Just plot the 1700 quake
    
    # extract 1700 quake
    thisX = X[,1]
    g = plotX(mesh, z=thisX, legendTitle="Spatial Effect X")
    plotDF = data.frame(x=DR$Lon[DR$event == earthquakes[1]],
                        y=DR$Lat[DR$event == earthquakes[1]])
    g = g +
      geom_point(data=plotDF, aes(x=x, y=y), colour="red") +
      ggtitle("1700 Earthquake")
    plot(g)
    rm(g)
    
    thisSlip = taper * exp(mu + sqrt(phi)*w[inla_mesh$idx$loc] + sqrt(1-phi)*thisX[inla_mesh$idx$loc])
    # plot tapered slips
    g = plotFault(fault, z=thisSlip, legendTitle="Tapered Slips (m)")
    g = g +
      ggtitle("1700 Earthquake")
    plot(g)
    rm(g)
  }
}

# Check the output of the multiQuake earthquake
# The different tau and kappa models
checkOptimisationSharedAniso = function(obj, opt0, allPar, mesh, fault, DR, plotAll=F){
  print("Last Parameters:")
  print(tail(obj$env$last.par, 8))
  print("Best Parameters")
  print(tail(obj$env$last.par.best, 8))
  print(opt0)
  
  if (all(tail(obj$env$last.par, 8) == opt0$par)){
    thispar = obj$env$last.par
    print("Last parameters")
  } else if (all(tail(obj$env$last.par.best, 8) == opt0$par)){
    thispar = obj$env$last.par.best
    print("Not last parameters")
  } else{
    stop("Parameters aren't last params")
  }
  
  # some nice colours for plotting
  CB_color_cycle = c('#377eb8', '#ff7f00', '#4daf4a',
                     '#f781bf', '#a65628', '#984ea3',
                     '#999999', '#e41a1c', '#dede00')
  
  # plot the iterations of parameters
  I = length(allPar[,1])
  for (j in 1:8){
    plot(x = c(1:I), y= allPar[,j],
         type='o',
         col = CB_color_cycle[j],
         xlab = "Iteration",
         ylab = colnames(allPar)[j])
  }
  
  # get the report from the TMB object
  r = obj$report(thispar)
  
  lambda = exp(opt0$par[1])
  mu = opt0$par[2]
  kappaX = exp(opt0$par[3])
  tauX =  exp(opt0$par[4])
  kappaW = exp(opt0$par[5])
  tauW =  exp(opt0$par[6])
  h1 = opt0$par[7]
  h2 = opt0$par[8]
  
  X = r$X
  w = r$w
  
  # Find the effective range
  effRangeX = sqrt(8*1)/kappaX
  print(paste("Effective Range X: ", round(effRangeX, 2), "Km"))
  
  marVarX = 1.0 / (4.0*3.14159265359*(kappaX^2)*(tauX^2))
  print(paste("Marginal Variance X: ", round(marVarX, 2)))
  
  effRangeW = sqrt(8*1)/kappaW
  print(paste("Effective Range W: ", round(effRangeW, 2), "Km"))
  
  marVarW = 1.0 / (4.0*3.14159265359*(kappaW^2)*(tauW^2))
  print(paste("Marginal Variance W: ", round(marVarW, 2)))
  
  # Check the anisotropy
  h = c(h1, h2)
  print(h)
  anisoDirection = 90 + rad2deg(atan(-h2/h1))
  print(anisoDirection)
  anisoStrength1 = sqrt(1 + h1^2 + h2^2)
  print(anisoStrength1)
  anisoStrength2 = 1 / sqrt(1 + h1^2 + h2^2)
  print(anisoStrength2)
  
  
  M = dim(X)[2] # number of earthquakes
  B = dim(X)[1] # number of basis functions
  
  # Check that X is different between earthquakes
  plotDF = data.frame(Quake=as.factor(rep(1:M, each=B)),
                      x = c(X),
                      triangle = rep(1:B, M))
  g = ggplot(plotDF, aes(x=triangle, y=x, colour=Quake))+
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
  
  # Calculate the taper
  taper = exp(-depths*lambda)
  
  # Plot the taper
  g = plotFault(fault, z=taper, legendTitle="Taper")
  plot(g)
  rm(g)
  
  # Plot the shared effect W
  g = plotX(mesh, z=w, legendTitle="Spatial Effect W")
  plot(g)
  rm(g)
  
  if (plotAll){
    # Check the underlying distributions
    earthquakes = c("T1","T2","T3","T4","T5","T6","T7","T8","T9","T10","T11","T12")
    
    for (i in 1:M){
      thisX = X[,i]
      g = plotX(mesh, z=thisX, legendTitle="Spatial Effect X")
      plotDF = data.frame(x=DR$Lon[DR$event == earthquakes[i]],
                          y=DR$Lat[DR$event == earthquakes[i]])
      g = g +
        geom_point(data=plotDF, aes(x=x, y=y), colour="red") +
        ggtitle(paste("Earthquake ", i))
      
      plot(g)
      rm(g)
      
      thisSlip = taper * exp(mu + w[inla_mesh$idx$loc] + thisX[inla_mesh$idx$loc])
      # plot tapered slips
      g = plotFault(fault, z=thisSlip, legendTitle="Tapered Slips (m)")
      g = g +
        ggtitle(paste("Earthquake ", i))
      plot(g)
      rm(g)
    }
  } else{
    # Just plot the 1700 quake
    
    # extract 1700 quake
    thisX = X[,1]
    g = plotX(mesh, z=thisX, legendTitle="Spatial Effect X")
    plotDF = data.frame(x=DR$Lon[DR$event == earthquakes[1]],
                        y=DR$Lat[DR$event == earthquakes[1]])
    g = g +
      geom_point(data=plotDF, aes(x=x, y=y), colour="red") +
      ggtitle("1700 Earthquake")
    plot(g)
    rm(g)
    
    thisSlip = taper * exp(mu + w[inla_mesh$idx$loc] + thisX[inla_mesh$idx$loc])
    # plot tapered slips
    g = plotFault(fault, z=thisSlip, legendTitle="Tapered Slips (m)")
    g = g +
      ggtitle("1700 Earthquake")
    plot(g)
    rm(g)
  }
}


summariseRandomEffects = function(xDraws, mesh, quake="1"){
  # Visualise the x draws and their standard deviations
  xSum = cbind(mean   = (apply(xDraws, 1, mean)),
               median = (apply(xDraws, 1, median)),
               sd     = (apply(xDraws, 1, sd)))
  
  # plot the mean x distribution
  g5 = plotX(mesh, z=xSum[,1], legendTitle="Spatial Random Effect\nMean")
  #g5 = g5 + ggtitle(paste("Earthquake ", quake))
  
  # plot the standard deviation of x distribution
  g6 = plotX(mesh, z=xSum[,3], legendTitle="Spatial Random Effect\nStandard Deviation", colourScale="plasma")
  #g6 = g6 + ggtitle(paste("Earthquake ", quake))
  
  plot(g5)
  plot(g6)
}
