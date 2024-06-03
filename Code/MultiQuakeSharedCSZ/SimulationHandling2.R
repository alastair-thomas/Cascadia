# Model results
setwd("~/Uni/NTNU/Masters Project/CSZ/R/Code/MultiQuakeSharedCSZ/Results")
load("Draws.RData")
load("Times.RData")
load("SD0.RData")
load("Best Params.RData")

# model data
setwd("~/Uni/NTNU/Masters Project/CSZ/R/Code/MultiQuakeSharedCSZ")
load("Fault.RData")
load("spdeMesh.RData")
load("Data.RData")
load("DR4.RData")

# model constants
S = length(fault)
B = dim(spdeMesh$loc)[1]
E = length(data$subsidence)
nSims = 1000

#-----------------------------Direct Outputs----------------------------------
print(bestPars)
print(optTime)
print(sdTime)
print(SD0)

#-----------------------------Fixed Parameters----------------------------------

# A dataframe for all the fixed parameters
fixedDraws = data.frame(Lambda = exp(results$logLambdaDraws),
                        Mu     = results$muDraws,
                        KappaX = exp(results$logKappaXDraws),
                        TauX   = exp(results$logTauXDraws),
                        KappaW = exp(results$logKappaWDraws),
                        TauW   = exp(results$logTauWDraws))

# Calculate the useful re-parameterisations of the fixed effects
rhoDrawsX    = sqrt(8) / fixedDraws$KappaX
sigma2DrawsX = 1.0 / (4*pi*(fixedDraws$KappaX^2)*(fixedDraws$TauX)^2)
rhoDrawsW    = sqrt(8) / fixedDraws$KappaW
sigma2DrawsW = 1.0 / (4*pi*(fixedDraws$KappaW^2)*(fixedDraws$TauW)^2)

# add to dataframe
fixedDraws$rhoX    = rhoDrawsX
fixedDraws$sigma2X = sigma2DrawsX
fixedDraws$rhoW    = rhoDrawsW
fixedDraws$sigma2W = sigma2DrawsW

# Summarise the fixed parameters
fixedSummary = cbind(mean   = (apply(fixedDraws, 2, mean)),
                     median = (apply(fixedDraws, 2, median)),
                     sd     = (apply(fixedDraws, 2, sd)),
                     lower  = (apply(fixedDraws, 2, quantile, .05)),
                     upper  = (apply(fixedDraws, 2, quantile, .95)))
print(round(fixedSummary, 6))

usefulDraws = data.frame(Lambda  = exp(results$logLambdaDraws),
                         Mu      = results$muDraws,
                         RhoX    = rhoDrawsX,
                         Sigma2X = sigma2DrawsX,
                         RhoW    = rhoDrawsW,
                         Sigma2W = sigma2DrawsW)

# call a function that gives histograms and basic statistics of fixed pars
plotFixedParameters(usefulDraws,
                    histLabel=c('\u03bb', '\u03bc',
                                '\u03c1 - Independent', '\u03c3^2 - Independent',
                                '\u03c1 - Shared', '\u03c3^2 - Shared'))

# Now calculate the prediction bounds for the taper
plotTaper(exp(results$logLambdaDraws))

#---------------------------Spatial Random Effects------------------------------


# Then all the Xs together
plotAllX(spdeMesh, results$xDraws, DR4)

# The the w field
plotOneX(spdeMesh, as.matrix(results$wDraws), DR4, xw="W")

#---------------------------Slip Distribution-----------------------------------

# Find the slip summary for the simulated data
# Create the slip distributions
slipDraws = list()
earthquakes = as.factor(c("T1","T2","T3","T4","T5","T6",
                          "T7","T8","T9","T10","T11","T12"))
for (j in 1:E){
  thisXDraws = as.matrix(results$xDraws[[j]])
  ## create the slip draws
  theseSlipDraws = matrix(data=0, nrow=S, ncol=nSims)
  for (n in 1:nSims){
    taper = exp(-exp(results$logLambdaDraws[n])*data$depth)
    theseSlips = taper * exp(results$muDraws[n] + (data$A %*% (thisXDraws[,n] + results$wDraws[,n])))
    theseSlipDraws[,n] = theseSlips[,1]
  }
  
  # save slip draws
  slipDraws[[j]] = theseSlipDraws
}

setwd("~/Uni/NTNU/Masters Project/CSZ/R/Code/MultiQuakeSharedCSZ/Results")
save(slipDraws, file="Slips.RData")

# Here I need to calculate the limits for plotting

# Now plot all the slip distributions
load("Slips.RData")
plotAllSlips(fault, slipDraws, DR4)

#---------------------------Magnitudes------------------------------------------

# Now get the magnitudes from all slips
nSims = 1000
magnitudes = matrix(0, nrow=nSims, ncol=E)
for (j in 1:E){
  for (n in 1:nSims){
    mag = getMomentFromSlip(slipDraws[[j]][,n], fault)
    magnitudes[n,j] = mag
  }
}

# plot all of them
plotAllMagnitudes(magnitudes)

#---------------------------Subsidences-----------------------------------------

# Now predict subsidence and plot against real data
# save the subsidence estimate for later model comparisons
G = data$okada
subsidences = list()
for (e in 1:E){
  thisG = G[[e]]
  theseSlipDraws = slipDraws[[e]]
  
  okadaSubSims = matrix(0, nrow=dim(thisG)[1], ncol=nSims)
  for (n in 1:nSims){
    thisSlip = theseSlipDraws[,n]
    okadaSubSims[,n] = thisG %*% thisSlip
  }
  
  subsidences[[e]] = -okadaSubSims
}
setwd("~/Uni/NTNU/Masters Project/CSZ/R/Code/MultiQuakeSharedCSZ/Results")
save(subsidences, file="Subsidence.RData")

load("Subsidence.RData")
# plot predicted subsidence vs data all
plotAllSubsidencePredicition(subsidences, DR4)

# Now plot over the entire fault
lons = seq(-128, -123, length.out=25)
lats = seq(40, 50, length.out=50)
lonLat = expand.grid(lons, lats)
lonLat = cbind(as.numeric(lonLat[,1]),
               as.numeric(lonLat[,2]))
# GWhole = getOkada(fault, lon=lonLat[,1], lat=lonLat[,2])
setwd("~/Uni/NTNU/Masters Project/CSZ/R/Code/MultiQuakeSharedCSZ")
load("GWhole.RData")

wholeSubsidence = list()
for (j in 1:E){
  theseSlipDraws = slipDraws[[j]]
  
  thisWholeSub = matrix(0, nrow=dim(GWhole)[1], ncol=nSims)
  
  for (n in 1:nSims){
    thisSlip = theseSlipDraws[,n]
    thisWholeSub[,n] = -(GWhole %*% thisSlip)
  }
  wholeSubsidence[[j]] = thisWholeSub
}

setwd("~/Uni/NTNU/Masters Project/CSZ/R/Code/MultiQuakeSharedCSZ/Results")
save(wholeSubsidence, file="Whole Subsidence.RData")

load("Whole Subsidence.RData")
plotAllSubsidenceGrid(wholeSubsidence, lonLat, DR4)

#---------------------------Future Quake----------------------------------------

# Predict a future earthquake
lambda = exp(median(results$logLambdaDraws))
mu     = median(results$muDraws)
kappaX = exp(median(results$logKappaXDraws))
tauX   = exp(median(results$logTauXDraws))
kappaW = exp(median(results$logKappaWDraws))
tauW   = exp(median(results$logTauWDraws))

# Add the TMB code to make the Q matrix
setwd("~/Uni/NTNU/Masters Project/CSZ/R/Code/MultiQuakeSharedCSZ")
compile("MakeQShared.cpp")
dyn.load(dynlib("MakeQShared"))

# Create objects needed
qData       = list(M0 = data$M0,
                   M1 = data$M1,
                   M2 = data$M2)
qParameters = list(logKappaX = log(kappaX),
                   logTauX   = log(tauX),
                   logKappaW = log(kappaW),
                   logTauW   = log(tauW))
# Make TMB object
obj = MakeADFun(qData,
                qParameters,
                DLL = "MakeQShared")
# extract Q matrix
rep = obj$report()
Qx = as.matrix(rep$Qx)

futureSlips = matrix(data=0, nrow=S, ncol=nSims)
taper = exp(-lambda*data$depth)

allX = inla.qsample(n=nSims, Q=Qx)
allW = results$wDraws
for (n in 1:nSims){
  X = allX[,n]
  w = allW[,n]
  theseSlips = taper * exp(mu + (data$A %*% (X+w)))
  futureSlips[,n] = theseSlips[,1]
}
setwd("~/Uni/NTNU/Masters Project/CSZ/R/Code/MultiQuakeSharedCSZ/Results")
save(futureSlips, file="Future Slips.RData")

rm(futureSlips)
load("Future Slips.RData")
# plot the future mean slip
plotFutureSlip(fault, futureSlips)


# Magnitude of future quake
futureMagnitude = rep(0, nSims)
for (n in 1:nSims){
  mag = getMomentFromSlip(futureSlips[,n], fault)
  futureMagnitude[n] = mag
}
# plot the future magnitude distribution
plotOneMagnitude(futureMagnitude)

setwd("~/Uni/NTNU/Masters Project/CSZ/R/Code/MultiQuakeSharedCSZ")
load("GWhole.RData")
# Now plot the subsidence across the whole fault
futureSubs = matrix(0, nrow=dim(GWhole)[1], ncol=nSims)
for (n in 1:nSims){
  futureSubs[,n] = -(GWhole %*% futureSlips[,n])
}

setwd("~/Uni/NTNU/Masters Project/CSZ/R/Code/MultiQuakeSharedCSZ/Results")
save(futureSubs, file="Future Subsidence.RData")

# now plot it.
load("Future Subsidence.RData")
plotFutureSubsidenceGrid(futureSubs, lonLat)
