#---------------------------MultiQuakeAnisoCSZ----------------------------------

setwd("~/Uni/NTNU/Masters Project/CSZ/R/Code/MultiQuakeAnisoCSZ/Results")
load("Draws.RData")
load("Times.RData")
load("SD0.RData")
load("Best Params.RData")

setwd("~/Uni/NTNU/Masters Project/CSZ/R/Code/MultiQuakeAnisoCSZ")
load("Fault.RData")
load("spdeMesh.RData")
load("Data.RData")
load("DR4.RData")

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
                        Kappa  = exp(results$logKappaDraws),
                        Tau    = exp(results$logTauDraws),
                        Psi    = exp(results$logPsiDraws))

# Now summarise the useful parameters
fixedDraws$rhoStrike = fixedDraws$Psi * (sqrt(8)/fixedDraws$Kappa)
fixedDraws$rhoDip    = (1 / fixedDraws$Psi) * (sqrt(8)/fixedDraws$Kappa)
fixedDraws$sigma2 = 1 / (4*pi*(fixedDraws$Kappa^2)*(fixedDraws$Tau^2))

# Calculate mean, median, and standard deviation for each column
fixedSummary = cbind(mean   = (apply(fixedDraws, 2, mean)),
                     median = (apply(fixedDraws, 2, median)),
                     sd     = (apply(fixedDraws, 2, sd)),
                     lower  = (apply(fixedDraws, 2, quantile, .05)),
                     upper  = (apply(fixedDraws, 2, quantile, .95)))
print(fixedSummary)

usefulDraws = data.frame(Lambda    = fixedDraws$Lambda,
                         Mu        = fixedDraws$Mu,
                         Psi       = fixedDraws$Psi,
                         rhoStrike = fixedDraws$rhoStrike,
                         rhoDip    = fixedDraws$rhoDip,
                         sigma2    = fixedDraws$sigma2)

# call a function that gives histograms and basic statistics of fixed pars
histLabel = c('\u03bb', '\u03bc', '\u03c6', '\u03C1 - Strike', '\u03C1 - Dip', '\u03c3^2')
plotFixedParameters(usefulDraws, histLabel=histLabel)

# Now calculate the prediction bounds for the taper
plotTaper(usefulDraws$Lambda)

# Now show the anisotropic factor range
plotCorrelation(usefulDraws$Psi)

#---------------------------Spatial Random Effects------------------------------

# Then all the Xs together
plotAllX(spdeMesh, results$xDraws, DR4)

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
    theseSlips = taper * exp(results$muDraws[n] + (data$A %*% thisXDraws[,n]))
    theseSlipDraws[,n] = theseSlips[,1]
  }
  
  # save slip draws
  slipDraws[[j]] = theseSlipDraws
}

setwd("~/Uni/NTNU/Masters Project/CSZ/R/Code/MultiQuakeAnisoCSZ/Results")
save(slipDraws, file="Slips.RData")

# Now plot all the slip distributions
load("Slips.RData")
plotAllSlips(fault, slipDraws, DR4)

#---------------------------Magnitudes------------------------------------------

# Now get the magnitudes from all slips
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
okadaSubs = list()
for (j in 1:E){
  thisG = G[[j]]
  theseSlipDraws = slipDraws[[j]]
  
  okadaSubSims = matrix(0, nrow=dim(thisG)[1], ncol=nSims)
  for (n in 1:nSims){
    thisSlip = theseSlipDraws[,n]
    okadaSubSims[,n] = thisG %*% thisSlip
  }
  
  okadaSubs[[j]] = -okadaSubSims
}
setwd("~/Uni/NTNU/Masters Project/CSZ/R/Code/MultiQuakeAnisoCSZ/Results")
save(okadaSubs, file="Subsidence.RData")

load("Subsidence.RData")
# plot predicted subsidence vs data all
plotAllSubsidencePredicition(okadaSubs, DR4)

# Now predict across the whole fault
lons = seq(-128, -123, length.out=25)
lats = seq(40, 50, length.out=50)
lonLat = expand.grid(lons, lats)
lonLat = cbind(as.numeric(lonLat[,1]),
               as.numeric(lonLat[,2]))
setwd("~/Uni/NTNU/Masters Project/CSZ/R/Code/MultiQuakeAnisoCSZ")
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

setwd("~/Uni/NTNU/Masters Project/CSZ/R/Code/MultiQuakeAnisoCSZ/Results")
save(wholeSubsidence, file="Whole Subsidence.RData")

load("Whole Subsidence.RData")
plotAllSubsidenceGrid(wholeSubsidence, lonLat, DR4)

#---------------------------Future Quake----------------------------------------

# Predict a future earthquake
lambda = exp(median(results$logLambdaDraws))
mu     = median(results$muDraws)
kappa  = exp(median(results$logKappaDraws))
tau    = exp(median(results$logTauDraws))
psi    = exp(median(results$logPsiDraws))

# Add the TMB code to make the Q matrix
setwd("~/Uni/NTNU/Masters Project/CSZ/R/Code/MultiQuakeAnisoCSZ")
compile("MakeQAniso.cpp",
        framework="TMBad")
dyn.load(dynlib("MakeQAniso"))

# Create objects needed
qData       = list(spde = data$spde)
qParameters = list(logKappa = log(kappa),
                   logTau   = log(tau),
                   logh     = c(log(psi), deg2rad(-12.39776)))
# Make TMB object
obj = MakeADFun(qData,
                qParameters,
                DLL     = "MakeQAniso")
# extract Q matrix
rep = obj$report()
Qh = as.matrix(rep$Q)

futureSlips = matrix(data=0, nrow=S, ncol=nSims)
taper = exp(-lambda*data$depth)
allX = inla.qsample(n=nSims, Q=Qh)
for (n in 1:nSims){
  X = allX[,n]
  theseSlips = taper * exp(mu + (data$A %*% X))
  futureSlips[,n] = theseSlips[,1]
}

# save data
setwd("~/Uni/NTNU/Masters Project/CSZ/R/Code/MultiQuakeAnisoCSZ/Results")
save(futureSlips, file="Future Slips.RData")

rm(futureSlips)
load("Future Slips.RData")
# plot the future mean slip
plotFutureSlip(fault, futureSlips)


# Magnitude of future quake
magnitudes = rep(0, nSims)
for (n in 1:nSims){
  mag = getMomentFromSlip(futureSlips[,n], fault)
  magnitudes[n] = mag
}
# plot the future magnitude distribution
plotOneMagnitude(magnitudes)

setwd("~/Uni/NTNU/Masters Project/CSZ/R/Code/MultiQuakeAnisoCSZ")
load("Gwhole.RData")
futureSubs = matrix(0, nrow=dim(GWhole)[1], ncol=nSims)
for (n in 1:nSims){
  futureSubs[,n] = -(GWhole %*% futureSlips[,n])
}
# save data
setwd("~/Uni/NTNU/Masters Project/CSZ/R/Code/MultiQuakeAnisoCSZ/Results")
save(futureSubs, file="Future Subsidence.RData")

rm(futureSubs)
load("Future Subsidence.RData")
# plot results
plotFutureSubsidenceGrid(futureSubs, lonLat)

