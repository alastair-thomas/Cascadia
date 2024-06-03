# load data needed
setwd("~/Uni/NTNU/Masters Project/CSZ/R/Code/MultiQuakeCSZ")
load("Fault.RData")
load("spdeMesh.RData")
load("Data.RData")
load("DR4.RData")

# load simulated draws and other results
setwd("~/Uni/NTNU/Masters Project/CSZ/R/Code/MultiQuakeCSZ/Results")
load("Draws.RData")
load("Times.RData")
load("SD0.RData")
load("Pars.RData")

results = allDraws
rm(allDraws)

# set constants needed
B     = dim(spdeMesh$loc)[1]    # number of basis functions
E     = length(data$subsidence) # number of earthquakes
nSims = length(results$muDraws) # number of simulations (1000)
S     = length(fault)           # number of subfaults

#-----------------------------Fixed Parameters----------------------------------

# A dataframe for all the fixed parameters
fixedDraws = data.frame(Lambda = exp(results$logLambdaDraws),
                        Mu     = results$muDraws,
                        Kappa  = exp(results$logKappaDraws),
                        Tau    = exp(results$logTauDraws))

# Now summarise the useful parameters
fixedDraws$rho = sqrt(8)/fixedDraws$Kappa
fixedDraws$sigma2 = 1 / (4*pi*(fixedDraws$Kappa^2)*(fixedDraws$Tau^2))


usefulDraws = data.frame(Lambda = fixedDraws$Lambda,
                         Mu     = fixedDraws$Mu,
                         rho    = fixedDraws$rho,
                         sigma2 = fixedDraws$sigma2)

# call a function that gives histograms and basic statistics of fixed pars
histLabel = c('\u03bb', '\u03bc', '\u03C1', '\u03c3^2')
plotFixedParameters(usefulDraws, histLabel=histLabel)

# Now calculate the prediction bounds for the taper
plotTaper(usefulDraws$Lambda)

#-----------------------------Random Spatial Effects----------------------------

# Plot the mean and sd of all the X distributions
plotAllX(spdeMesh, results$xDraws, DR4)

#-----------------------------Slip Distributions--------------------------------

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

setwd("~/Uni/NTNU/Masters Project/CSZ/R/Code/MultiQuakeCSZ/Results")
save(slipDraws, file="Slips.RData")

# Now plot all the slip distributions
rm(slipDraws)
load("Slips.RData")
plotAllSlips(fault, slipDraws, DR4)

#-----------------------------Magnitudes----------------------------------------

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

#-----------------------------Subsidence----------------------------------------

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
setwd("~/Uni/NTNU/Masters Project/CSZ/R/Code/MultiQuakeCSZ/Results")
save(okadaSubs, file="Subsidence.RData")

load("Subsidence.RData")
# plot predicted subsidence vs data all
plotAllSubsidencePredicition(okadaSubs, DR4)

# Now plot over the entire fault
lons = seq(-128, -123, length.out=25)
lats = seq(40, 50, length.out=50)
lonLat = expand.grid(lons, lats)
lonLat = cbind(as.numeric(lonLat[,1]),
               as.numeric(lonLat[,2]))

# plot the grid
gridPlot = plotBase(proj="UTM", labels=F, countryBoundary=F)
xy = projCSZ(lonLat, units="m")
plotDF = data.frame(x=xy[,1], y=xy[,2])
gridPlot = gridPlot +
  geom_point(data=plotDF, aes(x=x, y=y), size=1)
limits = data.frame(x = -c(128, 123),
                    y = c(40, 50))# Fix the limits
limits = st_as_sf(x=limits,
                  coords = c("x", "y"),
                  crs=st_crs("EPSG:4326"))
limits = toUTM(limits)
gridPlot = gridPlot + 
  coord_sf(xlim = c(st_coordinates(limits)[,1]),
           ylim = c(st_coordinates(limits)[,2]),
           crs=st_crs("EPSG:32610"))
plot(gridPlot)


#GWhole = getOkada(fault, lon=lonLat[,1], lat=lonLat[,2])

setwd("~/Uni/NTNU/Masters Project/CSZ/R/Code/MultiQuakeCSZ")
load("GWhole.RData")

wholeSubsidence = list()
for (j in 1:E){
  theseSlipDraws = slipDraws[[j]]
  
  Subs = matrix(0, nrow=dim(GWhole)[1], ncol=nSims)
  for (n in 1:nSims){
    thisSlip = theseSlipDraws[,n]
    Subs[,n] = -(GWhole %*% thisSlip)
  }
  wholeSubsidence[[j]] = Subs
}

setwd("~/Uni/NTNU/Masters Project/CSZ/R/Code/MultiQuakeCSZ/Results")
save(wholeSubsidence, file="Whole Subsidence.RData")

load("Whole Subsidence.RData")
plotAllSubsidenceGrid(wholeSubsidence, lonLat, DR4)
#-----------------------------Future Quake--------------------------------------

# Predict a future earthquake
lambda = exp(median(results$logLambdaDraws))
mu     = median(results$muDraws)
kappa  = exp(median(results$logKappaDraws))
tau    = exp(median(results$logTauDraws))

# get Q from TMB code
# Add the TMB code to make the Q matrix
setwd("~/Uni/NTNU/Masters Project/CSZ/R/Code/MultiQuakeCSZ")
load("Data.RData")

compile("MakeQ.cpp")
dyn.load(dynlib("MakeQ"))

# Create objects needed
qData = list(M0 = data$M0,
             M1 = data$M1,
             M2 = data$M2)
qParameters = list(logKappa = log(kappa),
                   logTau   = log(tau))
# Make TMB object
obj = MakeADFun(qData,
                qParameters,
                DLL = "MakeQ")
# extract Q matrix
rep = obj$report()
Qx = as.matrix(rep$Q)
dyn.unload(dynlib("MakeQ"))

futureSlips = matrix(data=0, nrow=S, ncol=nSims)
taper = exp(-lambda*data$depth)

# simulate from X distribution
allX = inla.qsample(n=nSims, Q=Qx)
for (n in 1:nSims){
  X = allX[,n]
  theseSlips = taper * exp(mu + (data$A %*% X))
  futureSlips[,n] = theseSlips[,1]
}
# save data
setwd("~/Uni/NTNU/Masters Project/CSZ/R/Code/MultiQuakeCSZ/Results")
save(futureSlips, file="Future Slips.RData")

# plot the future mean slip
load("Future Slips.RData")
plotFutureSlip(fault, futureSlips)


# Magnitude of future quake
magnitudes = rep(0, nSims)
for (n in 1:nSims){
  mag = getMomentFromSlip(futureSlips[,n], fault)
  magnitudes[n] = mag
}
# plot the future magnitude distribution
plotOneMagnitude(magnitudes)



# Predict the subsidence of a future earthquake
# Plot the subsidence across the entire fault
setwd("~/Uni/NTNU/Masters Project/CSZ/R/Code/MultiQuakeCSZ")
load("GWhole.RData")
futureSubs = matrix(0, nrow=dim(GWhole)[1], ncol=nSims)
for (n in 1:nSims){
  futureSubs[,n] = -(GWhole %*% futureSlips[,n])
}
# save data
setwd("~/Uni/NTNU/Masters Project/CSZ/R/Code/MultiQuakeCSZ/Results")
save(futureSubs, file="Future Subsidence.RData")

# plot results
load("Future Subsidence.RData")
plotFutureSubsidenceGrid(futureSubs, lonLat)
