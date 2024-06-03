# load data needed
setwd("~/Uni/NTNU/Masters Project/CSZ/R/Code/BayesianCSZ")
load("fault.RData")
load("spdeMesh.RData")
load("BayesianData.RData")
load("DR4.RData")

# load simulated draws and other results
setwd("~/Uni/NTNU/Masters Project/CSZ/R/Code/BayesianCSZ")
load("1700results.RData")

# set constants needed
B = dim(spdeMesh$loc)[1]        # number of basis functions
E = length(data$subsidence)     # number of earthquakes
nSims = length(results$muDraws) # number of simulations (1000)
S = length(fault)               # number of subfaults

# A dataframe for all the fixed parameters
fixedDraws = data.frame(Lambda = exp(results$logLambdaDraws),
                        Mu     = results$muDraws,
                        Kappa  = exp(results$logKappaDraws),
                        Tau    = exp(results$logTauDraws))

# call a function that gives histograms and basic statistics of fixed pars
# remember this saves a plot so need to change file name
summariseFixedParameters(fixedDraws, histLabel=c('\u03bb', '\u03bc', '\u03ba', '\u03c4'))

# Now calculate the prediction bounds for the taper
depths = seq(0,30,0.01)
lambdaQ = stats::quantile(exp(results$logLambdaDraws), probs=c(0.025, 0.5, 0.975))
tlower = exp(-lambdaQ[1]*depths)
tmedian = exp(-lambdaQ[2]*depths)
tupper = exp(-lambdaQ[3]*depths)
plotDF = data.frame(depth = depths,
                    lower = tlower,
                    taper = tmedian,
                    upper = tupper)
ggplot(plotDF) +
  geom_ribbon(aes(x=depth, ymin=lower, ymax=upper), fill="lightblue", alpha=0.75) +
  geom_line(aes(x=depth, y=taper), colour="black") +
  labs(x="Depth (km)", y="Taper")

# Calculate the useful re-parameterisations of the fixed effects
rhoDraws    = sqrt(8) / fixedDraws$Kappa
sigma2Draws = 1.0 / (4*pi*(fixedDraws$Kappa^2)*(fixedDraws$Tau)^2)

usefulParams = cbind(Rho   = rhoDraws,
                     Sigma = sigma2Draws)

usefulSummary = cbind(mean   = apply(usefulParams, 2, mean),
                      median = apply(usefulParams, 2, median),
                      sd     = apply(usefulParams, 2, sd))
print(usefulSummary)


# Draws plots of the spatial effect and standard deviation
# First the 1700 Quake by itself
plotOneX(spdeMesh, as.matrix(results$xDraws), DR4)

# Create the slip distributions
earthquakes = as.factor(c("T1","T2","T3","T4","T5","T6","T7","T8","T9","T10","T11","T12"))
## create the slip draws
slipDraws = matrix(data=0, nrow=S, ncol=nSims)
for (n in 1:nSims){
  taper = exp(-exp(results$logLambdaDraws[n])*data$depth)
  theseSlips = taper * exp(results$muDraws[n] + (data$A %*% results$xDraws[,n]))
  slipDraws[,n] = theseSlips[,1]
}


# Now plot 1700 quake
mean1700  = apply(slipDraws, 1, mean)
sd1700    = apply(slipDraws, 1, sd)
lower1700 = apply(slipDraws, 1, quantile, 0.025)
upper1700 = apply(slipDraws, 1, quantile, 0.975)
width1700 = upper1700 - lower1700

g = plotFault(fault, z=mean1700, legendTitle="Posterior Slip\nMean (m)")
plot(g)
g = plotFault(fault, z=sd1700, legendTitle="Posterior Slip\nStandard Deviation (m)", colourScale="plasma")
plot(g)
g = plotFault(fault, z=width1700, legendTitle="Posterior Slip\n95% PI Width (m)", colourScale="mako")
plot(g)

# Plot a few realisations
rs = c(348, 907, 172, 683)
for (r in rs){
  g = plotFault(fault, z=slipDraws[,r], legendTitle=paste("Realisation", r, "\nSlip (m)"))
  plot(g)
}


# Now get the magnitudes from all slips
nSims = 1000
magnitudes = rep(0, nSims)
for (n in 1:nSims){
  mag = getMomentFromSlip(slipDraws[,n], fault)
  magnitudes[n] = mag
}

# plot the 1700 magnitude
plotOneMagnitude(magnitudes)

# Now predict subsidence and plot against real data
# save the subsidence estimate for later model comparisons
G = data$okada
okadaSubSims = matrix(0, nrow=dim(G)[1], ncol=nSims)
for (n in 1:nSims){
  thisSlip = slipDraws[,n]
  okadaSubSims[,n] = -G %*% thisSlip
}
setwd("~/Uni/NTNU/Masters Project/CSZ/R/Code/BayesianCSZ")
save(okadaSubSims, file="Subsidence.RData")

# plot predicted subsidence vs data 1700
# colour by error
# type by sign error
plotOneSubsidencePrediction(okadaSubSims, DR4)

# Now predict subsidence across the whole fault
lons = seq(-128, -123, length.out=10)
lats = seq(40, 50, length.out=20)
lonLat = expand.grid(lons, lats)
Gwhole = getOkada(fault, lonLat[,1], lonLat[,2])
allSubs = matrix(0, nrow=length(lons)*length(lats), ncol=nSims)
for (n in 1:nSims){
  allSubs[,n] = -Gwhole %*% slipDraws[,n]
}
plotSubsidenceGrid3(allSubs, cbind(lonLat[,1],
                                   lonLat[,2]), DR4)




# Predict a future earthquake
setwd("~/Uni/NTNU/Masters Project/CSZ/R/Code/MultiQuakeCSZ/Backup Results")
load("Backup Pars.RData")
print(bestPars)
lambda = exp(bestPars[1])
mu  = bestPars[2]
kappa = exp(bestPars[3])
tau = exp(bestPars[4])
nSims = 1000

# get Q from TMB code
# Add the TMB code to make the Q matrix
setwd("~/Uni/NTNU/Masters Project/CSZ/R/Code/MultiQuakeCSZ")
load("backupData.RData")

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
                DLL     = "MakeQ")
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

## find the median and sd across draws, as well as 95% intervals
futureSlipSum = cbind(mean   = (apply(futureSlips, 1, mean)),
                      sd     = (apply(futureSlips, 1, sd)),
                      width  = (apply(futureSlips, 1, quantile, .975)) -
                        (apply(futureSlips, 1, quantile, .025)))
# plot the future mean slip
g = plotFault(fault,
              z=futureSlipSum[,1],
              legendTitle="Future Slip\nMean (m)")
plot(g)

# plot the future slip standard deviation
g = plotFault(fault,
              z=futureSlipSum[,2],
              legendTitle="Future Slip\nStandard Deviation (m)",
              colourScale="plasma")
plot(g)

# plot the future 95% width
g = plotFault(fault,
              z=futureSlipSum[,3],
              legendTitle="Future Slip\n95% PI Width (m)",
              colourScale="mako")
plot(g)


# Magnitude of future quake
nSims = 1000
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
load("G Whole Backup.RData")
futureSubs = matrix(0, nrow=dim(GFault)[1], ncol=nSims)
for (n in 1:nSims){
  futureSubs[,n] = -(GFault %*% futureSlips[,n])
}
futureSub = rowMeans(futureSubs)
g = plotSubsidenceGrid2(futureSub, cbind(lonLatGrid[,1],
                                         lonLatGrid[,2]))
plot(g)
