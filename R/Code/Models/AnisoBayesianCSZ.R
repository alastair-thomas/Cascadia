setwd("~/Uni/NTNU/Masters Project/CSZ/R/Code")

# Building the CSZ model using TMB
library(TMB)

# Add the negative log likelihood functions from the TMB code 
compile("BayesianCSZ.cpp")
dyn.load(dynlib("BayesianCSZ"))


# Which earthquake the model is run for
# Levels: T1 T10 T10R1 T11 T12 T2 T3 T3a T4 T4a T5
#         T5a T5b T6 T6a T7 T7a T8 T8a T9 T9a
# T1 == 1700 earthquake
earthquake = "T1"


# load in the subsidence data
setwd("~/Uni/NTNU/Masters Project/CSZ/R/Data/Subsidence")
load("DR3.RData")
setwd("~/Uni/NTNU/Masters Project/CSZ/R/Code")


# load in the fault geometry
# if not saved run the function "getFullFaultGeometry"
setwd("~/Uni/NTNU/Masters Project/CSZ/R/Data/Fault Geometry")
load("simpleFault.RData") # loads in triangular fault geometry as "mediumFault"
setwd("~/Uni/NTNU/Masters Project/CSZ/R/Code")


# load in the okada matrix
# If not created use the "getokada" function with the current fault data
setwd("~/Uni/NTNU/Masters Project/CSZ/R/Data/Okada")
load(paste0("OkadaMatrix", earthquake, ".RData"))
setwd("~/Uni/NTNU/Masters Project/CSZ/R/Code")


# load in the spde mesh and spde structure
# created based on the fault geometry
# if not saved then run "spde_mesh.R" to create the files
setwd("~/Uni/NTNU/Masters Project/CSZ/R/Data/SPDE")
load("spdeMesh.RData") # loads in two spde objects
setwd("~/Uni/NTNU/Masters Project/CSZ/R/Code")


## get the observed subsidence data for the correct earthquake
# Note change sign to minus since Okada model gives negative values for downwards movements.
subsidence = -DR3$subsidence[DR3$event == earthquake]
# get the variance in data
V = DR3$Uncertainty[DR3$event == earthquake]


N = length(subsidence) # number of data points
K = length(simpleFault) # number of subfaults
M = dim(inla_mesh$loc)[1] # number of points in spde mesh


## extract the depths of the centers of each subfault
depths = rep(0, K) # empty vector for depths
for (i in 1:K){
  # Use positive depths
  # as depth increases, taper goes from one to zero
  depths[i] = -simpleFault[[i]]$depth / (10^3) # change to kms
}


## SPDE part:
# Create the spde object with correct priors
# PC priors on SPDE
# P(range < a) = b
# P(sigma > c) = d
a = 150 # Estimated range of the spatial field
b = 0.5
c = 1 # an ``upper bound'' on the variance
d = 0.10
maternPri = c(a, b, c, d)

# create spde object
inla_spde = inla.spde2.pcmatern(inla_mesh,
                                alpha=2,
                                prior.range=c(a, b),
                                prior.sigma=c(c, d))
# get the matricies out
spdeMatrix = inla_spde$param.inla[c("M0","M1","M2")]

# shape found from taperPrior.R
shape = 93.757027137
scale = 0.001316695
taperPri = c(shape, scale)


# set up the data list
data = list(depth      = depths,
            subsidence = subsidence,
            V          = V,
            okada      = as.matrix(G),
            spde_idx   = inla_mesh$idx$loc - 1, # -1 because c++ starts array indexes at 0.
            M0         = spdeMatrix$M0,
            M1         = spdeMatrix$M1,
            M2         = spdeMatrix$M2,
            maternPri  = maternPri,
            taperPri   = taperPri)


# the parameters to optimise
# I am guessing are good starting parameters
parameters = list(x          = rep(0, M),
                  logLambda  = -2.1,
                  mu         = 3.6,
                  logKappa   = -4.4,
                  logTau     = 3.7)


## create the objective function using TMB
obj = MakeADFun(data,
                parameters,
                random="x",
                DLL="BayesianCSZ",
                hessian=TRUE)

## Optimize the model parameters
optTime = system.time(opt0 <- nlminb(start = obj$par,
                                      objective = obj$fn,
                                      gradient  = obj$gr))[3]

# Check the state of the optimisation
# Gets important outputs from opt0 and obj.
# Plots subsidences, spatial effects, tapered and untapetred slips
checkOptimisation(obj, opt0, inla_mesh, simpleFault, DR3, earthquake)

# If everything looks good continue and get results from the model.

## Get standard errors via SD report
SD0 = TMB::sdreport(obj, getJointPrecision=TRUE,
                    bias.correct = TRUE,
                    bias.correct.control = list(sd = TRUE))

# summarise the SD report
summary(SD0, "report")

# Check the joint precision is positive definite
print(is.positive.definite(as.matrix(SD0$jointPrecision)))
# TRUE

# Get draws from all the model parameters
nSims = 1000
draws = simulateCSZ(SD0, nSims=nSims)

## extract the draws for each parameter the draws
parnames = c(names(SD0[['par.fixed']]), names(SD0[['par.random']]))
logLambdaDraws  = draws[parnames == 'logLambda',]
muDraws = draws[parnames == 'mu',]
logKappaDraws = draws[parnames == 'logKappa',]
logTauDraws = draws[parnames == 'logTau',]
xDraws = draws[parnames == 'x',]


# A dataframe for all the fixed parameters
fixedDraws = data.frame(logLambda = logLambdaDraws,
                        Mu = muDraws,
                        logKappa = logKappaDraws,
                        logTau = logTauDraws)

# call a function that gives histograms and basic statistics of fixed pars
summariseFixedParameters(fixedDraws, histLabel=c('logKappa'='log(\u03ba)',
                                                 'logLambda'='log(\u03bb)',
                                                 'logTau'='log(\u03c4)',
                                                 'Mu'='\u03bc'))

# draws plots of the spatial effect and standard deviation
summariseRandomEffects(xDraws, inla_mesh)

## create the slip draws
slipDraws = matrix(data=0, nrow=K, ncol=nSims)

## loop and do all simulations for each subfault
for (i in 1:K){
  taper = exp(-exp(logLambdaDraws) * depths[i]) # 1000 sims og the kth fault
  idx = inla_mesh$idx$loc[i] # correct x index
  slipDraws[i,] = taper * exp(muDraws + xDraws[idx,]) # 1000 slips
}

## find the median and sd across draws, as well as 90% intervals
slipSum = cbind(mean = (apply(slipDraws, 1, mean)),
                median = (apply(slipDraws, 1, median)),
                sd     = (apply(slipDraws, 1, sd)),
                lower = (apply(slipDraws, 1, quantile, .05)),
                upper = (apply(slipDraws, 1, quantile, .95)))

# plot the posterior mean slip
g7 = plotFault(simpleFault, z=slipSum[,1], legendTitle="Posterior Slip\nMean (m)")
plot(g7)

# plot the posterior median slip
g8 = plotFault(simpleFault, z=slipSum[,2], legendTitle="Posterior Slip\nMedian (m)")
plot(g8)

# plot the posterior slip standard deviation
g9 = plotFault(simpleFault,
               z=slipSum[,3],
               legendTitle="Posterior Slip\nStandard Deviation (m)",
               colourScale="plasma")
plot(g9)

# Now I find the magnitude of all the simulated slips
magnitudes = rep(0, nSims)
for (i in 1:nSims){
  magnitudes[i] = getMomentFromSlip(slips = slipDraws[,i],
                                    fault = simpleFault)
}

# Plot the resulting magnitudes as a histogram
# Add line for mean
g10 = summariseMagnitudes(magnitudes)
plot(g10)

#badSlip = slipDraws[,234]

#testg = plotFault(simpleFault, z=badSlip, legendTitle="Slip (m)")




# now calculate subsidences for each data point
# I use the mean slip across each sub fault
subPred = G %*% slipSum[,1]

absError = abs(subsidence - subPred)
signError = sign(subsidence) == sign(subPred)

g11 = plotErrors(absError, signError)
plot(g11)

print(paste("Count of wrong sign: ", N - sum(signError)))

# Now calculate subsidence across a mesh of the whole fault

gridN = 25
lon = seq(-128, -123, length.out=gridN)
lat = seq(40, 50, length.out=gridN)
lonLat = expand.grid(x=lon, y=lat)

# load in the okada matrix for the whole fault
dir = paste0("~/Uni/NTNU/Masters Project/CSZ/R/Data/Okada/OkadaMatrixT1Whole.RData")
# should have already made Okada matrix for this earthquake
if (file.exists(dir)){
  setwd("~/Uni/NTNU/Masters Project/CSZ/R/Data/Okada")
  load("OkadaMatrixT1Whole.RData")
  G2 = G
  load("OkadaMatrixT1.RData")
  setwd("~/Uni/NTNU/Masters Project/CSZ/R/Code")
  
  
} else {
  
  # otherwise create the new okada unit vector
  # The summation step is skipped, as this happens when G is matrix multiplied with the slip vector
  G4 = getOkada(geom = mediumFault,
                lon  = lonLat$x,
                lat  = lonLat$y,
                earthquake = "T1Whole")
  
  setwd("~/Uni/NTNU/Masters Project/CSZ/R/Data/Okada")
  save(G4, file="OkadaMatrixMediumT1Whole.RData")
  setwd("~/Uni/NTNU/Masters Project/CSZ/R/Code")
  
}

# predict subsidence
subWholeFault = G2 %*% slipSum[,1]

# set up for plotting
subDF = data.frame(Lon=lonLat$x,
                   Lat=lonLat$y,
                   Sub=subWholeFault)

# here I use the Scico package which has diverging colour scales
# that are also safe to use for blind people.
g12 = plotBase(labels=FALSE, countryBoundary=FALSE)
g12 = g12 +
  geom_tile(data=subDF, aes(x=Lon, y=Lat, fill=Sub)) +
  scale_fill_scico(alpha=0.75, palette = 'vik', midpoint=0, name="Subsidence (m)") +
  #scale_fill_viridis_c(alpha=0.75, name="Subsidence (m)", option="viridis") +
  theme(legend.position = "right", legend.key.height = unit(3, 'cm')) +
  coord_sf(xlim=-c(128, 122), ylim=c(40, 50))

plot(g12)

# Now I calculate subsidence draws across the whole grid so I can get
# Uncertainties and prediction intervals

subsidenceDraws = matrix(0, nrow=gridN^2, ncol=nSims)
for (i in 1:nSims){
  subsidenceDraws[,i] = G2 %*% slipDraws[,i]
}

subsidenceSum = cbind(mean = (apply(subsidenceDraws, 1, mean)),
                      median = (apply(subsidenceDraws, 1, median)),
                      sd     = (apply(subsidenceDraws, 1, sd)),
                      lower = (apply(subsidenceDraws, 1, quantile, .05)),
                      upper = (apply(subsidenceDraws, 1, quantile, .95)))

subDF = cbind(subDF, subsidenceSum)

g13 = plotBase(labels=FALSE, countryBoundary=FALSE)
g13 = g13 +
  geom_tile(data=subDF, aes(x=Lon, y=Lat, fill=sd, alpha=sd)) +
  scale_fill_viridis_c(name="Subsidence\nStandard Deviation (m)",
                       option="plasma", direction=-1) +
  theme(legend.position = "right", legend.key.height = unit(1, 'cm')) +
  coord_sf(xlim=-c(128, 123), ylim=c(40, 50)) +
  guides(alpha="none")

plot(g13)

