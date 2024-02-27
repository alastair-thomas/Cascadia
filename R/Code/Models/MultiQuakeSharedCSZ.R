setwd("~/Uni/NTNU/Masters Project/CSZ/R/Code")

# Building the CSZ model using TMB
library(TMB)

# Add the negative log likelihood functions from the TMB code 
compile("MultiQuakeSharedCSZ.cpp")
dyn.load(dynlib("MultiQuakeSharedCSZ"))


# load in the subsidence data
setwd("~/Uni/NTNU/Masters Project/CSZ/R/Data/Subsidence")
load("DR3.RData")
setwd("~/Uni/NTNU/Masters Project/CSZ/R/Code")


# load in the fault geometry
# if not saved run the function "getFullFaultGeometry"
setwd("~/Uni/NTNU/Masters Project/CSZ/R/Data/Fault Geometry")
load("simpleFault.RData") # loads in triangular fault geometry as "mediumFault"
setwd("~/Uni/NTNU/Masters Project/CSZ/R/Code")
K = length(simpleFault) # number of subfaults

# Creat the correct data for each earthquake

G = list() # the list of Okada matrices
subsidences = list() # list of subsidence vectors
v = list() # list of uncertainties

earthquakes = c("T1","T2","T3","T4","T5","T6","T7","T8","T9","T10","T11","T12")
M = length(earthquakes) # number of earthquakes

for (i in 1:M){
  
  thisLon = DR3$Lon[DR3$event == earthquakes[i]]
  thisLat = DR3$Lat[DR3$event == earthquakes[i]]
  thisSub = DR3$subsidence[DR3$event == earthquakes[i]]
  thisUnc = DR3$Uncertainty[DR3$event == earthquakes[i]]
  
  dir = paste0("~/Uni/NTNU/Masters Project/CSZ/R/Data/Okada/Simple MultiQuake/G", i, ".RData")
  if (file.exists(dir)){
    load(dir)
  } else{
    thisG = getOkada(simpleFault, thisLon, thisLat)
    save(thisG, file=dir)
  }
  
  G[[i]] = as.matrix(thisG)
  subsidences[[i]] = -thisSub # remember the change of sign
  v[[i]] = thisUnc
}


# load in the spde mesh and spde structure
# created based on the fault geometry
# if not saved then run "spde_mesh.R" to create the files
setwd("~/Uni/NTNU/Masters Project/CSZ/R/Data/SPDE")
load("spdeMesh.RData") # loads in two spde objects
setwd("~/Uni/NTNU/Masters Project/CSZ/R/Code")
B = dim(inla_mesh$loc)[1] # number of basis functions

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
a = 200 # Estimated range of the spatial field
b = 0.5
c = 1 # an ``upper bound'' on the variance
d = 0.05
maternPri = c(a, b, c, d)

# create spde object
inla_spde = inla.spde2.pcmatern(inla_mesh,
                                alpha=2,
                                prior.range=c(a, b),
                                prior.sigma=c(c, d))
# get the matricies out
spdeMatrix = inla_spde$param.inla[c("M0","M1","M2")]

# shape found from taperPrior.R
# means that 95% of the time the taper drops to 0.05 by 20km.
shape = 93.757027137
scale = 0.001316695
taperPri = c(shape, scale)


# set up the data list
data = list(depth      = depths,
            subsidence = subsidences,
            v          = v,
            okada      = G,
            spdeIDX    = inla_mesh$idx$loc - 1, # -1 because c++ starts array indexes at 0.
            M0         = spdeMatrix$M0,
            M1         = spdeMatrix$M1,
            M2         = spdeMatrix$M2,
            maternPri  = maternPri,
            taperPri   = taperPri)


# the parameters to optimise
# I am guessing are good starting parameters
parameters = list(X          = matrix(0.5, nrow=B, ncol=M),
                  w          = rep(0, B),
                  logLambda  = -2.1,
                  mu         = 3.6,
                  logKappaX  = -4.4,
                  logTauX    = 3.7,
                  logKappaW  = -4.4,
                  logTauW    = 3.7)


## create the objective function using TMB
obj = MakeADFun(data,
                parameters,
                random="X",
                DLL="MultiQuakeSharedCSZ",
                hessian=TRUE)

## Optimize the model parameters
optTime = system.time(opt0 <- nlminb(start     = obj$par,
                                     objective = obj$fn,
                                     gradient  = obj$gr,
                                     control   = list(iter.max  = 100000,
                                                      eval.max  = 100000)))[3]

# maybe need some wrappers for inner and outer.
# I should add reports to compare values, make sure it's all correct

print(optTime)

# Check the state of the optimisation
# Gets important outputs from opt0 and obj.
# Plots subsidences, spatial effects, tapered and untapetred slips
checkOptimisation(obj, opt0, inla_mesh, simpleFault, DR3, earthquake)

# If everything looks good continue and get results from the model.

## Get standard errors via SD report
sdTime = system.time(SD0 <- TMB::sdreport(obj, getJointPrecision=TRUE,
                                          bias.correct = TRUE,
                                          bias.correct.control = list(sd = TRUE)))[3]

# summarise the SD report
summary(SD0, "report")

# Check the joint precision is positive definite
print(is.positive.definite(as.matrix(SD0$jointPrecision)))
# TRUE

# Get draws from all the model parameters
draws = simulateCSZ(SD0, nSims=1000)

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
slipDraws = matrix(data=0, nrow=K, ncol=1000)

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


# now calculate subsidences for each data point
# I use the mean slip across each sub fault
subPred = G %*% slipSum[,1]

absError = abs(subsidence - subPred)
signError = sign(subsidence) == sign(subPred)

g10 = plotErrors(absError, signError)
plot(g10)

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
g11 = plotBase(labels=FALSE, countryBoundary=FALSE)
g11 = g11 +
  geom_tile(data=subDF, aes(x=Lon, y=Lat, fill=Sub)) +
  scale_fill_scico(alpha=0.75, palette = 'vik', midpoint=0, name="Subsidence (m)") +
  #scale_fill_viridis_c(alpha=0.75, name="Subsidence (m)", option="viridis") +
  theme(legend.position = "right", legend.key.height = unit(3, 'cm')) +
  coord_sf(xlim=-c(128, 122), ylim=c(40, 50))

plot(g11)

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

g12 = plotBase(labels=FALSE, countryBoundary=FALSE)
g12 = g12 +
  geom_tile(data=subDF, aes(x=Lon, y=Lat, fill=sd, alpha=sd)) +
  scale_fill_viridis_c(name="Subsidence\nStandard Deviation (m)",
                       option="plasma", direction=-1) +
  theme(legend.position = "right", legend.key.height = unit(1, 'cm')) +
  coord_sf(xlim=-c(128, 123), ylim=c(40, 50)) +
  guides(alpha="none")

plot(g12)

