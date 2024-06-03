source("~/Uni/NTNU/Masters Project/CSZ/R/Code/Setup.R")
setwd("~/Uni/NTNU/Masters Project/CSZ/R/Code/BayesianCSZ")


# Which earthquake the model is run for
# Levels: T1 T10 T10R1 T11 T12 T2 T3 T3a T4 T4a T5
#         T5a T5b T6 T6a T7 T7a T8 T8a T9 T9a
# T1 == 1700 earthquake
earthquake = "T1"

load("DR4.RData")

# Create the meshes
# First the subfaults
fault    = getFullFaultGeom(max.edge=c(15, 1000), cutoff=5)
S        = length(fault) # number of subfaults
print(S)
plotFault(fault, z=rep(NA, length(fault)))

# Then the SPDE mesh based on the subfaults
spdeStuff = getSPDEMesh(fault, max.edge=c(15, 1000), cutoff=3)
spdeMesh  = spdeStuff$mesh
A         = spdeStuff$A
B         = dim(spdeMesh$loc)[1] # number of basis functions
print(B)
print(spdeMesh)
plotMesh(spdeMesh)

# Visualise both to check they are okay
plotBothMesh(fault, spdeMesh)

# Okada matrix
G = getOkada(fault,
             DR4$Lon[DR4$event == earthquake],
             DR4$Lat[DR4$event == earthquake])

# subsidence data
subsidence = -DR4$subsidence[DR4$event == earthquake]

# uncertainty data
v = DR4$Uncertainty[DR4$event == earthquake]

## extract the depths of the centers of each subfault
depths = rep(0, S) # empty vector for depths
for (i in 1:S){
  # Use positive depths
  # as depth increases, taper goes from one to zero
  depths[i] = -fault[[i]]$depth / (10^3) # change to kms
}

## SPDE part:
# Create the spde object with correct priors
# PC priors on SPDE
# P(range < a) = b
# P(sigma > c) = d
a = 40 # Estimated range of the spatial field
b = 0.05
c = 5 # an ``upper bound'' on the variance
d = 0.05
maternPri = c(a, b, c, d)

# create spde object
spdeINLA = inla.spde2.pcmatern(spdeMesh,
                               alpha=2,
                               prior.range=c(a, b),
                               prior.sigma=c(c, d))
# get the matricies out
spdeMatrix = spdeINLA$param.inla[c("M0","M1","M2")]

# shape found from taperPrior.R
shape = 93.757027137
scale = 0.001316695
taperPri = c(shape, scale)


# set up the data list
data = list(depth      = depths,
            subsidence = subsidence,
            V          = v,
            okada      = as.matrix(G),
            A          = A,
            M0         = spdeMatrix$M0,
            M1         = spdeMatrix$M1,
            M2         = spdeMatrix$M2,
            maternPri  = maternPri,
            taperPri   = taperPri)


# the parameters to optimise
# I am guessing are good starting parameters
parameters = list(x          = rep(0.1, B),
                  logLambda  = -2.1,
                  mu         = 3.6,
                  logKappa   = -4.4,
                  logTau     = 3.7)

save(fault, file="fault.RData")
save(spdeMesh, file="spdeMesh.RData")
save(data, file="BayesianData.RData")
save(parameters, file="BayesianParameters.RData")


load("BayesianData.RData")
load("BayesianParameters.RData")

data = list(depth      = data$depth,
            subsidence = data$subsidence,
            V          = data$V,
            okada      = data$okada,
            A          = data$spde_idx,
            M0         = data$M0,
            M1         = data$M1,
            M2         = data$M2,
            maternPri  = data$maternPri,
            taperPri   = data$taperPri)

# Add the negative log likelihood functions from the TMB code 
TMB::compile("BayesianCSZ.cpp",
             framework="TMBad")
base::dyn.load(dynlib("BayesianCSZ"))

## create the objective function using TMB
TMB::config(tmbad.sparse_hessian_compress = 1)
obj = MakeADFun(data,
                parameters,
                random="x",
                DLL="BayesianCSZ",
                hessian=TRUE)

## Optimize the model parameters
optTime = system.time(opt0 <- nlminb(start     = obj$par,
                                     objective = obj$fn,
                                     gradient  = obj$gr,
                                     control   = list(iter.max  = 100000,
                                                      eval.max  = 100000)))[3]

# Check the state of the optimisation
# Gets important outputs from opt0 and obj.
# Plots subsidences, spatial effects, tapered and untapetred slips
checkOptimisation(obj, opt0, spdeMesh, fault, DR3, earthquake)

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
summariseRandomEffects(xDraws, spdeMesh)

## create the slip draws
slipDraws = matrix(data=0, nrow=S, ncol=nSims)

## loop and do all simulations for each subfault
for (i in 1:S){
  taper = exp(-exp(logLambdaDraws) * depths[i]) # 1000 sims og the kth fault
  idx = spdeMesh$idx$loc[i] # correct x index
  slipDraws[i,] = taper * exp(muDraws + xDraws[idx,]) # 1000 slips
}

## find the median and sd across draws, as well as 90% intervals
slipSum = cbind(mean = (apply(slipDraws, 1, mean)),
                median = (apply(slipDraws, 1, median)),
                sd     = (apply(slipDraws, 1, sd)),
                lower = (apply(slipDraws, 1, quantile, .05)),
                upper = (apply(slipDraws, 1, quantile, .95)))
slipSum = cbind(slipSum,
                width = slipSum[,5] - slipSum[,4])

# plot the posterior mean slip
g7 = plotFault(fault, z=slipSum[,1], legendTitle="Posterior Slip\nMean (m)")
plot(g7)

# plot the posterior median slip
g8 = plotFault(fault, z=slipSum[,2], legendTitle="Posterior Slip\nMedian (m)")
plot(g8)

# plot the posterior slip standard deviation
g9 = plotFault(fault,
               z=slipSum[,3],
               legendTitle="Posterior Slip\nStandard Deviation (m)",
               colourScale="plasma")
plot(g9)

# plot the posterior slip width
g9b = plotFault(fault,
               z=slipSum[,6],
               legendTitle="Posterior Slip\n95% Width (m)",
               colourScale="mako")
plot(g9b)

# Now I find the magnitude of all the simulated slips
magnitudes = rep(0, nSims)
for (i in 1:nSims){
  magnitudes[i] = getMomentFromSlip(slips = slipDraws[,i],
                                    fault = fault)
}

# Plot the resulting magnitudes as a histogram
# Add line for mean
g10 = summariseMagnitudes(magnitudes)
plot(g10)




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

