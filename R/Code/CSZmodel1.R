setwd("~/Uni/NTNU/Masters Project/CSZ/R/Code")

# Building the CSZ model using TMB
library(TMB)

# add the negative log likelihood from the C++ code 
compile("CSZmodel1.cpp")
dyn.load(dynlib("CSZmodel1"))

# which earthquake the model is run for
# Levels: T1 T10 T10R1 T11 T12 T2 T3 T3a T4 T4a T5 T5a T5b T6 T6a T7 T7a T8 T8a T9 T9a
# takes a while to make T1 okada matrix, can use T3 or T4 for quicker run times.
earthquake = "T1"

# load in the subsidence data
setwd("~/Uni/NTNU/Masters Project/CSZ/R/Data/Subsidence")
load("DR3.RData")
setwd("~/Uni/NTNU/Masters Project/CSZ/R/Code")

# load in the fault geometry
# if not saved run the function "getFaultGeometry"
setwd("~/Uni/NTNU/Masters Project/CSZ/R/Data/Fault Geometry")
load("simpleFault.RData") # loads in triangular fault geometry as "simpleFault"
setwd("~/Uni/NTNU/Masters Project/CSZ/R/Code")

# load in the okada matrix
dir = paste0("~/Uni/NTNU/Masters Project/CSZ/R/Data/Okada/OkadaMatrix", earthquake, ".RData")
# should have already made Okada matrix for this earthquake
if (file.exists(dir)){
  setwd("~/Uni/NTNU/Masters Project/CSZ/R/Data/Okada")
  load(paste0("OkadaMatrix", earthquake, ".RData"))
  setwd("~/Uni/NTNU/Masters Project/CSZ/R/Code")
} else {
  # otherwise create the new okada unit vector
  # The summation step is skipped, as this happens when G is matrix multiplied with the slip vector
  G = getOkada(geom = simpleFault,
               lon  = DR3$Lon[DR3$event == earthquake],
               lat  = DR3$Lat[DR3$event == earthquake],
               earthquake = earthquake)
}


# load in the spde mesh and spde structure
# created based on the fault geometry
# if not saved then run "spde_mesh.R" to create the files
setwd("~/Uni/NTNU/Masters Project/CSZ/R/Data/SPDE")
load("spdeMesh.RData") # loads in two spde objects
setwd("~/Uni/NTNU/Masters Project/CSZ/R/Code")


## get the observed subsidence data for the correct earthquake
# Note change sign to minus since Okada model gives negative values for downwards movements.
subsidence = -DR3$subsidence[DR3$event == earthquake]
# I assume that these are the sigma values to be used?
sigma = DR3$Uncertainty[DR3$event == earthquake]


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


## SPDE part: builds 3 components of Q (precision matrix)
# copied from lindgren example
spde = inla_spde$param.inla[c("M0","M1","M2")]


# set up the data list
data = list(depth      = depths,
            subsidence = subsidence,
            sigma      = sigma,
            okada      = as.matrix(G),
            spde_idx   = inla_mesh$idx$loc - 1, # -1 because c++ starts array indexes at 0.
            spde       = spde)


# the parameters to optimise
# I am guessing ar good starting parameters
parameters = list(x          = rep(0, M),
                  log_lambda = -2.1,
                  mu         = 3.7,
                  log_kappa  = -4.6,
                  log_tau    = -1.9)


## create the objective function using TMB
obj = MakeADFun(data,
                parameters,
                random="x",
                DLL="CSZmodel1",
                hessian=TRUE)


## minimize the objective function for the parameters.
#opt = optim(obj$par, obj$fn, obj$gr, method="BFGS")

## it seems like nlminb works well!
opt0 = nlminb(start       =    obj$par,
              objective   =    obj$fn,
              gradient    =    obj$gr)

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
  
r = obj$report(thispar)

effRange = sqrt(8*1)/exp(opt0$par[3])
print(paste("Effective Range: ", round(effRange, 2), "Km"))

# check predicted subsidences
lats  = DR3$Lat[DR3$event == earthquake]
xlim = range(c(r$okadaSubsidence, -DR3$subsidence[DR3$event == earthquake]))
plot(r$okadaSubsidence, lats, xlim=xlim, pch=19, col='blue', cex=.2)
points(-DR3$subsidence[DR3$event == earthquake], lats, pch=19, col='black', cex=.2)

# check underlying spatial field
g = plotFault(simpleFault, z=r$x[inla_mesh$idx$loc], legendTitle="Normalized log slips")
plot(g)

# check average untapered slip
untaperedSlips = r$untaperedSlips
g2 = plotFault(simpleFault, z=r$untaperedSlips, legendTitle="Untapered Slips (m)")
plot(g2)

# plot tapered slips
taperedSlips = r$taperedSlips
g3 = plotFault(simpleFault, z=r$taperedSlips, legendTitle="Tapered Slips (m)")
plot(g3)

# check the mean of a log normal
mu = opt0$par[2]
tau = exp(opt0$par[4])
u = exp(mu + tau/2)
print(is.positive.definite(round(as.matrix(r$Q), 10)))
Q = round(as.matrix(r$Q), 10)
Q1 = solve(Q)
EUntaperedSlip = exp(mu)*exp(diag(Q1) / 2)
print(range(EUntaperedSlip))
print(mean(EUntaperedSlip))

## Get standard errors
SD0 = TMB::sdreport(obj, getJointPrecision=TRUE,
                    bias.correct = TRUE,
                    bias.correct.control = list(sd = TRUE))


# here I should check that I get a PD covariance structure
print(is.positive.definite(as.matrix(SD0$jointPrecision)))
# TRUE!!!

# simulate draws
rmvnorm_prec = function(mu, chol_prec, nSims) {
  z = matrix(rnorm(length(mu) * nSims), ncol=nSims)
  L = chol_prec #Cholesky(prec, super=TRUE)
  z = Matrix::solve(L, z, system = "Lt") ## z = Lt^-1 %*% z
  z = Matrix::solve(L, z, system = "Pt") ## z = Pt    %*% z
  z = as.matrix(z)
  mu + z
}

# take samples from fitted model
mu = c(SD0$par.fixed, SD0$par.random)
nSims = 1000
L = Cholesky(SD0[['jointPrecision']], super = T)
draws = rmvnorm_prec(mu = mu , chol_prec = L, nSims = nSims)

## summarize the draws
parnames = c(names(SD0[['par.fixed']]), names(SD0[['par.random']]))
logLambdaDraws  = draws[parnames == 'log_lambda',]
muDraws = draws[parnames == 'mu',]
logKappaDraws = draws[parnames == 'log_kappa',]
logTauDraws = draws[parnames == 'log_tau',]
xDraws = draws[parnames == 'x',]

# Now I get the distributional results for the parameters
# Histograms - DONE
# Median + SD
# 95% Prediction Interval

draws = data.frame(logLambda = logLambdaDraws,
                   Mu = muDraws,
                   logKappa = logKappaDraws,
                   logTau = logTauDraws)

draws = as.data.frame(pivot_longer(draws, cols = everything(), names_to = "Parameter", values_to = "Value"))

# Create a 2x2 grid of histograms
# cairo_pdf used so that the greek symbols are saved nicely
g = ggplot(draws, aes(x = Value, fill=Parameter)) +
  geom_histogram(bins=30, position = "identity", alpha = 0.8) +
  facet_wrap(~Parameter,
             scales = "free",
             labeller = as_labeller(c('logKappa'='log(\u03ba)',
                                      'logLambda'='log(\u03bb)',
                                      'logTau'='log(\u03c4)',
                                      'Mu'='\u03bc'))) +
  guides(fill = "none") +
  xlab("Parameter Value") +
  ylab("Count")

plot(g)
# With Cairo
ggsave(g, filename = "Parameter Hists.pdf", 
       device = cairo_pdf, width=9, height=9)

draws = data.frame(logLamda = logLambdaDraws,
                   Mu = muDraws,
                   logKappa = logKappaDraws,
                   logTau = logTauDraws)

parSum = cbind(mean=(apply(draws, 2, mean)),
               median = (apply(draws, 2, median)),
               sd     = (apply(draws, 2, sd)),
               lower = (apply(draws, 2, quantile, .05)),
               upper = (apply(draws, 2, quantile, .95)))


# create the slip draws
slipDraws = matrix(data=0, nrow=K, ncol=nSims)

# loop and do all simulations for each subfault
for (i in 1:K){
  taper = exp(-exp(logLambdaDraws) * depths[i]) # singular value
  idx = inla_mesh$idx$loc[i] # also singular value
  slipDraws[i,] = taper * exp(muDraws + xDraws[idx,]) # I need to simulate each mu
}

## find the median and sd across draws, as well as 90% intervals
slipSum = cbind(mean = (apply(slipDraws, 1, mean)),
                median = (apply(slipDraws, 1, median)),
                sd     = (apply(slipDraws, 1, sd)),
                lower = (apply(slipDraws, 1, quantile, .05)),
                upper = (apply(slipDraws, 1, quantile, .95)))

# plot the posterior mean slip
g4 = plotFault(simpleFault, z=slipSum[,1], legendTitle="Mean\nPosterior Slip (m)")
plot(g4)

# plot the posterior median slip
g5 = plotFault(simpleFault, z=slipSum[,2], legendTitle="Median\nPosterior Slip (m)")
plot(g5)

# plot the posterior slip standard deviation
g6 = plotFault(simpleFault, z=slipSum[,3], legendTitle="Stanadard Deviation\nPosterior Slip (m)")
plot(g6)


# now calculate subsidences for each data point
# I use the mean slip across each sub fault
subPred = G %*% slipSum[,1]

absError = abs(subsidence - subPred)
signError = sign(subsidence) == sign(subPred)


plotErrors = function(error, signError, scale=2){
  # gets the base map of the CSZ
  g = plotBase(scale=scale)
  
  Sites = data.frame(Site=DR3$Site[DR3$event == earthquake],
                     Lat=DR3$Lat[DR3$event == earthquake],
                     Lon=DR3$Lon[DR3$event == earthquake])
  Sites$Error = error
  Sites$sign = signError
  
  
  g = g +
    geom_point(data=Sites, aes(x=Lon, y=Lat, colour=Error, shape=sign), size=scale*1.5) +
    scale_colour_gradient(low="red", high="green", name="Absolute Error (m)", trans = "reverse") +
    scale_shape_manual(name = "Same Sign?", values = c(15, 20), labels = c("No", "Yes")) + 
    coord_sf(xlim=-c(128, 122), ylim=c(40, 50))
  
  return(g)
}

g7 = plotErrors(absError, signError)
plot(g7)

print(paste("Count of wrong sign: ", N - sum(signError)))

# Now calculate subsidence across a mesh of the whole fault
# I need to create a new G here

lon = seq(-128, -123, length.out=25)
lat = seq(40, 50, length.out=25)
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
  G2 = getOkada(geom = simpleFault,
                lon  = lonLat$x,
                lat  = lonLat$y,
                earthquake = "T1Whole")
}

subWholeFault = G2 %*% slipSum[,1]

subDF = data.frame(Lon=lonLat$x,
                   Lat=lonLat$y,
                   Sub=subWholeFault)

g8 = plotBase(scale=2, labels=FALSE, countryBoundary=FALSE)
g8 = g8 +
  geom_raster(data=subDF, aes(x=Lon, y=Lat, fill=Sub), alpha=0.75) +
  scale_fill_gradientn(colours=c("red", "white", "blue"), name = "Vertical\nDisplacement (m)",
                       values=scales::rescale(c(min(subDF$Sub), 0, max(subDF$Sub)))) +
  theme(legend.position = "right", legend.key.height = unit(3, 'cm')) +
  coord_sf(xlim=-c(128, 122), ylim=c(40, 50))
plot(g8)
