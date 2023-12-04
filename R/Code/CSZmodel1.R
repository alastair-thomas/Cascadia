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
load("DR1.RData")
setwd("~/Uni/NTNU/Masters Project/CSZ/R/Code")


# load in the fault geometry
# if not saved run the function "getFaultGeometry"
setwd("~/Uni/NTNU/Masters Project/CSZ/R/Data/Fault Geometry")
load("fault.RData") # loads in triangular fault geometry as triGeomFault
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
  G = getOkada(geom = triGeomFull,
               lon  = dr1$Lon[dr1$event == earthquake],
               lat  = dr1$Lat[dr1$event == earthquake],
               earthquake = earthquake)
}


# load in the spde mesh and spde structure
# created based on the fault geometry
# if not saved then run "spde_mesh.R" to create the files
setwd("~/Uni/NTNU/Masters Project/CSZ/R/Data/SPDE")
load("spde_mesh.RData") # loads in triangular fault geometry as triGeomFault
setwd("~/Uni/NTNU/Masters Project/CSZ/R/Code")


## get the observed subsidence data for the correct earthquake
subsidence = dr1$subsidence[dr1$event == earthquake]
# I assume that these are the sigma values to be used?
sigma = dr1$Uncertainty[dr1$event == earthquake]


N = length(subsidence) # number of data points
K = length(triGeomFull) # number of subfaults
M = dim(inla_mesh$loc)[1] # number of points in spde mesh


## extract the depths of the centers of each subfault
depths = rep(0, K) # empty vector for depths
for (i in 1:K){
  # Use positive depths
  # as depth increases, taper goes from one to zero
  depths[i] = -triGeomFull[[i]]$depth
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
parameters = list(x          = rep(-1, M),
                  log_lambda = log(1/10000),
                  mu         = 3,
                  log_kappa  = 2.5,
                  log_tau    = -2.0)


## create the objective function using TMB
# is random="x"?
# it was in the lindgren paper
# it doesn't work with random="x"
obj = MakeADFun(data, parameters, random="x", DLL="CSZmodel1")


## minimize the objective function for the parameters.
# could use other minimisation function here
# maybe add in lower and upper bounds - need to know what they are...!

L = c(rep(-Inf, M), -Inf, -Inf,   2.0, -3.0)
U = c(rep(Inf, M),  Inf, Inf, 3.0, -1.0)

opt = optim(obj$par, obj$fn, obj$gr, method="L-BFGS-B", lower=L, upper=U)

# Calculate standard deviations, and extract rho
#Rep = sdreport(obj, getJointPrecision=TRUE, bias.correct=TRUE)
#rho_est = summary(Rep,"report")


## Now extract and plot slips

n = length(opt$par)
# need taper function
taper = function(depth, lambda=exp(opt$par[[n-3]])){
  return(exp(-lambda*depth))
}

s = rep(0, K)

t = taper(depths)

mu = opt$par[[n-2]]

for (i in 1:K){
  s[i] = t[i] * exp(mu + opt$par[[inla_mesh$idx$loc[i]]])
}

lon = rep(0, K)
lat = rep(0, K)

for (i in 1:K){
  lon[i]=triGeomFull[[i]]$lon
  lat[i]=triGeomFull[[i]]$lat
}

data = data.frame(x=lon, y=lat, Slip=s)

ggplot(data)+
  geom_point(aes(x=x, y=y, color=Slip))+
  scale_color_viridis_c()

# plotting subsidence data
data2 = data.frame(x=dr1$Lon[dr1$event == earthquake],
                   y=dr1$Lat[dr1$event == earthquake],
                   z=dr1$subsidence[dr1$event==earthquake])
ggplot(data2)+
  geom_point(aes(x=x, y=y, color=z))+
  scale_color_viridis_c()

# calcualte subsidence estimates
# subsidence is a normal variable
# Can calculate the Gs subsidence vector?

y = G %*% s

res = y - dr1$subsidence[dr1$event==earthquake]

plot(res)
plot(dr1$subsidence[dr1$event==earthquake], y)

# we get negative subsidence
