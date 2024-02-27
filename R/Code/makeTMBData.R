# This script creates the data object to be passed to TMB
# It works for the BayesianCSZ Model
# Currently earthquake T1
# Currently on a very detailed fault geometry

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
load("fullFault.RData") # loads in triangular fault geometry as "mediumFault"
setwd("~/Uni/NTNU/Masters Project/CSZ/R/Code")


# load in the okada matrix
# If not created use the "getokada" function with the current fault data
setwd("~/Uni/NTNU/Masters Project/CSZ/R/Data/Okada")
load("G_fullFault_T1.RData")
setwd("~/Uni/NTNU/Masters Project/CSZ/R/Code")


# load in the spde mesh and spde structure
# created based on the fault geometry
# if not saved then run "spde_mesh.R" to create the files
setwd("~/Uni/NTNU/Masters Project/CSZ/R/Data/SPDE")
load("spdeMesh_fullFault.RData") # loads in two spde objects
setwd("~/Uni/NTNU/Masters Project/CSZ/R/Code")


## get the observed subsidence data for the correct earthquake
# Note change sign to minus since Okada model gives negative values for downwards movements.
subsidence = -DR3$subsidence[DR3$event == earthquake]
# get the variance in data
V = DR3$Uncertainty[DR3$event == earthquake]


N = length(subsidence) # number of data points
K = length(fullFault) # number of subfaults
M = dim(spdeMesh$loc)[1] # number of points in spde mesh


## extract the depths of the centers of each subfault
depths = rep(0, K) # empty vector for depths
for (i in 1:K){
  # Use positive depths
  # as depth increases, taper goes from one to zero
  depths[i] = -fullFault[[i]]$depth / (10^3) # change to kms
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
inla_spde = inla.spde2.pcmatern(spdeMesh,
                                alpha=2,
                                prior.range=c(a, b),
                                prior.sigma=c(c, d))
# get the matricies out
spdeMatrix = inla_spde$param.inla[c("M0","M1","M2")]

shape = 1.413327
scale = 1
taperPri = c(shape, scale)


# set up the data list
data = list(depth      = depths,
            subsidence = subsidence,
            V          = V,
            okada      = as.matrix(G),
            spde_idx   = spdeMesh$idx$loc - 1, # -1 because c++ starts array indexes at 0.
            M0         = spdeMatrix$M0,
            M1         = spdeMatrix$M1,
            M2         = spdeMatrix$M2,
            maternPri  = maternPri,
            taperPri   = taperPri)

