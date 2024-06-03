setwd("~/Uni/NTNU/Masters Project/CSZ/R/Code/MultiQuakeAnisoCSZ")

# load in the data
load("DR4.RData")

# Create the meshes
# First the subfaults
fault    = getFullFaultGeom(max.edge=c(75, 1000), cutoff=25)
S        = length(fault) # number of subfaults
print(S)
plotFault(fault, z=rep(NA, length(fault)))

## extract the depths of the centers of each subfault
depths = rep(0, S) # empty vector for depths
for (i in 1:S){
  # Use positive depths
  # as depth increases, taper goes from one to zero
  depths[i] = -fault[[i]]$depth / (10^3) # change to kms
}
plotFault(fault, z=depths, legendTitle="Depth (km)", colourScale="magma")

# Then the SPDE mesh based on the subfaults
spdeStuff = getSPDEMesh(fault, max.edge=c(75, 1000), cutoff=15)
spdeMesh  = spdeStuff$mesh
A         = spdeStuff$A # projection matrix
B         = dim(spdeMesh$loc)[1] # number of basis functions
print(B)
print(spdeMesh)
plotMesh(spdeMesh)

# Visualise both to check they are okay
plotBothMesh(fault, spdeMesh)

# Create the correct data for each earthquake
G = list() # the list of Okada matrices
subsidences = list() # list of subsidence vectors
v = list() # list of uncertainties

earthquakes = c("T1","T2","T3","T4","T5","T6","T7","T8","T9","T10","T11","T12")
E = length(earthquakes) # number of earthquakes
for (e in 1:E){
  
  thisLon = DR4$Lon[DR4$event == earthquakes[e]]
  thisLat = DR4$Lat[DR4$event == earthquakes[e]]
  thisSub = DR4$subsidence[DR4$event == earthquakes[e]]
  thisUnc = DR4$Uncertainty[DR4$event == earthquakes[e]]
  thisG = getOkada(fault, thisLon, thisLat)
  
  G[[e]] = as.matrix(thisG)
  subsidences[[e]] = thisSub
  v[[e]] = thisUnc
}


## SPDE part:
# Create the spde object with correct priors
# PC priors on SPDE
# P(range < a) = b
# P(sigma > c) = d
a = 100 # Estimated minimum range of the spatial field
b = 0.50
c = 1 # an ``upper bound'' on the variance
d = 0.01
maternPri = c(a, b, c, d)

# create spde object
spdeINLA = inla.spde2.pcmatern(spdeMesh,
                               alpha=2,
                               prior.range=c(a, b),
                               prior.sigma=c(c, d))

# shape found from taperPrior.R
# means that 95% of the time the taper drops to 0.05 by 20km.
shape = 93.757027137
scale = 0.001316695
taperPri = c(shape, scale)

# Anisotropic prior found from anisoPrior.R
anisoPri = 0.38


# set up the data list
data = list(depth       = depths,
            subsidence  = subsidences,
            v           = v,
            okada       = G,
            A           = A,
            maternPri   = maternPri,
            taperPri    = taperPri,
            anisoPri    = anisoPri)

# spdeMesh$idx$loc - 1

# ---------- Start code that prepare objects for anisotropy. Should not be modified when you make your own example
Dset = 1:2
# Triangle info
TV = spdeMesh$graph$tv           # Triangle to vertex indexing
V0 = spdeMesh$loc[TV[,1],Dset]   # V = vertices for each triangle
V1 = spdeMesh$loc[TV[,2],Dset]
V2 = spdeMesh$loc[TV[,3],Dset]
E0 = V2 - V1                      # E = edge for each triangle
E1 = V0 - V2
E2 = V1 - V0  
# Calculate Areas
TmpFn = function(Vec1, Vec2) abs(det( rbind(Vec1, Vec2) ))
Tri_Area = rep(NA, nrow(E0))
for(i in 1:length(Tri_Area)) Tri_Area[i] = TmpFn( E0[i,],E1[i,] )/2   # T = area of each triangle

data$spde <- list(
  "n_s"      = spdeINLA$n.spde,
  "n_tri"    = nrow(TV),
  "Tri_Area" = Tri_Area,
  "E0"       = E0,
  "E1"       = E1,
  "E2"       = E2,
  "TV"       = TV - 1,
  "G0"       = spdeINLA$param.inla$M0,
  "G0_inv"   = as(diag(1/diag(spdeINLA$param.inla$M0)), "TsparseMatrix"))
# ---------- End code that prepare objects for anisotropy.


# change to what you want to save as
setwd("~/Uni/NTNU/Masters Project/CSZ/R/Code/MultiQuakeAnisoCSZ")
save(fault, file="Fault.RData")
save(spdeMesh, file="spdeMesh.RData")
save(data, file="Data.RData")
