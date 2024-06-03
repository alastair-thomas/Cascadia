# set the correct working directory
setwd("~/Uni/NTNU/Masters Project/CSZ/R/Code/MultiQuakeCSZ")
# load in the subsidence data
load("DR4.RData")

# Create the meshes
# First the subfaults
fault    = getFullFaultGeom(max.edge=c(75, 1000), cutoff=25)
S        = length(fault) # number of subfaults
print(S)

depths = rep(0, S)
for (i in 1:S){
  depths[i] = -fault[[i]]$depth / 1000
}
plotFault(fault, z=depths, legendTitle="Depth (km)", colourScale="mako")
plotFault(fault, z=rep(NA, length(fault)))

# Then the SPDE mesh based on the subfaults
spdeStuff = getSPDEMesh(fault, max.edge=c(75, 1000), cutoff=15)
spdeMesh = spdeStuff$mesh
A        = spdeStuff$A # projection matrix
B        = dim(spdeMesh$loc)[1] # number of basis functions

# show results of SPDE mesh
print(spdeMesh)

#plot SPDE Mesh
plotMesh(spdeMesh)
plotBothMesh(fault, spdeMesh)

# Create the correct data for each earthquake
G = list() # the list of Okada matrices
subsidences = list() # list of subsidence vectors
v = list() # list of uncertainties

earthquakes = c("T1","T2","T3","T4","T5","T6","T7","T8","T9","T10","T11","T12")
E = length(earthquakes) # number of earthquakes

for (e in 1:E){
  
  print(e)
  
  thisLon = DR4$Lon[DR4$event == earthquakes[e]]
  thisLat = DR4$Lat[DR4$event == earthquakes[e]]
  thisSub = DR4$subsidence[DR4$event == earthquakes[e]]
  thisUnc = DR4$Uncertainty[DR4$event == earthquakes[e]]
  thisG = getOkada(fault, thisLon, thisLat)
  
  G[[e]] = as.matrix(thisG)
  subsidences[[e]] = thisSub # I change the sign of okada model
  v[[e]] = thisUnc
}


## SPDE part:
# Create the spde object with correct priors
# PC priors on SPDE
# P(range < a) = b
# P(sigma > c) = d
a = 100 # Estimated range of the spatial field
b = 0.5
c = 1 # an ``upper bound'' on the variance
d = 0.01
maternPri = c(a, b, c, d)

# create spde object
spdeINLA = inla.spde2.pcmatern(spdeMesh,
                               alpha=2,
                               prior.range=c(a, b),
                               prior.sigma=c(c, d))
# get the matrices out
spdeMatrix = spdeINLA$param.inla[c("M0","M1","M2")]

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
            A          = A,
            M0         = spdeMatrix$M0,
            M1         = spdeMatrix$M1,
            M2         = spdeMatrix$M2,
            maternPri  = maternPri,
            taperPri   = taperPri)

setwd("~/Uni/NTNU/Masters Project/CSZ/R/Code/MultiQuakeCSZ")
save(data, file="Data.RData")
save(fault, file="Fault.RData")
save(spdeMesh, file="spdeMesh.RData")
