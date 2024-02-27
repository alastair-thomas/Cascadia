# load in the fault geometry
# if not saved run the function "getFullFaultGeometry"
setwd("~/Uni/NTNU/Masters Project/CSZ/R/Data/Fault Geometry")
load("simpleFault.RData") # loads in triangular fault geometry as "mediumFault"
setwd("~/Uni/NTNU/Masters Project/CSZ/R/Code")

dips = rep(0, length(simpleFault))
strikes = rep(0, length(simpleFault))
for (i in 1:length(simpleFault)){
  thisDip = simpleFault[[i]]$dip
  thisStrike = simpleFault[[i]]$strike
  
  if (thisStrike > 100){
    thisStrike = thisStrike - 360
  }
  
  dips[i] = thisDip
  strikes[i] = thisStrike
}


g1a = plotFault(simpleFault, z=dips, legendTitle="Dip (°)")
plot(g1a)

g1b = plotGrid(varName="dip", legendTitle="Dip (°)")
plot(g1b)

g2a = plotFault(simpleFault, z=strikes, legendTitle="Strike (°)")
plot(g2a)

g2b = plotGrid(varName="strike", legendTitle="Strike (°)")
plot(g2b)
