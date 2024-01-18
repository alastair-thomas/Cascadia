# Script to load in new subsidence data from Humboldt Bay
# Data from Padgett2021 "Timing and amount of southern Cascadia earthquake subsidence over
# the past 1700 years at northern Humboldt Bay, California, USA"

# A=T1
# C=T3
# D=T5
#
# MD.3:  (-124.1114, 40.8578)
# MD.6:  (-124.1041, 40.8601)
# MD.13: (-124.1025, 40.8604)

Region = rep("California", 3)
Site = rep("Humboldt Bay", 3)
Lat = c(40.858, 40.860, 40.860)
Lon = c(-124.111, -124.104, -124.103)
event = c("T1", "T3", "T5")
subsidence = c(0.85, 0.42, 0.79)
Uncertainty = c(0.46, 0.37, 0.47)
source = rep(49, 3)
quality = rep(1, 3) # is it quality 1? It uses the bayesian transfer function?

dr2 = data.frame(Region=Region,
                 Site=Site,
                 Lat=Lat,
                 Lon=Lon,
                 event=event,
                 subsidence=subsidence,
                 Uncertainty=Uncertainty,
                 source=as.factor(source),
                 quality=as.factor(quality))

DR3 = bind_rows(dr1, dr2)
save(DR3, file="DR3.RData")
