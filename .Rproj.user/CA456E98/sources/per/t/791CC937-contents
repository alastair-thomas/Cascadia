
# load old rectangular CSZ fault geometry
loadOldCSZgeom = function() {
  # fault geometry data (convert from km to m and to longitude format consistent with other data)
  faultGeom = read.csv("CSZe01.csv")
  faultGeom$longitude = faultGeom$longitude - 360
  kmCols = c(4, 6, 7)
  faultGeom[,kmCols] = faultGeom[,kmCols] * 10^3
  
  csz = divideFault(faultGeom)
  
  csz
}

# loads subsidence data
loadSubsidenceDat = function() {
  
  # subsidence data
  load("data/DR1.RData", envir=globalenv())
  event = dr1$event
  
  # get unique CSZ earthquake event names
  events = as.character(event)
  uniqueEvents = unique(events)
  sortI = rev(c(1, 12, 11, 14, 2, 15, 3, 4, 18, 5, 13, 6, 19, 7, 20, 10, 21, 8, 9, 16, 17)) # T1 is most recent, so reverse (T1 is last)
  uniqueEvents = uniqueEvents[sortI]
  
  # subset uniqueEvents so it only contains major events
  eventSubset = c(1, 2, 4, 6, 8, 10, 12, 15, 17, 19, 20, 21)
  majorEvents = uniqueEvents[eventSubset]
  
  # sort dr1 by event
  dr1$event = factor(dr1$event, levels=uniqueEvents)
  inds = order(dr1$event)
  dr1 = dr1[inds,]
  events = as.character(dr1$event)
  
  # divide uncertainty by 2 because Uncertainty seems to represent 1.96 sigma
  dr1$Uncertainty = dr1$Uncertainty/qnorm(0.975)
  
  allSubDat = dr1
  majorSubDat = allSubDat[as.character(allSubDat$event) %in% majorEvents,]
  
  list(allSubDat=allSubDat, majorSubDat=majorSubDat)
}