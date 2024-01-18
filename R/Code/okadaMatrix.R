getUnitMatrix = function(fault, x, y, 
                         rake=rep(90, length(fault))) {
  
  slip=rep(1, length(fault))
  
  if(length(rake) == 1) {
    rake = rep(rake, length(fault))
  }
  
  # add slip, rake to subfaults
  for(i in 1:length(fault)) {
    fault[[i]]$rake = rake[i]
    fault[[i]]$slip = slip[i]
  }
  
  # get topographic deformations for all subfaults
  allTopos = lapply(fault, okadaSubfaultTri, x=x, y=y)
  
  return(allTopos)
}

# returns the Okada matrix for the specified data locations and subfault geometry
getOkada = function(geom, lon, lat, earthquake){
  allTopos = getUnitMatrix(geom, x=lon, y=lat)
  
  # allTopos[[i]]$dZ in an NxN matrix which gives the effect of the ith subfault
  # on each of the locations in a mesh of lon, lat
  # Extracting diagonal gives for the pairs of (lon, lat) gives.
  # could sum it up, then have just an okada vector
  
  K = length(allTopos) # number of subfaults
  N = length(lon) # number of data points
  
  # the actual okada matrix
  G = matrix(data=0,
             nrow=N,
             ncol=K)
  
  for (i in 1:K){
    # extract the diagonals of the matrix as these represent the data points
    G[,i] = diag(allTopos[[i]]$dZ)
  }
  
  # save matrix for next time
  setwd("~/Uni/NTNU/Masters Project/CSZ/R/Data/Okada")
  save(G, file=paste0("OkadaMatrix", earthquake, ".RData"))
  setwd("~/Uni/NTNU/Masters Project/CSZ/R/Code")
  
  return(G)
}

testOkadaSign = function(fault, lon, lat){
  dtopo = okadaTri(fault, lon, lat, slip=rep(1, length(fault)))
  
  return(dtopo)
}

