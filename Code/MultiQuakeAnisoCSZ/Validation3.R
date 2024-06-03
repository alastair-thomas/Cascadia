library(TMB)
library(Matrix)

simulateCSZ = function(SD, nSims=1000){
  # take samples from fitted model
  mu = c(SD$par.fixed, SD$par.random)
  L = Matrix::Cholesky(SD[['jointPrecision']], super = T)
  draws = rmvnorm_prec(mu = mu ,
                       chol_prec = L,
                       nSims = nSims)
  return(draws)
}

# simulate draws
rmvnorm_prec = function(mu, chol_prec, nSims) {
  z = matrix(rnorm(length(mu) * nSims), ncol=nSims)
  L = chol_prec #Cholesky(prec, super=TRUE)
  z = Matrix::solve(L, z, system = "Lt") ## z = Lt^-1 %*% z
  z = Matrix::solve(L, z, system = "Pt") ## z = Pt    %*% z
  z = as.matrix(z)
  return(mu + z)
}

# degrees to radians
deg2rad = function(deg){
  return(2*pi*deg / 360)
}

# set correct working directory
setwd("~/Uni/NTNU/Masters Project/CSZ/R/Code/MultiQuakeAnisoCSZ")

# Add the negative log likelihood functions from the TMB code 
compile("MultiQuakeAnisoCSZ.cpp",
        framework="TMBad")

# load in the subsidence data
load("DR4.RData")

# load in the model data
load("backupData.RData")

# load in the fault geometry
load("backupFault.RData")
S = length(fault)

# load in the spdeMesh
load("backupMesh.RData")
B = dim(spdeMesh$loc)[1]

# number of simulations to draw
nSims = 1000

# Number of megathrust events
E = length(data$subsidence)

foldDraws = list()

# do all the validation folds at once
for (fold in 1:E){
  
  # see update
  print(fold)
  
  # Create fold data
  thisData = data
  thisData$okada[[fold]] = NULL
  thisData$subsidence[[fold]] = NULL
  thisData$v[[fold]] = NULL
  
  # parameters to optimise
  thisParameters = list(X          = matrix(0.1, nrow=B, ncol=(E-1)),
                        logLambda  = -2.1,
                        mu         = 3.6,
                        logKappa   = -3.8,
                        logTau     = 2.7,
                        logh       = c(0,deg2rad(-12.39776)))
  
  # TMB stuff
  dyn.load(TMB::dynlib("MultiQuakeAnisoCSZ"))
  #TMB::config(tmbad.sparse_hessian_compress = 1)
  obj = TMB::MakeADFun(thisData,
                       thisParameters,
                       map = list(logh = as.factor(c(0,NA))),
                       random=c("X"),
                       DLL="MultiQuakeAnisoCSZ",
                       hessian=TRUE)
  
  ## Optimize the model parameters
  opt0 = nlminb(start     = obj$par,
                objective = obj$fn,
                gradient  = obj$gr,
                control   = list(iter.max  = 100000,
                                 eval.max  = 100000))
  
  ## Get standard errors via SD report
  SD0 = TMB::sdreport(obj, getJointPrecision=TRUE,
                           bias.correct = TRUE,
                           bias.correct.control = list(sd = TRUE))
  
  # unload TMB code
  dyn.unload(TMB::dynlib("MultiQuakeAnisoCSZ"))
  
  # Get draws from all the model parameters
  thisDraws = simulateCSZ(SD=SD0, nSims=nSims)
  parnames = c(names(SD0[['par.fixed']]), names(SD0[['par.random']]))
  
  
  foldDraws[[fold]] = list(logLambda = thisDraws[parnames == 'logLambda',],
                           mu        = thisDraws[parnames == 'mu',],
                           logKappa  = thisDraws[parnames == 'logKappa',],
                           logTau    = thisDraws[parnames == 'logTau',],
                           logh      = thisDraws[parnames == 'logh',])
  
}

setwd("~/Uni/NTNU/Masters Project/CSZ/R/Code/MultiQuakeAnisoCSZ/Backup Results")
save(foldDraws, file="Backup Validation Draws.RData")