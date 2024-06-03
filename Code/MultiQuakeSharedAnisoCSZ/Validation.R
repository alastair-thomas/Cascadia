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

# set correct working directory
setwd("~/Uni/NTNU/Masters Project/CSZ/R/Code/MultiQuakeSharedAnisoCSZ")

# Add the negative log likelihood functions from the TMB code 
compile("MultiQuakeSharedAnisoCSZ.cpp",
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
  
  thisParameters = list(X          = matrix(0.1, nrow=B, ncol=(E-1)),
                        w          = rep(0.1, B),
                        logLambda  = -2.1,
                        mu         = 3.6,
                        logKappaX  = -4,
                        logTauX    = 5.5,
                        logKappaW  = -4.2,
                        logTauW    = 3.4,
                        logh       = c(0,deg2rad(-12.39776)))
  
  # TMB stuff
  dyn.load(TMB::dynlib("MultiQuakeSharedAnisoCSZ"))
  TMB::config(DLL="MultiQuakeSharedAnisoCSZ",
              tmbad.sparse_hessian_compress = 1)
  obj = TMB::MakeADFun(thisData,
                       thisParameters,
                       map = list(logh = as.factor(c(0,NA))),
                       random=c("X", "w"),
                       DLL="MultiQuakeSharedAnisoCSZ",
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
  dyn.unload(TMB::dynlib("MultiQuakeSharedAnisoCSZ"))
  
  # Get draws from all the model parameters
  thisDraws = simulateCSZ(SD=SD0, nSims=nSims)
  parnames = c(names(SD0[['par.fixed']]), names(SD0[['par.random']]))
  
  
  foldDraws[[fold]] = list(logLambda = thisDraws[parnames == 'logLambda',],
                           mu        = thisDraws[parnames == 'mu',],
                           logKappaX = thisDraws[parnames == 'logKappaX',],
                           logTauX   = thisDraws[parnames == 'logTauX',],
                           logKappaW = thisDraws[parnames == 'logKappaW',],
                           logTauW   = thisDraws[parnames == 'logTauW',],
                           logh      = thisDraws[parnames == 'logh',])
  
}

setwd("~/Uni/NTNU/Masters Project/CSZ/R/Code/MultiQuakeSharedAnisoCSZ/Backup Results")
save(foldDraws, file="Backup Validation Draws.RData")