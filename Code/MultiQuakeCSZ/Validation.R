# This script runs the 12-fold earthquake predicting validation
# It is designed to be self sufficient so can be run on faster machines
# Need the data.Rdata file which can be made from the makeInput.R script

# set working directory
setwd("~/Uni/NTNU/Masters Project/CSZ/R/Code/MultiQuakeCSZ")

# load libraries needed
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

# Add the negative log likelihood functions from the TMB code 
TMB::compile("MultiQuakeCSZ.cpp",
             framework="TMBad")

# load in the subsidence data
load("DR4.RData")

# load in the model data
load("Data.RData")

# load in the fault geometry
load("Fault.RData")
S = length(fault)

# load in the spdeMesh
load("Mesh.RData")
B = dim(spdeMesh$loc)[1]

# number of simulations to draw
nSims = 1000

# Number of megathrust events
E = length(data$subsidence)

# the results to collect
foldDraws = list()
for (fold in 1:E) {
  print(fold)
  
  # Create fold data
  thisData = data
  thisData$okada[[fold]] = NULL
  thisData$subsidence[[fold]] = NULL
  thisData$v[[fold]] = NULL
  
  # created fold data
  thisParameters = list(X          = matrix(0.0, nrow=B, ncol=(E-1)),
                        logLambda  = -2.1,
                        mu         = 3.6,
                        logKappa   = -4.4,
                        logTau     = 3.7)
  
  # load the tmb code
  dyn.load(TMB::dynlib("MultiQuakeCSZ"))
  
  # use a memory efficient matrix compression (slower)
  TMB::config(tmbad.sparse_hessian_compress = 1)
  
  # make the objective function
  obj = TMB::MakeADFun(thisData,
                       thisParameters,
                       random="X",
                       DLL="MultiQuakeCSZ",
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
  dyn.unload(TMB::dynlib("MultiQuakeCSZ"))
  
  # Get draws from all the model parameters
  thisDraws = simulateCSZ(SD=SD0, nSims=nSims)
  parnames = c(names(SD0[['par.fixed']]), names(SD0[['par.random']]))
  
  # save results to a list
  foldDraws[[fold]] = list(logLambda = thisDraws[parnames == 'logLambda',],
                           mu        = thisDraws[parnames == 'mu',],
                           logKappa  = thisDraws[parnames == 'logKappa',],
                           logTau    = thisDraws[parnames == 'logTau',])
}

# save all results
setwd("~/Uni/NTNU/Masters Project/CSZ/R/Code/MultiQuakeCSZ/Results")
save(foldDraws, file="Validation Draws.RData")
