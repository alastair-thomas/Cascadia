library(TMB)

setwd("~/Uni/NTNU/Masters Project/CSZ/R/Code/MultiQuakeAnisoCSZ")

# load in the data
load("DR4.RData")
load("spdeMesh.RData")
load("Fault.RData")
load("Data.RData")

# set constants
B = dim(spdeMesh$loc)[1]
S = length(fault)
E = length(data$subsidence)

# parameters to optimise
parameters = list(X          = matrix(0.1, nrow=B, ncol=E),
                  logLambda  = -2.1,
                  mu         = 3.6,
                  logKappa   = -3.8,
                  logTau     = 2.7,
                  logh       = c(0,deg2rad(-12.39776)))

# Add the negative log likelihood functions from the TMB code
compile("MultiQuakeAnisoCSZ.cpp",
        framework="TMBad")
dyn.load(dynlib("MultiQuakeAnisoCSZ"))

## create the objective function using TMB
TMB::config(DLL="MultiQuakeAnisoCSZ",
            tmbad.sparse_hessian_compress = 1)
obj = MakeADFun(data,
                parameters,
                map = list(logh = as.factor(c(0,NA))),
                random  = "X",
                DLL     = "MultiQuakeAnisoCSZ",
                hessian = TRUE)

## Optimize the model parameters
optTime = system.time(opt0 <- nlminb(start     = obj$par,
                                     objective = obj$fn,
                                     gradient  = obj$gr,
                                     control   = list(iter.max  = 100000,
                                                      eval.max  = 100000)))[3]

bestPars = tail(obj$env$last.par.best, 5)
print(bestPars)
print(opt0)
print(optTime)

## Get standard errors via SD report
sdTime = system.time(SD0 <- TMB::sdreport(obj, getJointPrecision=TRUE,
                                          bias.correct = TRUE,
                                          bias.correct.control = list(sd = TRUE)))[3]

print(sdTime)
print(SD0)

# Get draws from all the model parameters
nSims = 1000
draws = simulateCSZ(SD0, nSims=nSims)

## extract the draws for each parameter the draws
parnames = c(names(SD0[['par.random']]), names(SD0[['par.fixed']]))

print(all(parnames == rownames(SD0$jointOrecision)))

# first randomm effects
tempxDraws = draws[parnames == 'X',]
xDraws = list()
for (e in 1:E){
  print(paste(((e-1)*B + 1), e*B))
  xDraws[[e]] = tempxDraws[((e-1)*B + 1):(e*B),]
}

# now fixed parameters
logLambdaDraws  = draws[parnames == 'logLambda',]
muDraws = draws[parnames == 'mu',]
logKappaDraws = draws[parnames == 'logKappa',]
logTauDraws = draws[parnames == 'logTau',]
logPsiDraws = draws[parnames == 'logh',]

# save everything
results = list(logLambdaDraws = logLambdaDraws,
               muDraws        = muDraws,
               logKappaDraws  = logKappaDraws,
               logTauDraws    = logTauDraws,
               logPsiDraws    = logPsiDraws,
               xDraws         = xDraws)

setwd("~/Uni/NTNU/Masters Project/CSZ/R/Code/MultiQuakeAnisoCSZ/Results")
save(results, file="Draws.RData")
save(optTime, sdTime, file="Times.RData")
save(bestPars, file="Best Params.RData")
save(opt0, file="Opt0.RData")
save(SD0, file="SD0.RData")
