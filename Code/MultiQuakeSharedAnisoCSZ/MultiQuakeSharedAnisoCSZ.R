setwd("~/Uni/NTNU/Masters Project/CSZ/R/Code/MultiQuakeSharedAnisoCSZ")

# Building the CSZ model using TMB
library(TMB)

# load in the data
load("DR4.RData")
load("spdeMesh.RData")
load("Fault.RData")
load("Data.RData")

# set constants
B = dim(spdeMesh$loc)[1]
S = length(fault)
E = length(data$subsidence)

# Add the negative log likelihood functions from the TMB code 
compile("MultiQuakeSharedAnisoCSZ.cpp",
        framework="TMBad")
dyn.load(dynlib("MultiQuakeSharedAnisoCSZ"))

# the parameters to optimise
# I am guessing are good starting parameters
parameters = list(X          = matrix(0.1, nrow=B, ncol=E),
                  w          = rep(0.1, B),
                  logLambda  = -2.1,
                  mu         = 3.6,
                  logKappaX  = -4,
                  logTauX    = 5.5,
                  logKappaW  = -4.2,
                  logTauW    = 3.4,
                  logh       = c(0,deg2rad(-12.39776)))


## create the objective function using TMB
TMB::config(DLL = "MultiQuakeSharedAnisoCSZ",
            tmbad.sparse_hessian_compress = 1)
obj = MakeADFun(data,
                parameters,
                map = list(logh = as.factor(c(0,NA))),
                random=c("X", "w"),
                DLL="MultiQuakeSharedAnisoCSZ",
                hessian=TRUE)

## Optimize the model parameters
optTime = system.time(opt0 <- nlminb(start     = obj$par,
                                     objective = obj$fn,
                                     gradient  = obj$gr,
                                     control   = list(iter.max  = 100000,
                                                     eval.max  = 100000)))[3]

# outputs of optimisation
print(opt0)
print(optTime)

# best parameters
bestPars = tail(obj$env$last.par.best, 7)
print(bestPars)

## Get standard errors via SD report
sdTime = system.time(SD0 <- TMB::sdreport(obj, getJointPrecision=TRUE,
                                          bias.correct = TRUE,
                                          bias.correct.control = list(sd = TRUE)))[3]
dyn.unload(dynlib("MultiQuakeSharedAnisoCSZ"))

# time to do sdreport
print(sdTime)
# summarise the SD report
print(SD0)

nSims = 1000
# Get draws from all the model parameters
draws = simulateCSZ(SD0, nSims=nSims)

## extract the draws for each parameter the draws
parnames = c(names(SD0[['par.random']]), names(SD0[['par.fixed']]))

print(all(parnames == rownames(SD0$jointPrecision)))

# random effects
tempxDraws         = draws[parnames == 'X',]
xDraws = list()
for (e in 1:E){
  print(paste(((e-1)*B + 1), e*B))
  xDraws[[e]] = tempxDraws[((e-1)*B + 1):(e*B),]
}
wDraws         = draws[parnames == 'w',]

# fixed effects
logLambdaDraws = draws[parnames == 'logLambda',]
muDraws        = draws[parnames == 'mu',]
logKappaXDraws = draws[parnames == 'logKappaX',]
logTauXDraws   = draws[parnames == 'logTauX',]
logKappaWDraws = draws[parnames == 'logKappaW',]
logTauWDraws   = draws[parnames == 'logTauW',]
logPsiDraws      = draws[parnames == 'logh',]

# save everything
results = list(logLambdaDraws = logLambdaDraws,
               muDraws        = muDraws,
               logKappaXDraws = logKappaXDraws,
               logTauXDraws   = logTauXDraws,
               logKappaWDraws = logKappaWDraws,
               logTauWDraws   = logTauWDraws,
               logPsiDraws    = logPsiDraws,
               xDraws         = xDraws,
               wDraws         = wDraws)

setwd("~/Uni/NTNU/Masters Project/CSZ/R/Code/MultiQuakeSharedAnisoCSZ/Results")
save(results, file="Draws.RData")
save(optTime, sdTime, file="Times.RData")
save(bestPars, file="Best Params.RData")
save(opt0, file="Opt0.RData")
save(SD0, file="SD0.RData")

