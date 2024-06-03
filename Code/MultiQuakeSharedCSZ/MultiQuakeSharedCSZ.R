source("~/Uni/NTNU/Masters Project/CSZ/R/Code/Setup.R")
setwd("~/Uni/NTNU/Masters Project/CSZ/R/Code/MultiQuakeSharedCSZ")

# load in the subsidence data
load("DR4.RData")
load("spdeMesh.RData")
load("Fault.Rdata")
load("Data.RData")

# set constants
B = dim(spdeMesh$loc)[1]
S = length(fault)
nSims = 1000
E = 12

# Add the negative log likelihood functions from the TMB code 
compile("MultiQuakeSharedCSZ.cpp",
        framework="TMBad")
dyn.load(dynlib("MultiQuakeSharedCSZ"))

# the parameters to optimise
# I am guessing are good starting parameters
parameters = list(X          = matrix(0.1, nrow=B, ncol=E),
                  w          = rep(0.1, B),
                  logLambda  = -2.1,
                  mu         = 3.6,
                  logKappaX  = -4,
                  logTauX    = 5.5,
                  logKappaW  = -4.2,
                  logTauW    = 3.4)


## create the objective function using TMB
TMB::config(tmbad.sparse_hessian_compress = 1)
obj = MakeADFun(data,
                parameters,
                random=c("X", "w"),
                DLL="MultiQuakeSharedCSZ",
                hessian=TRUE)

## Optimize the model parameters
optTime = system.time(opt0 <- nlminb(start     = obj$par,
                                     objective = obj$fn,
                                     gradient  = obj$gr,
                                     control   = list(iter.max  = 100000,
                                                     eval.max  = 100000)))[3]
# Time to optimise
print(optTime)
print(opt0)

bestPars = tail(obj$env$last.par.best, 6)

# If everything looks good continue and get results from the model.

## Get standard errors via SD report
sdTime = system.time(SD0 <- TMB::sdreport(obj, getJointPrecision=TRUE,
                                          bias.correct = TRUE,
                                          bias.correct.control = list(sd = TRUE)))[3]

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

# extract individual random effects into a list
tempXDraws = draws[parnames == 'X',]
xDraws = list()
for (e in 1:E){
  print(paste(((e-1)*B + 1), e*B))
  xDraws[[e]] = tempXDraws[((e-1)*B + 1):(e*B),]
}

# extract shared random effects
wDraws          = draws[parnames == 'w',]

# extract fixed effects
logLambdaDraws  = draws[parnames == 'logLambda',]
muDraws         = draws[parnames == 'mu',]
logKappaXDraws  = draws[parnames == 'logKappaX',]
logTauXDraws    = draws[parnames == 'logTauX',]
logKappaWDraws  = draws[parnames == 'logKappaW',]
logTauWDraws    = draws[parnames == 'logTauW',]



# save everything
results = list(logLambdaDraws = logLambdaDraws,
               muDraws        = muDraws,
               logKappaXDraws = logKappaXDraws,
               logTauXDraws   = logTauXDraws,
               logKappaWDraws = logKappaWDraws,
               logTauWDraws   = logTauWDraws,
               xDraws         = xDraws,
               wDraws         = wDraws)

setwd("~/Uni/NTNU/Masters Project/CSZ/R/Code/MultiQuakeSharedCSZ/Results")
save(results, file="Draws.RData")
save(optTime, sdTime, file="Times.RData")
save(bestPars, file="Best Params.RData")
save(opt0, file="Opt0.RData")
save(SD0, file="SD0.RData")
