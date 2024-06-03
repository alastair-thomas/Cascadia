setwd("~/Uni/NTNU/Masters Project/CSZ/R/Code/MultiQuakeCSZ")
load("DR4.RData") # load in the dataset

# This loads in the required data
# All is created using the "makeTMBFiles.R" script
load("Data.RData") # load in the TMB data
load("Fault.RData") # load in the subfault triangulation
load("spdeMesh.RData") # load in the spde Mesh

# Add the negative log likelihood functions from the TMB code 
compile("MultiQuakeCSZ.cpp",
        framework="TMBad")
dyn.load(dynlib("MultiQuakeCSZ"))

B = dim(spdeMesh$loc)[1]
S = length(fault)
E = length(data$subsidence)

# The parameters to optimise
# Guessing are good starting parameters based on previous optimisations
parameters = list(X          = matrix(0.1, nrow=B, ncol=E),
                  logLambda  = -2.1,
                  mu         = 3.6,
                  logKappa   = -4.4,
                  logTau     = 3.7)

## create the objective function using TMB
TMB::config(tmbad.sparse_hessian_compress = 1)
obj = MakeADFun(data,
                parameters,
                random="X",
                DLL="MultiQuakeCSZ",
                hessian=TRUE)

## Optimize the model parameters
optTime = system.time(opt0 <- nlminb(start     = obj$par,
                                     objective = obj$fn,
                                     gradient  = obj$gr,
                                     control   = list(iter.max  = 100000,
                                                      eval.max  = 100000)))[3]

bestPars = tail(obj$env$last.par.best, 4)

print(bestPars)# best parameters found
print(optTime) # print the time it takes
print(opt0)    # display the optimisation object

## Get standard errors via SD report
sdTime = system.time(SD0 <- TMB::sdreport(obj, getJointPrecision=TRUE,
                                          bias.correct = TRUE,
                                          bias.correct.control = list(sd = TRUE)))[3]

print(sdTime) # show time to run SD report
print(SD0)    # display the error results

bestX   = matrix(summary(SD0, select="random")[,1], ncol=12)
bestXsd = matrix(summary(SD0, select="random")[,2], ncol=12)

plotAllX2(spdeMesh, bestX, DR4)
plotAllX3(spdeMesh, bestXsd, DR4)

# Get draws from all the model parameters
nSims = 1000
draws = simulateCSZ(SD0, nSims=nSims)

## extract the draws for each parameter the draws
parnames = c(names(SD0$par.random), names(SD0$par.fixed))

#extract the random effects
xDraws = draws[parnames == 'X',]
xDraws2 = list()
for (e in 1:E){
  print(paste(((e-1)*B + 1), e*B))
  xDraws2[[e]] = xDraws[((e-1)*B + 1):(e*B),]
}

# extract the fixed parameters
logLambdaDraws = draws[parnames == 'logLambda',]
muDraws        = draws[parnames == 'mu',]
logKappaDraws  = draws[parnames == 'logKappa',]
logTauDraws    = draws[parnames == 'logTau',]


# save everything
allDraws = list(xDraws         = xDraws2,
                logLambdaDraws = logLambdaDraws,
                muDraws        = muDraws,
                logKappaDraws  = logKappaDraws,
                logTauDraws    = logTauDraws)


setwd("~/Uni/NTNU/Masters Project/CSZ/R/Code/MultiQuakeCSZ/Results")
save(allDraws,        file="Draws.RData")
save(opt0,            file="Opt0.RData")
save(optTime, sdTime, file="Times.RData")
save(SD0,             file="SD0.RData")
save(bestPars,        file="Pars.RData")
