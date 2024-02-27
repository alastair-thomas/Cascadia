
library(TMB)
# load in the data needed
load("data_fullFault_T1.RData")
load("par_fullFault_T1.RData")

# Add the negative log likelihood functions from the TMB code 
compile("BayesianCSZ.cpp")
dyn.load(dynlib("BayesianCSZ"))

## create the objective function using TMB
obj = MakeADFun(data,
                parameters,
                random="x",
                DLL="BayesianCSZ",
                hessian=TRUE)

## Optimize the model parameters
opt0 = nlminb(start     = obj$par,
              objective = obj$fn,
              gradient  = obj$gr)

# save the results
save(obj, opt0, file="optFull.RData")