library(TMB)
compile("linreg_parallel.cpp")
dyn.load(dynlib("linreg_parallel"))
openmp(max=TRUE) ## Set max available threads

set.seed(123)
x <- seq(0, 10, length=50001)
data <- list(Y=rnorm(length(x)) + x, x=x)
parameters <- list(a=0, b=0, logSigma=0)
obj <- MakeADFun(data, parameters, DLL="linreg_parallel")
opt <- do.call("optim", obj)
