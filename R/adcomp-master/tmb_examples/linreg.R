library(TMB)
compile("linreg.cpp")
dyn.load(dynlib("linreg"))
set.seed(123)

## creates data Y = a + b*x
createData <- function(a, b, n){
  x <- seq(1, 100, length.out=n)
  Y <- a + b*x
  Y <- Y + rnorm(n)
  
  return(list(x=x, Y=Y))
}

#data <- list(Y = rnorm(10) + 1:10, x=1:10)
data <- createData(a=6, b=4, n=100)
parameters <- list(a=0, b=0, logSigma=0)
obj <- MakeADFun(data, parameters, DLL="linreg")
obj$hessian <- TRUE
opt <- do.call("optim", obj)
opt
opt$hessian ## <-- FD hessian from optim
obj$he()    ## <-- Analytical hessian
sdreport(obj)
