
# The depth at which 100Alpha% of the time the taper is less than tMin
depthDrop = 25
alpha = 0.95
tMin = 0.05

# calculate the critical value
cv = -1/depthDrop * log(tMin)

# a function that depends on the shape parameter of gamma prior
# fix the rate to be 1
# return error between CDF at critical value and 1-alpha
# can then use in a minimisation
fn = function(shape){
  return(abs(pgamma(cv, shape=shape, scale=1) - (1-alpha)))
}

# find the best shape value for prior on lambda
# The condition is this:
# P(taper < 0.05 | depth = 25km) = 0.95
# Assume lambda has a gamma prior
#
opt1 = optim(par=c(2), fn)

# extract the result
shapePrior = opt1$par[1]
print(shapePrior)

# Now test that 95% of the time the taper is below 0.05 at depth=25km.
nSims = 1000
depths = seq(0, 30, 0.1)
testLambda = rgamma(n=nSims, shape=shapePrior, rate=1)
count = 0

for (i in 1:nSims){
  taper = exp(-testLambda[i] * depths)
  idx = which(taper <= 0.05)[1]
  
  # if the taper actually got to 0.05
  if (is.na(idx) == FALSE){
    d = depths[idx]
    if (d <= 25){
      count = count + 1
    }
  }
}

testPerc = 100*count/nSims
print(testPerc)
