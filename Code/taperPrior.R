## Prior on the taper function parameter lambda

depthLower = 20
depthUpper = 30
alpha = 0.95 # with what confidence
tMin = 0.05 # The value taper should be lower than

# calculate the critical values
lowerCV = -log(tMin) / depthUpper
upperCV = -log(tMin) / depthLower

# Assume λ ~ Gamma(k, theta)
f = function(x){
  k = x[1]
  theta = x[2]
  
  Q1 = qgamma(.025, shape=k, scale=theta)
  Q2 = qgamma(.975, shape=k, scale=theta)
    
  loss = (Q1 - lowerCV)^2 + (Q2 - upperCV)^2
  
  return(loss)
}

testPrior = function(shape, scale){
  nSims = 100000
  depths = seq(0, 30, 0.1)
  testLambda = rgamma(n=nSims, shape=shape, scale=scale)
  count = 0
  for (i in 1:nSims){
    taper = exp(-testLambda[i] * depths)
    idx = which(taper <= tMin)[1]
    
    # if the taper actually got to 0.05
    if (is.na(idx) == FALSE){
      d = depths[idx]
      # if depth of 0.05 value is in in 20km<=d<=30km
      if ((depthLower <= d)&(d <= depthUpper)){
        count = count + 1
      }
    }
  }
  
  testPerc = 100*count/nSims
  return(testPerc)
}

# find the best parameters
priorPars = optim(c(1,1), f, method="Nelder-Mead")
k = priorPars$par[1]
shape = priorPars$par[2]
print(k)
print(shape)

# test that they work
# should get ~alpha*100
test = testPrior(k, shape)
print(test)

# plot the taper given by the mean from the distribution
lambda = k*shape
print(lambda)
depths = seq(0, 30, 0.01)
taper = exp(-lambda*depths)
plotDF = data.frame(depth = depths,
                    taper = taper)
ggplot(plotDF)+
  geom_line(aes(x=depth, y=taper), linewidth=1) +
  labs(x="Depth (km)", y="Taper")
