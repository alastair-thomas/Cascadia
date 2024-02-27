## Prior on the taper function parameter lambda

depthLower = 20
depthUpper = 30
alpha = 0.95 # with what confidence
tMin = 0.05 # The value taper should be lower than

# calculate the critical values
lowerCV = -log(tMin) / depthUpper
upperCV = -log(tMin) / depthLower

# Assume λ ~ Gamma(k, theta)
# The function we actually want the roots of
f = function(x){
  k = x[1]
  theta = x[2]
  
  Q2 = qgamma(.025, shape=k, scale=theta)
  Q1 = qgamma(.975, shape=k, scale=theta)
  
  lowerVal = -log(.05)/30
  upperVal = -log(.05)/20
    
  loss = (Q2 - lowerVal)^2 + (Q1 - upperVal)^2
  
  #f1 = (1 - pgamma(lowerCV, shape=k, scale=theta)) - ((1-alpha)/2)
  #f2 = pgamma(upperCV, shape=k, scale=theta) - (alpha + (1-alpha)/2)
  return(loss)
  #f1 = Q2 - lowerVal
  #f2 = Q1 - upperVal
  #return(c(f1, f2))
}

testPrior = optim(c(1,1), f, method="Nelder-Mead")

# The find the root so that the shape which satisfies (1) is found
# This doesn't work currently
#testRoot = multiroot(f, start=c(1,1), maxiter=10000, positive=TRUE)

# function to help manually find good starting conditions
testParams = function(x){
  # Plot the PDF
  testLambda = seq(0,1,0.001)[-1]
  plot(x=testLambda, y=dgamma(testLambda, shape=x[1], scale=x[2]))
  
  # Print useful outputs
  print((1 - pgamma(lowerCV, shape=x[1], scale=x[2])))
  print(pgamma(upperCV, shape=x[1], scale=x[2]))
  print(f(x))
}

plotParamsError = function(){
  x = seq(0.05,0.15,0.001)[-1]
  y = seq(0.05,0.15,0.001)[-1]
  XY = expand.grid(x,y)
  
  z = rep(0, dim(XY)[1])
  for (i in 1:dim(XY)[1]){
    z[i] = abs(f(c(XY$Var1[i], XY$Var2[i]))[1]) +
      abs(f(c(XY$Var1[i], XY$Var2[i]))[2])
  }
  
  data = data.frame(x=XY$Var1, y=XY$Var2, Error=z)
  
  g = ggplot(data, aes(x=x, y=y, fill=Error))+
    geom_tile() +
    scale_fill_viridis_c() +
    labs(x="Shape", y="Scale")
  plot(g)
}


# The find the root so that the shape which satisfies (1) is found
# This doesn't work currently
root = multiroot(f, start=c(0.1,0.1), maxiter=10000, positive=TRUE)

# extract the result
taperPrior = root$root
print(taperPrior)

testPrior = function(shape, scale){
  nSims = 10000
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
  return(abs(testPerc - 95))
}

plotTesting = function(){
  x = seq(0,5,0.1)[-1]
  y = seq(0,5,0.1)[-1]
  XY = expand.grid(x,y)
  
  z = rep(0, dim(XY)[1])
  for (i in 1:dim(XY)[1]){
    z[i] = testPrior(XY$Var1, XY$Var2)
  }
  
  data = data.frame(shape=XY$Var1, scale=XY$Var2, Error=z)
  
  g = ggplot(data, aes(x=shape, y=scale, fill=Error))+
    geom_tile() +
    scale_fill_viridis_c() +
    labs(x="Shape", y="Scale")
  plot(g)
}
