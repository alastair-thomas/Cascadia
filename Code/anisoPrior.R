lowerPsi = 0.5
upperPsi = 2

f = function(v){
  Q1 = qnorm(.025, mean=0, sd=v)
  Q2 = qnorm(.975, mean=0, sd=v)
  
  loss = (Q1 - lowerPsi)^2 + (Q2 - upperPsi)^2
}

# find the best parameters
priorPars = optimize(f, interval=c(0, 1))
print(priorPars)
# results is 0.3826601

testPrior = function(v){
  nSims = 100000
  testLogPhi = rnorm(n=nSims, mean=0, sd=v)
  testPhi = exp(testLogPhi)
  
  count = 0
  for (i in 1:nSims){
    # if the taper actually got to 0.05
    if ((testPhi[i] > lowerPsi)&(testPhi[i] < upperPsi)){
      count = count + 1
    }
  }
  
  testPerc = 100*count/nSims
  return(testPerc)
}

test = testPrior(0.3826601)
print(test)




testPriorStretch = function(v){
  nSims = 100
  u = c(0,0)
  v1 = as.numeric(v[1])
  v2 = as.numeric(v[2])
  S = matrix(c(v1,v2,
               v2,v1), nrow=2)
  
  print(S)
  
  print(matrixcalc::is.positive.semi.definite(S))
  
  logh = rmvnorm(nSims, mu=u, Sigma=S)
  
  h0 = exp(logh[,1])
  h1 = exp(logh[,2])
  
  countStretch = 0
  countAngle = 0
  
  for (i in 1:nSims){
    H = matrix(c(h0[i],      log(h1[i]),
                 log(h1[i]), (1+log(h1[i])^2)/h0[i]),
               nrow=2)
    # Compute eigenvalues and eigenvectors
    eigen_result = eigen(H)
    eigenvalues = eigen_result$values
    eigenvectors = eigen_result$vectors
    
    stretch = sqrt(eigenvalues[1]/eigenvalues[2])
    angle = rad2deg(atan(eigenvectors[1,1]/eigenvectors[2,1]))
    
    if (stretch<5){
      countStretch = countStretch + 1
    }
    
    if ((angle<23.08573)&(angle>-44.6502)){
      countAngle = countAngle + 1
    }
    
  }
  
  percStretch = 100*countStretch/nSims
  percAngle = 100*countAngle/nSims
  
  return(c(percStretch, percAngle))
}

v0 = seq(0.3, 5, 0.1)
v1 = seq(0.3, 5, 0.1)
testGrid = expand.grid(v0, v1)
nTests = dim(testGrid)[1]

testPAngle = c()
testPStretch = c()
testV1 = c()
testV2 = c()

for (i in 1:nTests){
  if (testGrid[i,][1] > testGrid[i,][2]){
    testV1 = c(testV1, testGrid[i,][1])
    testV2 = c(testV2, testGrid[i,][2])
    
    res = testPriorStretch(testGrid[i,])
    
    testPStretch = c(testPStretch, res[1])
    testPAngle   = c(testPAngle, res[2])
  }
}

results = data.frame(v1 = unlist(testV1),
                     v2 = unlist(testV2),
                     PA = unlist(testPAngle),
                     PS = unlist(testPStretch))

ggplot(results)+
  geom_raster(aes(x=v1, y=v2, fill=PS)) +
  scale_fill_viridis()

ggplot(results)+
  geom_raster(aes(x=v1, y=v2, fill=PA)) +
  scale_fill_viridis()

print(results[((results$PStretch<95.5) & (results$PStretch>94.5)),])






testPsiPrior = function(v){
  nSims = 1000
  
  logPsi = rnorm(nSims, mean=0, sd=v)
  
  psi = exp(logPsi)
  
  countStretch = 0
  
  for (i in 1:nSims){
    if ((0.5 < psi[i]) & (psi[i] < 2)){
      countStretch = countStretch + 1
    }
  }
  
  percStretch = 100*countStretch/nSims
  
  return(percStretch)
}

v = seq(0.1,5,0.001)
P = rep(0, length(v))
for (i in 1:length(v)){
  P[i] = testPsiPrior(v[i])
}
results = data.frame(P=P, v=v)
ggplot(plotDF)+
  geom_point(aes(x=v, y=P))

print(results[results$P < 95.5 & results$P > 94.5, ])
# result is 0.36