load("backupMesh.RData")
load("backupData.RData")
load("backupFault.RData")
load("Backup Validation Draws.RData")

B = dim(spdeMesh$loc)[1]    # number of basis functions
S = length(fault)           # number of subfaults
E = length(data$subsidence) # number of earthquakes

# Calculate the subsidence for the left out earthquakes
foldSubsidence = list()
for (fold in 1:E){
  
  theseSubsidences = matrix(0, nrow=dim(data$okada[[fold]])[1],
                               ncol=nSims)
  
  for (n in 1:nSims){
    print(paste(fold, n))
    logLambda = foldDraws[[fold]]$logLambda[n]
    mu        = foldDraws[[fold]]$mu[n]
    logKappaX = foldDraws[[fold]]$logKappaX[n]
    logTauX   = foldDraws[[fold]]$logTauX[n]
    logKappaW = foldDraws[[fold]]$logKappaW[n]
    logTauW   = foldDraws[[fold]]$logTauW[n]
    Qx = inla.spde2.precision(spdeINLA,
                              theta=c(logTauX, logKappaX))
    Qw = inla.spde2.precision(spdeINLA,
                              theta=c(logTauW, logKappaW))
    X = inla.qsample(Q=Qx)
    w = inla.qsample(Q=Qw)
    
    taper = exp(-exp(logLambda)*data$depth)
    slips = taper * exp(mu + data$A %*% X + data$A %*% w)
    theseSubsidences[,n] = -(data$okada[[fold]] %*% slips)[,1]
  }
  foldSubsidence[[fold]] = theseSubsidences
}


# Calculate the scores
# mse, mae, crps, empirical coverage, width
# could be worth running this loop in parallel
foldErrors = list()
for (fold in 1:E){
  foldErrors[[fold]] = list()
  
  yObs  = -(data$subsidence[[fold]])
  foldErrors[[fold]]$data = yObs
  
  yPred = foldSubsidence[[fold]]
  
  yObsMat = matrix(rep(yObs, nSims), ncol=nSims)
  
  allSE = (yObsMat - yPred)^2
  allAE = abs(yObsMat - yPred)
  
  foldErrors[[fold]]$SE = apply(allSE, 1, mean)
  foldErrors[[fold]]$AE = apply(allAE, 1, mean)
  
  foldErrors[[fold]]$CRPS = scoringRules::crps_sample(yObs, yPred)
  
  predQuant = cbind(lower = apply(yPred, 1, quantile, 0.025),
                    upper = apply(yPred, 1, quantile, 0.975))
  
  foldErrors[[fold]]$Coverage = as.numeric(((predQuant[,1] < yObs) & (yObs < predQuant[,2])))
  foldErrors[[fold]]$Width = predQuant[,2] - predQuant[,1]
  foldErrors[[fold]]$Interval = scoringutils::interval_score(yObs,
                                                             lower=predQuant[,1],
                                                             upper=predQuant[,2],
                                                             interval_range=95)
}


# make a dataframe of all scores
events = c("T1","T2","T3","T4","T5","T6","T7","T8","T9","T10","T11","T12")
for (fold in 1:E){
  Nj = length(data$subsidence[[fold]])
  
  totTime = foldOptTimes[fold] + foldSDTimes[fold]
  
  if(fold == 1){
    errorsMat = cbind(Event    = rep(events[fold], Nj),
                      SE       = foldSE[[fold]],
                      AE       = foldAE[[fold]],
                      CRPS     = foldCRPS[[fold]],
                      Coverage = foldCoverage[[fold]],
                      Width    = foldWidth[[fold]],
                      Time     = rep(totTime, Nj))
  }
  else{
    errorsMat = rbind(errorsMat,
                      cbind(Event    = rep(events[fold], Nj),
                            SE       = foldSE[[fold]],
                            AE       = foldAE[[fold]],
                            CRPS     = foldCRPS[[fold]],
                            Coverage = foldCoverage[[fold]],
                            Width    = foldWidth[[fold]],
                            Time     = rep(totTime, Nj)))
  }
}
scoreDF = data.frame(errorsMat)
scoreDF$Event = factor(scoreDF$Event,
                       levels=events,
                       ordered=T)
scoreDF$SE = as.numeric(scoreDF$SE)
scoreDF$AE = as.numeric(scoreDF$AE)
scoreDF$CRPS = as.numeric(scoreDF$CRPS)
scoreDF$Coverage = as.numeric(scoreDF$Coverage)
scoreDF$Width = as.numeric(scoreDF$Width)
scoreDF$Time = as.numeric(scoreDF$Time)

# Now plot some results
ggplot(scoreDF)+
  geom_boxplot(aes(x=Event, y=CRPS, fill=Time))+
  scale_fill_viridis_c()+
  labs(fill="Time (s)")

save(scoreDF, file="Validation Scores Model 2.RData")
write.csv(scoreDF, "Validation Scores Model 2.csv")

