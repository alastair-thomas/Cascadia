load("mixScores.RData")
load("mixFullScores.RData")

df1 = scoreDF
df2 = subset(fullScoreDF, select=-data)

scores = rbind(df1, df2)
scores$model = as.factor(c(rep("Validation", dim(df1)[1]),
                           rep("Full Model", dim(df2)[1])))

scores$Event = as.factor(scores$Event)

ggplot(scores, aes(x=Event,y=CRPS,fill=model)) +
  geom_boxplot(outlier.shape = NA) +
  stat_summary(fun=mean,
               geom='point',
               shape=20,
               size=3,
               aes(colour=model),
               position = position_dodge(width = 0.75)) +
  scale_fill_scico_d(palette="berlin") +
  scale_colour_scico_d(palette="roma", direction = -1) +
  labs(fill="Model", colour="Model Mean")

Nj = unname(table(df1$Event))
events = unique(df1$Event)
eventMeans = data.frame(Event = events,
                        SE = tapply(df1$SE, df1$Event, mean),
                        AE = tapply(df1$AE, df1$Event, mean),
                        CRPS = tapply(df1$CRPS, df1$Event, mean),
                        Coverage = (unname(tapply(df1$Coverage, df1$Event, sum))/Nj)*100,
                        Width = tapply(df1$Width, df1$Event, mean),
                        Interval = tapply(df1$Interval, df1$Event, mean))


# Calculate the scores
# mse, mae, crps, empirical coverage, width
foldErrors = list()
foldSubsidence = okadaSubs
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
}


# make a dataframe of all scores
events = c("T1","T2","T3","T4","T5","T6","T7","T8","T9","T10","T11","T12")
for (fold in 1:E){
  Nj = length(foldErrors[[fold]]$data)
  
  if(fold == 1){
    errorsMat = cbind(Event    = rep(events[fold], Nj),
                      data     = foldErrors[[fold]]$data,
                      SE       = foldErrors[[fold]]$SE,
                      AE       = foldErrors[[fold]]$AE,
                      CRPS     = foldErrors[[fold]]$CRPS,
                      Coverage = foldErrors[[fold]]$Coverage,
                      Width    = foldErrors[[fold]]$Width)
  }
  else{
    errorsMat = rbind(errorsMat,
                      cbind(Event    = rep(events[fold], Nj),
                            data     = foldErrors[[fold]]$data,
                            SE       = foldErrors[[fold]]$SE,
                            AE       = foldErrors[[fold]]$AE,
                            CRPS     = foldErrors[[fold]]$CRPS,
                            Coverage = foldErrors[[fold]]$Coverage,
                            Width    = foldErrors[[fold]]$Width))
  }
}

fullScoreDF = data.frame(errorsMat)
fullScoreDF$Event = factor(fullScoreDF$Event,
                           levels=events,
                           ordered=T)

fullScoreDF$SE = as.numeric(fullScoreDF$SE)
fullScoreDF$AE = as.numeric(fullScoreDF$AE)
fullScoreDF$CRPS = as.numeric(fullScoreDF$CRPS)
fullScoreDF$Coverage = as.numeric(fullScoreDF$Coverage)
fullScoreDF$Width = as.numeric(fullScoreDF$Width)


save(fullScoreDF, file="mixFullScores.RData")

# 
# 
# 
# 
# 
# 
# 
# 
# 
# # Now calculate predictive errors
# # MSE, MAE, CRPS, Empirical Coverage 95% level, 95% credible interval width
# 
# # First collate all the simulated subsidences
# okadaSubs = list()
# for (e in 1:E){
#   thisG = G[[e]]
#   theseSlipDraws = slipDraws[[e]]
#   
#   okadaSubSims = matrix(0, nrow=dim(thisG)[1], ncol=nSims)
#   
#   for (n in 1:nSims){
#     thisSlip = theseSlipDraws[,n]
#     okadaSubSims[,n] = thisG %*% thisSlip
#   }
#   
#   okadaSubs[[e]] = okadaSubSims
# }
# 
# # Calculate the error
# error = list()
# for (e in 1:E){
#   # extract real subsidence
#   thisSub = subsidences[[e]]
#   
#   # average over all the simulations
#   thisOkadaSub = rowMeans(okadaSubs[[e]])
#   
#   # calculate error for each data point
#   error[[e]] = thisSub - thisOkadaSub
# }
# 
# # Plot the errors by earthquake and latitude
# earthquakes = as.factor(c("T1","T2","T3","T4","T5","T6","T7","T8","T9","T10","T11","T12"))
# allErrors = c()
# allQuakes = c()
# allLats   = c()
# for (e in 1:E){
#   allErrors = c(allErrors, error[[e]])
#   allQuakes = c(allQuakes, rep(earthquakes[e], length(error[[e]])))
#   allLats   = c(allLats, DR3$Lat[DR3$event == earthquakes[e]])
# }
# plotDF = data.frame(Quake = as.factor(allQuakes),
#                     Lat   = allLats,
#                     Error = allErrors)
# g11 = ggplot(plotDF) +
#   geom_point(aes(x=Error, y=Lat, colour=Quake)) +
#   scale_colour_viridis_d(option="plasma") +
#   xlim(c(-2, 2))
# plot(g11)
# 
# # Finally calculate the error stats for this model
# squareErrors   = c()
# absoluteErrors = c()
# for (e in 1:E){
#   # extract real subsidence
#   thisSub = subsidences[[e]]
#   
#   # average over all the simulations
#   thisOkadaSub = rowMeans(okadaSubs[[e]])
#   
#   squareErrors   = c(squareErrors, (thisSub - thisOkadaSub)^2)
#   absoluteErrors = c(absoluteErrors, abs(thisSub - thisOkadaSub))
# }
# mse = mean(squareErrors)
# mae = mean(absoluteErrors)


setwd("~/Uni/NTNU/Masters Project/CSZ/R/Code/MultiQuakeCSZ")
load("backupFault.RData")
load("backupMesh.RData")
load("backupData.RData")

setwd("~/Uni/NTNU/Masters Project/CSZ/R/Code/MultiQuakeCSZ/Backup Results")
load("Backup Validation Draws.RData")


# Calculate the subsidence for the left out earthquakes
foldSubsidence = list()
for (fold in 1:E){
  
  theseSubsidences = matrix(0, nrow = dim(data$okada[[fold]])[1],
                            ncol = nSims)
  
  for (n in 1:nSims){
    print(paste(fold, n))
    logLambda = foldDraws[[fold]]$logLambda[n]
    mu        = foldDraws[[fold]]$mu[n]
    logKappa = foldDraws[[fold]]$logKappa[n]
    logTau   = foldDraws[[fold]]$logTau[n]
    
    Qx = inla.spde2.precision(spdeINLA,
                              theta=c(logTau, logKappa))
    X = inla.qsample(Q=Qx)
    
    taper = exp(-exp(logLambda)*data$depth)
    slips = taper * exp(mu + data$A %*% X)
    theseSubsidences[,n] = -(data$okada[[fold]] %*% slips)[,1]
  }
  foldSubsidence[[fold]] = theseSubsidences
}

save(foldSubsidence, file="mixSubsidencesValidation.RData")

# Calculate the scores
# mse, mae, crps, empirical coverage, width
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


save(scoreDF, file="Validation Scores Model 1.RData")
write.csv(scoreDF, "Validation Scores Model 1.csv")

