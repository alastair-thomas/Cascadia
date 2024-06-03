
getScoresFromSubsidence = function(foldSubsidence, data, nSims=1000){
  
  E = 12
  
  # Calculate the scores
  # mse, mae, crps, empirical coverage, width
  foldErrors = list()
  for (fold in 1:E){
    foldErrors[[fold]] = list()
    
    yObs  = -(data$subsidence[[fold]])
    foldErrors[[fold]]$estimate = yObs
    
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
    Nj = length(foldErrors[[fold]]$estimate)
    
    if(fold == 1){
      errorsMat = cbind(Event    = rep(events[fold], Nj),
                        data     = foldErrors[[fold]]$estimate,
                        SE       = foldErrors[[fold]]$SE,
                        AE       = foldErrors[[fold]]$AE,
                        CRPS     = foldErrors[[fold]]$CRPS,
                        Coverage = foldErrors[[fold]]$Coverage,
                        Width    = foldErrors[[fold]]$Width,
                        Interval = foldErrors[[fold]]$Interval)
    }
    else{
      errorsMat = rbind(errorsMat,
                        cbind(Event    = rep(events[fold], Nj),
                              data     = foldErrors[[fold]]$estimate,
                              SE       = foldErrors[[fold]]$SE,
                              AE       = foldErrors[[fold]]$AE,
                              CRPS     = foldErrors[[fold]]$CRPS,
                              Coverage = foldErrors[[fold]]$Coverage,
                              Width    = foldErrors[[fold]]$Width,
                              Interval = foldErrors[[fold]]$Interval))
    }
  }
  
  scoreDF = data.frame(errorsMat)
  
  scoreDF$Event = factor(scoreDF$Event,
                         levels=events,
                         ordered=T)
  scoreDF$SE       = as.numeric(scoreDF$SE)
  scoreDF$AE       = as.numeric(scoreDF$AE)
  scoreDF$CRPS     = as.numeric(scoreDF$CRPS)
  scoreDF$Coverage = as.numeric(scoreDF$Coverage)
  scoreDF$Width    = as.numeric(scoreDF$Width)
  scoreDF$Interval = as.numeric(scoreDF$Interval)
  
  return(scoreDF)
}

#---------------------------MultiQuakeCSZ---------------------------------------

# load data required
setwd("~/Uni/NTNU/Masters Project/CSZ/R/Code/MultiQuakeCSZ")
load("Fault.RData")
load("spdeMesh.RData")
load("Data.RData")

# load validation Draws
setwd("~/Uni/NTNU/Masters Project/CSZ/R/Code/MultiQuakeCSZ/Results")
load("Validation Draws.RData")

# set constants
E = length(data$subsidence)
S = length(fault)
B = dim(spdeMesh$loc)[1]
nSims = 1000

# Add the TMB code to make the Q matrix
setwd("~/Uni/NTNU/Masters Project/CSZ/R/Code/MultiQuakeCSZ")
compile("MakeQ.cpp")

# data object needed to make Q
qData = list(M0 = data$M0,
             M1 = data$M1,
             M2 = data$M2)

# Calculate the subsidence for the left out earthquakes
foldSubsidence1 = list()
for (fold in 1:E){
  
  theseSubsidences = matrix(0, nrow = dim(data$okada[[fold]])[1],
                               ncol = nSims)
  
  for (n in 1:nSims){
    print(paste(fold, n))
    
    logLambda = foldDraws1[[fold]]$logLambda[n]
    mu        = foldDraws1[[fold]]$mu[n]
    logKappa  = foldDraws1[[fold]]$logKappa[n]
    logTau    = foldDraws1[[fold]]$logTau[n]
    
    # Parameters to make the Q matrix with
    qParameters = list(logKappa = logKappa,
                       logTau   = logTau)
    # Make TMB object
    dyn.load(dynlib("MakeQ"))
    obj = MakeADFun(qData,
                    qParameters,
                    DLL = "MakeQ")
    
    # extract Q matrix
    rep = obj$report()
    Qx  = as.matrix(rep$Q)
    dyn.unload(dynlib("MakeQ"))
    
    X = inla.qsample(n=1, Q=Qx)
    
    taper = exp(-exp(logLambda)*data$depth)
    slips = taper * exp(mu + data$A %*% X)
    theseSubsidences[,n] = -(data$okada[[fold]] %*% slips)[,1]
  }
  foldSubsidence1[[fold]] = theseSubsidences
}

# save fold subsidences since takes so long to make
setwd("~/Uni/NTNU/Masters Project/CSZ/R/Code/MultiQuakeCSZ/Backup Results")
save(foldSubsidence1, file="Backup Validation Subsidences.RData")

# get the validation model scores
validScoreDF1 = getScoresFromSubsidence(foldSubsidence1, data)
setwd("~/Uni/NTNU/Masters Project/CSZ/R/Code/MultiQuakeCSZ/Backup Results")
save(validScoreDF1, file="Backup Validation Scores.RData")

# get the scores for the full model
load("Backup Subsidence.RData")
fullScoreDF1 = getScoresFromSubsidence(okadaSubs, data)

plotInnerOuterError(fullScoreDF1, validScoreDF1)
#-------------------------------------------------------------------------------
#---------------------------MultiQuakeSharedCSZ---------------------------------

# load data required
setwd("~/Uni/NTNU/Masters Project/CSZ/R/Code/MultiQuakeSharedCSZ")
load("backupFault.RData")
load("backupMesh.RData")
load("backupData.RData")

# load validation draws
setwd("~/Uni/NTNU/Masters Project/CSZ/R/Code/MultiQuakeSharedCSZ/Backup Results")
load("Backup Validation Params.RData")

# set constants
E = length(data$subsidence)
S = length(fault)
nSims = 1000

# Add the TMB code to make the Q matrices
setwd("~/Uni/NTNU/Masters Project/CSZ/R/Code/MultiQuakeSharedCSZ")
compile("MakeQShared.cpp")

# data object needed to make Q
qData = list(M0 = data$M0,
             M1 = data$M1,
             M2 = data$M2)

# Calculate the subsidence for the left out earthquakes
foldSubsidence2 = list()

for (fold in 1:E){
  
  theseSubsidences = matrix(0, nrow = dim(data$okada[[fold]])[1],
                               ncol = nSims)
  
  # find the best pars
  logLambda = foldParams[[fold]]$logLambda
  mu        = foldParams[[fold]]$mu
  logKappaX = foldParams[[fold]]$logKappaX
  logTauX   = foldParams[[fold]]$logTauX
  logKappaW = foldParams[[fold]]$logKappaW
  logTauW   = foldParams[[fold]]$logTauW
  
  qParameters = list(logKappaX = logKappaX,
                     logTauX   = logTauX,
                     logKappaW = logKappaW,
                     logTauW   = logTauW)
  
  # Make TMB object
  dyn.load(dynlib("MakeQShared"))
  obj = MakeADFun(qData,
                  qParameters,
                  DLL = "MakeQShared")
  
  # extract Q matrix
  rep = obj$report()
  Qx  = as.matrix(rep$Qx)
  Qw  = as.matrix(rep$Qw)
  dyn.unload(dynlib("MakeQShared"))
  
  allX = inla.qsample(n=nSims, Q=Qx)
  allw = inla.qsample(n=nSims, Q=Qw)
  
  for (n in 1:nSims){
    print(paste(fold, n))
    
    X = allX[,n]
    w = allw[,n]
    
    taper = exp(-exp(logLambda)*data$depth)
    slips = taper * exp(mu + data$A %*% X + data$A %*% w)
    theseSubsidences[,n] = -(data$okada[[fold]] %*% slips)[,1]
  }
  foldSubsidence2[[fold]] = theseSubsidences
}

# save fold subsidences since takes so long to make
setwd("~/Uni/NTNU/Masters Project/CSZ/R/Code/MultiQuakeSharedCSZ/Backup Results")
save(foldSubsidence2, file="Backup Validation Subsidences.RData")

# get the validation model scores
validScoreDF2 = getScoresFromSubsidence(foldSubsidence2, data)
setwd("~/Uni/NTNU/Masters Project/CSZ/R/Code/MultiQuakeSharedCSZ/Backup Results")
save(validScoreDF2, file="Backup Validation Scores.RData")

# get the scores for the full model
setwd("~/Uni/NTNU/Masters Project/CSZ/R/Code/MultiQuakeSharedCSZ/Backup Results")
load("Backup Subsidence.RData")
fullScoreDF2 = getScoresFromSubsidence(okadaSubs, data)

plotInnerOuterError(fullScoreDF2, validScoreDF2)
#-------------------------------------------------------------------------------
#---------------------------MultiQuakeAnisoCSZ----------------------------------

# load data required
setwd("~/Uni/NTNU/Masters Project/CSZ/R/Code/MultiQuakeAnisoCSZ")
load("backupFault.RData")
load("backupMesh.RData")
load("backupData.RData")

# load validation draws
setwd("~/Uni/NTNU/Masters Project/CSZ/R/Code/MultiQuakeAnisoCSZ/Backup Results")
load("Backup Validation Params.RData")

# set constants
E = length(data$subsidence)
S = length(fault)
nSims = 1000

# Add the TMB code to make the Q matrix
setwd("~/Uni/NTNU/Masters Project/CSZ/R/Code/MultiQuakeAnisoCSZ")
compile("MakeQAniso.cpp")


# data object needed to make Q
qData = list(spde = data$spde)

# Calculate the subsidence for the left out earthquakes
foldSubsidence3 = list()
for (fold in 1:E){
  
  theseSubsidences = matrix(0, nrow = dim(data$okada[[fold]])[1],
                               ncol = nSims)
  
  logLambda = foldParams[[fold]]$logLambda
  mu        = foldParams[[fold]]$mu
  logKappa  = foldParams[[fold]]$logKappa
  logTau    = foldParams[[fold]]$logTau
  logPsi    = foldParams[[fold]]$logh
  
  qParameters = list(logKappa = logKappa,
                     logTau   = logTau,
                     logh     = c(logPsi, deg2rad(-12.39776)))
  # Make TMB object
  dyn.load(dynlib("MakeQAniso"))
  obj = MakeADFun(qData,
                  qParameters,
                  map = list(logKappa = as.factor(NA),
                             logTau   = as.factor(NA),
                             logh     = as.factor(c(NA,NA))),
                  DLL = "MakeQAniso")
  
  # extract Q matrix
  rep = obj$report()
  Qx  = as.matrix(rep$Q)
  dyn.unload(dynlib("MakeQAniso"))
  
  # simulate
  allX = inla.qsample(n=nSims, Q=Qx)
  
  for (n in 1:nSims){
    print(paste(fold, n))
    
    X = allX[,n]
    
    taper = exp(-exp(logLambda)*data$depth)
    slips = taper*exp(mu + data$A %*% X)
    theseSubsidences[,n] = -(data$okada[[fold]] %*% slips)[,1]
  }
  foldSubsidence3[[fold]] = theseSubsidences
}

# save fold subsidences since takes so long to make
setwd("~/Uni/NTNU/Masters Project/CSZ/R/Code/MultiQuakeAnisoCSZ/Backup Results")
save(foldSubsidence3, file="Backup Validation Subsidences.RData")


# get the validation model scores
validScoreDF3 = getScoresFromSubsidence(foldSubsidence3, data)
setwd("~/Uni/NTNU/Masters Project/CSZ/R/Code/MultiQuakeAnisoCSZ/Backup Results")
save(validScoreDF3, file="Backup Validation Scores.RData")

# get the scores for the full model
setwd("~/Uni/NTNU/Masters Project/CSZ/R/Code/MultiQuakeAnisoCSZ/Backup Results")
load("Backup Subsidence.RData")
fullScoreDF3 = getScoresFromSubsidence(okadaSubs, data)

plotInnerOuterError(fullScoreDF3, validScoreDF3)
#-------------------------------------------------------------------------------
#--------------------------MultiQuakeSharedAnisoCSZ-----------------------------

# load data required
setwd("~/Uni/NTNU/Masters Project/CSZ/R/Code/MultiQuakeSharedAnisoCSZ")
load("backupFault.RData")
load("backupMesh.RData")
load("backupData.RData")

# load validation draws
setwd("~/Uni/NTNU/Masters Project/CSZ/R/Code/MultiQuakeSharedAnisoCSZ/Backup Results")
load("Backup Validation Params.RData")

# set constants
E = length(data$subsidence)
S = length(fault)
nSims = 1000

# Add the TMB code to make the Q matrix
setwd("~/Uni/NTNU/Masters Project/CSZ/R/Code/MultiQuakeSharedAnisoCSZ")
compile("MakeQSharedAniso.cpp")

# data object needed to make Q
qData = list(spde = data$spde)

# Calculate the subsidence for the left out earthquakes
foldSubsidence4 = list()
for (fold in 1:E){
  
  theseSubsidences = matrix(0, nrow = dim(data$okada[[fold]])[1],
                               ncol = nSims)
  
  logLambda = foldParams[[fold]]$logLambda
  mu        = foldParams[[fold]]$mu
  logKappaX = foldParams[[fold]]$logKappaX
  logTauX   = foldParams[[fold]]$logTauX
  logKappaW = foldParams[[fold]]$logKappaW
  logTauW   = foldParams[[fold]]$logTauW
  logPsi    = foldParams[[fold]]$logh
  
  qParameters = list(logKappaX = logKappaX,
                     logTauX   = logTauX,
                     logKappaW = logKappaW,
                     logTauW   = logTauW,
                     logh      = c(logPsi, deg2rad(-12.39776)))
  # Make TMB object
  dyn.load(dynlib("MakeQSharedAniso"))
  obj = MakeADFun(qData,
                  qParameters,
                  map = list(logKappaX = as.factor(NA),
                             logTauX   = as.factor(NA),
                             logKappaW = as.factor(NA),
                             logTauW   = as.factor(NA),
                             logh      = as.factor(c(NA,NA))),
                  DLL = "MakeQSharedAniso")
  # extract Q matrix
  rep = obj$report()
  Qx = as.matrix(rep$Qx)
  Qw = as.matrix(rep$Qw)
  dyn.unload(dynlib("MakeQSharedAniso"))
  
  allX = inla.qsample(n=nSims, Q=Qx)
  allw = inla.qsample(n=nSims, Q=Qw)
  
  for (n in 1:nSims){
    print(paste(fold, n))
    
    X = allX[,n]
    w = allw[,n]
    
    taper = exp(-exp(logLambda)*data$depth)
    slips = taper * exp(mu + data$A %*% X + data$A %*% w)
    theseSubsidences[,n] = -(data$okada[[fold]] %*% slips)[,1]
  }
  foldSubsidence4[[fold]] = theseSubsidences
}

# save fold subsidences since takes so long to make
setwd("~/Uni/NTNU/Masters Project/CSZ/R/Code/MultiQuakeSharedAnisoCSZ/Backup Results")
save(foldSubsidence4, file="Backup Validation Subsidences.RData")

# get the validation model scores
validScoreDF4 = getScoresFromSubsidence(foldSubsidence4, data)
setwd("~/Uni/NTNU/Masters Project/CSZ/R/Code/MultiQuakeSharedAnisoCSZ/Backup Results")
save(validScoreDF4, file="Backup Validation Scores.RData")

# get the scores for the full model
setwd("~/Uni/NTNU/Masters Project/CSZ/R/Code/MultiQuakeSharedAnisoCSZ/Backup Results")
load("Backup Subsidence.RData")
fullScoreDF4 = getScoresFromSubsidence(okadaSubs, data)

plotInnerOuterError(fullScoreDF4, validScoreDF4)
#-------------------------------------------------------------------------------
#---------------------------Subsidence Summary Plots----------------------------

# Plot all the validation subsidences
setwd("~/Uni/NTNU/Masters Project/CSZ/R/Code/MultiQuakeCSZ/Backup Results")
load("Backup Validation Subsidences.RData")

setwd("~/Uni/NTNU/Masters Project/CSZ/R/Code/MultiQuakeSharedCSZ/Backup Results")
load("Backup Validation Subsidences.RData")

setwd("~/Uni/NTNU/Masters Project/CSZ/R/Code/MultiQuakeAnisoCSZ/Backup Results")
load("Backup Validation Subsidences.RData")

setwd("~/Uni/NTNU/Masters Project/CSZ/R/Code/MultiQuakeSharedAnisoCSZ/Backup Results")
load("Backup Validation Subsidences.RData")

setwd("~/Uni/NTNU/Masters Project/CSZ/R/Code/MultiQuakeCSZ")
load("DR4.RData")

E = 12

plotDF = data.frame()
earthquakes = c("T1","T2","T3","T4","T5","T6","T7","T8","T9","T10","T11","T12")
for (e in 1:E){
  Nj = dim(foldSubsidence1[[e]])[1]
  
  # mean subsidence
  meanSub1 = rowMeans(foldSubsidence1[[e]])
  meanSub2 = rowMeans(foldSubsidence2[[e]])
  meanSub3 = rowMeans(foldSubsidence3[[e]])
  meanSub4 = rowMeans(foldSubsidence4[[e]])
  
  thisDF = data.frame(subsidence = c(DR4$subsidence[DR4$event == earthquakes[e]],
                                     meanSub1,
                                     meanSub2,
                                     meanSub3,
                                     meanSub4),
                      lat        = rep(DR4$Lat[DR4$event == earthquakes[e]], 5),
                      quake      = rep(earthquakes[e], 5*Nj),
                      model      = c(rep("Estimate", Nj),
                                     rep("Model 1", Nj),
                                     rep("Model 2", Nj),
                                     rep("Model 3", Nj),
                                     rep("Model 4", Nj)))
  
  plotDF = rbind(plotDF, thisDF)
}

plotDF$quake = factor(plotDF$quake,
                      levels = earthquakes,
                      ordered = T)

g = ggplot(plotDF) +
  geom_point(aes(x=subsidence, y=lat, colour=model)) +
  facet_wrap(~quake, nrow=2) +
  labs(x = "Subsidence (m)",
       y = "Latitude (°)",
       colour= "") +
  scale_colour_viridis_d()
plot(g)

# Now plot all the full model subsidence for comparison
setwd("~/Uni/NTNU/Masters Project/CSZ/R/Code/MultiQuakeCSZ/Backup Results")
load("Backup Subsidence.RData")
subsidence1 = okadaSubs

setwd("~/Uni/NTNU/Masters Project/CSZ/R/Code/MultiQuakeSharedCSZ/Backup Results")
load("Backup Subsidence.RData")
subsidence2 = okadaSubs

setwd("~/Uni/NTNU/Masters Project/CSZ/R/Code/MultiQuakeAnisoCSZ/Backup Results")
load("Backup Subsidence.RData")
subsidence3 = okadaSubs

setwd("~/Uni/NTNU/Masters Project/CSZ/R/Code/MultiQuakeSharedAnisoCSZ/Backup Results")
load("Backup Subsidence.RData")
subsidence4 = okadaSubs

setwd("~/Uni/NTNU/Masters Project/CSZ/R/Code/MultiQuakeCSZ")
load("DR4.RData")

E = 12

plotDF = data.frame()
earthquakes = c("T1","T2","T3","T4","T5","T6","T7","T8","T9","T10","T11","T12")
for (e in 1:E){
  Nj = dim(foldSubsidence1[[e]])[1]
  
  # mean subsidence
  meanSub1 = rowMeans(subsidence1[[e]])
  meanSub2 = rowMeans(subsidence2[[e]])
  meanSub3 = rowMeans(subsidence3[[e]])
  meanSub4 = rowMeans(subsidence4[[e]])
  
  thisDF = data.frame(subsidence = c(DR4$subsidence[DR4$event == earthquakes[e]],
                                     meanSub1,
                                     meanSub2,
                                     meanSub3,
                                     meanSub4),
                      lat        = rep(DR4$Lat[DR4$event == earthquakes[e]], 5),
                      quake      = rep(earthquakes[e], 5*Nj),
                      model      = c(rep("Estimate", Nj),
                                     rep("Model 1", Nj),
                                     rep("Model 2", Nj),
                                     rep("Model 3", Nj),
                                     rep("Model 4", Nj)))
  
  plotDF = rbind(plotDF, thisDF)
}

plotDF$quake = factor(plotDF$quake,
                      levels = earthquakes,
                      ordered = T)

g = ggplot(plotDF) +
  geom_point(aes(x=subsidence, y=lat, colour=model)) +
  facet_wrap(~quake, nrow=2) +
  labs(x = "Subsidence (m)",
       y = "Latitude (°)",
       colour= "") +
  scale_colour_viridis_d()
plot(g)

#-------------------------------------------------------------------------------
#---------------------------Score Summary Plots---------------------------------

setwd("~/Uni/NTNU/Masters Project/CSZ/R/Code/MultiQuakeCSZ/Backup Results")
load("Backup Validation Scores.RData")

setwd("~/Uni/NTNU/Masters Project/CSZ/R/Code/MultiQuakeSharedCSZ/Backup Results")
load("Backup Validation Scores.RData")

setwd("~/Uni/NTNU/Masters Project/CSZ/R/Code/MultiQuakeAnisoCSZ/Backup Results")
load("Backup Validation Scores.RData")

setwd("~/Uni/NTNU/Masters Project/CSZ/R/Code/MultiQuakeSharedAnisoCSZ/Backup Results")
load("Backup Validation Scores.RData")

allScores = rbind(validScoreDF1,
                  validScoreDF2,
                  validScoreDF3,
                  validScoreDF4)
allScores$model = c(rep("Model 1", dim(validScoreDF1)[1]),
                    rep("Model 2", dim(validScoreDF2)[1]),
                    rep("Model 3", dim(validScoreDF3)[1]),
                    rep("Model 4", dim(validScoreDF4)[1]))

# Plot CRPS
ggplot(allScores, aes(x=Event,y=CRPS,fill=model)) +
  geom_boxplot(outlier.shape = NA) +
  stat_summary(fun=mean,
               geom='point',
               shape=20,
               size=3,
               aes(colour=model),
               position = position_dodge(width = 0.75)) +
  scale_fill_viridis_d(option="plasma") +
  scale_colour_viridis_d(option="viridis", direction = -1) +
  labs(fill="Model", colour="Model Mean")

# Plot Squared Error
ggplot(allScores, aes(x=Event,y=SE,fill=model)) +
  geom_boxplot(outlier.shape = NA) +
  stat_summary(fun=mean,
               geom='point',
               shape=20,
               size=3,
               aes(colour=model),
               position = position_dodge(width = 0.75)) +
  scale_fill_viridis_d(option="plasma") +
  scale_colour_viridis_d(option="viridis", direction = -1) +
  labs(y="Squared Error",fill="Model", colour="Model Mean")

# Plot Absolute Error
ggplot(allScores, aes(x=Event,y=AE,fill=model)) +
  geom_boxplot(outlier.shape = NA) +
  stat_summary(fun=mean,
               geom='point',
               shape=20,
               size=3,
               aes(colour=model),
               position = position_dodge(width = 0.75)) +
  scale_fill_viridis_d(option="plasma") +
  scale_colour_viridis_d(option="viridis", direction = -1) +
  labs(y="Absolute Error",fill="Model", colour="Model Mean")

# Plot Interval Score
ggplot(allScores, aes(x=Event,y=Interval,fill=model)) +
  geom_boxplot(outlier.shape = NA) +
  stat_summary(fun=mean,
               geom='point',
               shape=20,
               size=3,
               aes(colour=model),
               position = position_dodge(width = 0.75)) +
  scale_fill_viridis_d(option="plasma") +
  scale_colour_viridis_d(option="viridis", direction = -1) +
  labs(y="Interval Score",fill="Model", colour="Model Mean")

# Find event averages Model 1
validScoreDF1 = base::subset(validScoreDF1, select=-Logs)
eventScore1 = aggregate(cbind(SE, AE, CRPS, DSS, Coverage, Width, Interval) ~ Event,
                        data = validScoreDF1,
                        FUN  = mean)
print(eventScore1)

# Find event averages Model 2
validScoreDF2 = base::subset(validScoreDF2, select=-Logs)
eventScore2 = aggregate(cbind(SE, AE, CRPS, DSS, Coverage, Width, Interval) ~ Event,
                        data = validScoreDF2,
                        FUN  = mean)
print(eventScore2)

# Find event averages Model 3
validScoreDF3 = base::subset(validScoreDF3, select=-Logs)
eventScore3 = aggregate(cbind(SE, AE, CRPS, DSS, Coverage, Width, Interval) ~ Event,
                        data = validScoreDF3,
                        FUN  = mean)
print(eventScore3)

# Find event averages Model 4
validScoreDF4 = base::subset(validScoreDF4, select=-Logs)
eventScore4 = aggregate(cbind(SE, AE, CRPS, DSS, Coverage, Width, Interval) ~ Event,
                        data = validScoreDF4,
                        FUN  = mean)
print(eventScore4)

# Now calculate overall score by model
allEventScores = rbind(eventScore1,
                       eventScore2,
                       eventScore3,
                       eventScore4)
allEventScores$model = c(rep("Model 1", dim(eventScore1)[1]),
                         rep("Model 2", dim(eventScore2)[1]),
                         rep("Model 3", dim(eventScore3)[1]),
                         rep("Model 4", dim(eventScore4)[1]))

finalScores = aggregate(cbind(SE, AE, CRPS, DSS, Coverage, Width, Interval) ~ model,
                        data = allEventScores,
                        FUN  = mean)
print(finalScores)
