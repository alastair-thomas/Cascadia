
getScoresFromSubsidence = function(foldSubsidence, data, nSims=1000){
  
  E = 12
  
  # Calculate the scores
  # mse, mae, crps, empirical coverage, width
  foldErrors = list()
  for (fold in 1:E){
    foldErrors[[fold]] = list()
    
    yObs  = data$subsidence[[fold]]
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

# load validation subsidence
setwd("~/Uni/NTNU/Masters Project/CSZ/R/Code/MultiQuakeCSZ/Results")
load("Validation Subsidences.RData")

# set constants
E = length(data$subsidence)
S = length(fault)
B = dim(spdeMesh$loc)[1]
nSims = 1000

# get the validation model scores
validScoreDF1 = getScoresFromSubsidence(foldSubsidence1, data)
setwd("~/Uni/NTNU/Masters Project/CSZ/R/Code/MultiQuakeCSZ/Results")
save(validScoreDF1, file="Validation Scores.RData")

# get the scores for the full model
load("Subsidence.RData")
subsidence1 = okadaSubs
fullScoreDF1 = getScoresFromSubsidence(subsidence1, data)
setwd("~/Uni/NTNU/Masters Project/CSZ/R/Code/MultiQuakeCSZ/Results")
save(fullScoreDF1, file="Full Scores.RData")

plotInnerOuterError(fullScoreDF1, validScoreDF1)

#-------------------------------------------------------------------------------
#---------------------------MultiQuakeSharedCSZ---------------------------------

# load data required
setwd("~/Uni/NTNU/Masters Project/CSZ/R/Code/MultiQuakeSharedCSZ")
load("Fault.RData")
load("spdeMesh.RData")
load("Data.RData")

# load validation subsidence
setwd("~/Uni/NTNU/Masters Project/CSZ/R/Code/MultiQuakeSharedCSZ/Results")
load("Validation Subsidences.RData")

# set constants
E = length(data$subsidence)
S = length(fault)
B = dim(spdeMesh$loc)[1]
nSims = 1000

# get the validation model scores
validScoreDF2 = getScoresFromSubsidence(foldSubsidence2, data)
setwd("~/Uni/NTNU/Masters Project/CSZ/R/Code/MultiQuakeSharedCSZ/Results")
save(validScoreDF2, file="Validation Scores.RData")

# get the scores for the full model
setwd("~/Uni/NTNU/Masters Project/CSZ/R/Code/MultiQuakeSharedCSZ/Results")
load("Subsidence.RData")
subsidence2 = subsidences
fullScoreDF2 = getScoresFromSubsidence(subsidence2, data)
setwd("~/Uni/NTNU/Masters Project/CSZ/R/Code/MultiQuakeSharedCSZ/Results")
save(fullScoreDF2, file="Full Scores.RData")

plotInnerOuterError(fullScoreDF2, validScoreDF2)

#-------------------------------------------------------------------------------
#---------------------------MultiQuakeAnisoCSZ----------------------------------

# load data required
setwd("~/Uni/NTNU/Masters Project/CSZ/R/Code/MultiQuakeAnisoCSZ")
load("Fault.RData")
load("spdeMesh.RData")
load("Data.RData")

# load validation draws
setwd("~/Uni/NTNU/Masters Project/CSZ/R/Code/MultiQuakeAnisoCSZ/Results")
load("Validation Subsidences.RData")

# set constants
E = length(data$subsidence)
S = length(fault)
B = dim(spdeMesh$loc)[1]
nSims = 1000

# get the validation model scores
validScoreDF3 = getScoresFromSubsidence(foldSubsidence3, data)
setwd("~/Uni/NTNU/Masters Project/CSZ/R/Code/MultiQuakeAnisoCSZ/Results")
save(validScoreDF3, file="Validation Scores.RData")

# get the scores for the full model
setwd("~/Uni/NTNU/Masters Project/CSZ/R/Code/MultiQuakeAnisoCSZ/Results")
load("Subsidence.RData")
subsidence3 = okadaSubs
fullScoreDF3 = getScoresFromSubsidence(subsidence3, data)
setwd("~/Uni/NTNU/Masters Project/CSZ/R/Code/MultiQuakeAnisoCSZ/Results")
save(fullScoreDF3, file="Full Scores.RData")

plotInnerOuterError(fullScoreDF3, validScoreDF3)

#-------------------------------------------------------------------------------
#--------------------------MultiQuakeSharedAnisoCSZ-----------------------------

# load data required
setwd("~/Uni/NTNU/Masters Project/CSZ/R/Code/MultiQuakeSharedAnisoCSZ")
load("Fault.RData")
load("spdeMesh.RData")
load("Data.RData")

# load validation draws
setwd("~/Uni/NTNU/Masters Project/CSZ/R/Code/MultiQuakeSharedAnisoCSZ/Results")
load("Validation Subsidences.RData")

# set constants
E = length(data$subsidence)
S = length(fault)
B = dim(spdeMesh$loc)[1]
nSims = 1000

for (j in 1:12){
  print(range(rowMeans(foldSubsidence4[[j]])))
}
# get the validation model scores
validScoreDF4 = getScoresFromSubsidence(foldSubsidence4, data)
setwd("~/Uni/NTNU/Masters Project/CSZ/R/Code/MultiQuakeSharedAnisoCSZ/Results")
save(validScoreDF4, file="Validation Scores.RData")

# get the scores for the full model
setwd("~/Uni/NTNU/Masters Project/CSZ/R/Code/MultiQuakeSharedAnisoCSZ/Results")
load("Subsidence.RData")
subsidence4 = okadaSubs
fullScoreDF4 = getScoresFromSubsidence(subsidence4, data)
setwd("~/Uni/NTNU/Masters Project/CSZ/R/Code/MultiQuakeSharedAnisoCSZ/Results")
save(fullScoreDF4, file="Full Scores.RData")

plotInnerOuterError(fullScoreDF4, validScoreDF4)

#-------------------------------------------------------------------------------

#---------------------------Subsidence Summary Plots----------------------------

# Plot all the validation subsidences
setwd("~/Uni/NTNU/Masters Project/CSZ/R/Code/MultiQuakeCSZ/Results")
load("Validation Subsidences.RData")

setwd("~/Uni/NTNU/Masters Project/CSZ/R/Code/MultiQuakeSharedCSZ/Results")
load("Validation Subsidences.RData")

setwd("~/Uni/NTNU/Masters Project/CSZ/R/Code/MultiQuakeAnisoCSZ/Results")
load("Validation Subsidences.RData")

setwd("~/Uni/NTNU/Masters Project/CSZ/R/Code/MultiQuakeSharedAnisoCSZ/Results")
load("Validation Subsidences.RData")

setwd("~/Uni/NTNU/Masters Project/CSZ/R/Code/MultiQuakeCSZ")
load("DR4.RData")

E = 12

earthquakes = c("T1","T2","T3","T4","T5","T6",
                "T7","T8","T9","T10","T11","T12")
subDF1 = data.frame()
subDF2 = data.frame()
subDF3 = data.frame()
subDF4 = data.frame()

for (e in 1:E){
  Nj = dim(foldSubsidence1[[e]])[1]
  
  # mean subsidence
  meanSub1 = rowMeans(foldSubsidence1[[e]])
  meanSub2 = rowMeans(foldSubsidence2[[e]])
  meanSub3 = rowMeans(foldSubsidence3[[e]])
  meanSub4 = rowMeans(foldSubsidence4[[e]])
  
  # medSub1 = apply(foldSubsidence1[[e]], 1, median)
  # medSub2 = apply(foldSubsidence2[[e]], 1, median)
  # medSub3 = apply(foldSubsidence3[[e]], 1, median)
  # medSub4 = apply(foldSubsidence4[[e]], 1, median)
  
  lowSub1 = apply(foldSubsidence1[[e]], 1, quantile, 0.025)
  lowSub2 = apply(foldSubsidence2[[e]], 1, quantile, 0.025)
  lowSub3 = apply(foldSubsidence3[[e]], 1, quantile, 0.025)
  lowSub4 = apply(foldSubsidence4[[e]], 1, quantile, 0.025)
  
  uppSub1 = apply(foldSubsidence1[[e]], 1, quantile, 0.975)
  uppSub2 = apply(foldSubsidence2[[e]], 1, quantile, 0.975)
  uppSub3 = apply(foldSubsidence3[[e]], 1, quantile, 0.975)
  uppSub4 = apply(foldSubsidence4[[e]], 1, quantile, 0.975)
  
  # res1 = DR4$subsidence[DR4$event == earthquakes[e]] - meanSub1
  # res2 = DR4$subsidence[DR4$event == earthquakes[e]] - meanSub2
  # res3 = DR4$subsidence[DR4$event == earthquakes[e]] - meanSub3
  # res4 = DR4$subsidence[DR4$event == earthquakes[e]] - meanSub4
  
  thisDF1 = data.frame(subsidence = c(DR4$subsidence[DR4$event == earthquakes[e]],
                                      meanSub1,
                                      lowSub1,
                                      uppSub1),
                      lat         = rep(DR4$Lat[DR4$event == earthquakes[e]], 4),
                      quake       = rep(earthquakes[e], 4*Nj),
                      dataType    = c(rep("Estimate", Nj),
                                      rep("Mean", Nj),
                                      rep("Lower", Nj),
                                      rep("Upper", Nj)))
  
  thisDF2 = data.frame(subsidence = c(DR4$subsidence[DR4$event == earthquakes[e]],
                                      meanSub2,
                                      lowSub2,
                                      uppSub2),
                       lat         = rep(DR4$Lat[DR4$event == earthquakes[e]], 4),
                       quake       = rep(earthquakes[e], 4*Nj),
                       dataType    = c(rep("Estimate", Nj),
                                       rep("Mean", Nj),
                                       rep("Lower", Nj),
                                       rep("Upper", Nj)))
  
  thisDF3 = data.frame(subsidence = c(DR4$subsidence[DR4$event == earthquakes[e]],
                                      meanSub3,
                                      lowSub3,
                                      uppSub3),
                       lat         = rep(DR4$Lat[DR4$event == earthquakes[e]], 4),
                       quake       = rep(earthquakes[e], 4*Nj),
                       dataType    = c(rep("Estimate", Nj),
                                       rep("Mean", Nj),
                                       rep("Lower", Nj),
                                       rep("Upper", Nj)))
  
  thisDF4 = data.frame(subsidence = c(DR4$subsidence[DR4$event == earthquakes[e]],
                                      meanSub4,
                                      lowSub4,
                                      uppSub4),
                       lat         = rep(DR4$Lat[DR4$event == earthquakes[e]], 4),
                       quake       = rep(earthquakes[e], 4*Nj),
                       dataType    = c(rep("Estimate", Nj),
                                       rep("Mean", Nj),
                                       rep("Lower", Nj),
                                       rep("Upper", Nj)))
  
  subDF1 = rbind(subDF1, thisDF1)
  subDF2 = rbind(subDF2, thisDF2)
  subDF3 = rbind(subDF3, thisDF3)
  subDF4 = rbind(subDF4, thisDF4)
}

subDF1$quake = factor(subDF1$quake,
                      levels = earthquakes,
                      ordered = T)
subDF2$quake = factor(subDF2$quake,
                      levels = earthquakes,
                      ordered = T)
subDF3$quake = factor(subDF3$quake,
                      levels = earthquakes,
                      ordered = T)
subDF4$quake = factor(subDF4$quake,
                      levels = earthquakes,
                      ordered = T)

plotDF1 = subDF1 %>%
  filter(dataType %in% c("Mean", "Estimate"))
g = ggplot(plotDF1) +
  geom_point(aes(x=subsidence, y=lat, colour=dataType, shape=dataType)) +
  scale_colour_manual(values = c("Mean"="#332288", "Estimate"="#44AA99")) +
  scale_shape_manual(values = c("Mean"=17, "Estimate"=1)) +
  facet_wrap(~quake, nrow=2) +
  labs(x = "Subsidence (m)",
       y = "Latitude (°)") +
  theme(legend.position = "none")
plot(g)

plotDF2 = subDF2 %>%
  filter(dataType %in% c("Mean", "Estimate"))
g = ggplot(plotDF2) +
  geom_point(aes(x=subsidence, y=lat, colour=dataType, shape=dataType)) +
  scale_colour_manual(values = c("Mean"="#AA4499", "Estimate"="#44AA99")) +
  scale_shape_manual(values = c("Mean"=17, "Estimate"=1)) +
  facet_wrap(~quake, nrow=2) +
  labs(x = "Subsidence (m)",
       y = "Latitude (°)") +
  theme(legend.position = "none")
plot(g)

plotDF3 = subDF3 %>%
  filter(dataType %in% c("Mean", "Estimate"))
g = ggplot(plotDF3) +
  geom_point(aes(x=subsidence, y=lat, colour=dataType, shape=dataType)) +
  scale_colour_manual(values = c("Mean"="#882255", "Estimate"="#44AA99")) +
  scale_shape_manual(values = c("Mean"=17, "Estimate"=1)) +
  facet_wrap(~quake, nrow=2) +
  labs(x = "Subsidence (m)",
       y = "Latitude (°)") +
  theme(legend.position = "none")
plot(g)

plotDF4 = subDF4 %>%
  filter(dataType %in% c("Mean", "Estimate"))
g = ggplot(plotDF4) +
  geom_point(aes(x=subsidence, y=lat, colour=dataType, shape=dataType)) +
  scale_colour_manual(values = c("Mean"="#661100", "Estimate"="#44AA99")) +
  scale_shape_manual(values = c("Mean"=17, "Estimate"=1)) +
  facet_wrap(~quake, nrow=2) +
  labs(x = "Subsidence (m)",
       y = "Latitude (°)") +
  theme(legend.position = "none")
plot(g)


#-------------------------------------------------------------------------------
#---------------------------Score Summary Plots---------------------------------

setwd("~/Uni/NTNU/Masters Project/CSZ/R/Code/MultiQuakeCSZ/Results")
load("Validation Scores.RData")

setwd("~/Uni/NTNU/Masters Project/CSZ/R/Code/MultiQuakeSharedCSZ/Results")
load("Validation Scores.RData")

setwd("~/Uni/NTNU/Masters Project/CSZ/R/Code/MultiQuakeAnisoCSZ/Results")
load("Validation Scores.RData")

setwd("~/Uni/NTNU/Masters Project/CSZ/R/Code/MultiQuakeSharedAnisoCSZ/Results")
load("Validation Scores.RData")

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
  labs(y="Squared Error",fill="Model", colour="Model Mean") +
  ylim(c(0, 5))

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
  labs(y="Interval Score",fill="Model", colour="Model Mean") +
  ylim(c(0, 0.2))

# Find event averages Model 1
eventScore1 = aggregate(cbind(SE, AE, CRPS, Coverage, Width, Interval) ~ Event,
                        data = validScoreDF1,
                        FUN  = mean)
print(eventScore1)

# Find event averages Model 2
eventScore2 = aggregate(cbind(SE, AE, CRPS, Coverage, Width, Interval) ~ Event,
                        data = validScoreDF2,
                        FUN  = mean)
print(eventScore2)

# Find event averages Model 3
eventScore3 = aggregate(cbind(SE, AE, CRPS, Coverage, Width, Interval) ~ Event,
                        data = validScoreDF3,
                        FUN  = mean)
print(eventScore3)

# Find event averages Model 4
eventScore4 = aggregate(cbind(SE, AE, CRPS, Coverage, Width, Interval) ~ Event,
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

finalScores = aggregate(cbind(SE, AE, CRPS, Coverage, Width, Interval) ~ model,
                        data = allEventScores,
                        FUN  = mean)
print(finalScores)
