# This script finds the minimum and maximum of the spatial variables for colour plotting
# The variables are:
#   - Slip Draws: Mean + SD + Width
#   - Whole Subsidence plots
#   - Future Slips: Mean + SD + Width
#   - Future whole subsidence

#-----------------------------------Slips---------------------------------------

# model 1
setwd("~/Uni/NTNU/Masters Project/CSZ/R/Code/MultiQuakeCSZ/Results")
load("Slips.RData")
slipDraws1 = slipDraws
rm(slipDraws)

# model 2
setwd("~/Uni/NTNU/Masters Project/CSZ/R/Code/MultiQuakeSharedCSZ/Results")
load("Slips.RData")
slipDraws2 = slipDraws
rm(slipDraws)

# model 3
setwd("~/Uni/NTNU/Masters Project/CSZ/R/Code/MultiQuakeAnisoCSZ/Results")
load("Slips.RData")
slipDraws3 = slipDraws
rm(slipDraws)

# model 4
setwd("~/Uni/NTNU/Masters Project/CSZ/R/Code/MultiQuakeSharedAnisoCSZ/Results")
load("Slips.RData")
slipDraws4 = slipDraws
rm(slipDraws)

slipMeanLimits  = c(1000, -1000)
slipSDLimits    = c(1000, -1000)
slipWidthLimits = c(1000, -1000)

E = length(slipDraws1)

for (j in 1:E){
  # means
  maxMeanEvent = max(c(max(apply(slipDraws1[[j]], 1, mean)),
                       max(apply(slipDraws2[[j]], 1, mean)),
                       max(apply(slipDraws3[[j]], 1, mean)),
                       max(apply(slipDraws4[[j]], 1, mean))))
  minMeanEvent = min(c(min(apply(slipDraws1[[j]], 1, mean)),
                       min(apply(slipDraws2[[j]], 1, mean)),
                       min(apply(slipDraws3[[j]], 1, mean)),
                       min(apply(slipDraws4[[j]], 1, mean))))
  if (minMeanEvent < slipMeanLimits[1]){
    slipMeanLimits[1] = minMeanEvent
  }
  if (maxMeanEvent > slipMeanLimits[2]){
    slipMeanLimits[2] = maxMeanEvent
  }
  
  # SDs
  maxSDEvent = max(c(max(apply(slipDraws1[[j]], 1, sd)),
                       max(apply(slipDraws2[[j]], 1, sd)),
                       max(apply(slipDraws3[[j]], 1, sd)),
                       max(apply(slipDraws4[[j]], 1, sd))))
  minSDEvent = min(c(min(apply(slipDraws1[[j]], 1, sd)),
                       min(apply(slipDraws2[[j]], 1, sd)),
                       min(apply(slipDraws3[[j]], 1, sd)),
                       min(apply(slipDraws4[[j]], 1, sd))))
  if (minSDEvent < slipSDLimits[1]){
    slipSDLimits[1] = minSDEvent
  }
  if (maxSDEvent > slipSDLimits[2]){
    slipSDLimits[2] = maxSDEvent
  }
  
  # Widths
  maxWidthEvent = max(c(max(apply(slipDraws1[[j]], 1, quantile, 0.975) -
                            apply(slipDraws1[[j]], 1, quantile, 0.025)),
                        max(apply(slipDraws2[[j]], 1, quantile, 0.975) -
                            apply(slipDraws2[[j]], 1, quantile, 0.025)),
                        max(apply(slipDraws3[[j]], 1, quantile, 0.975) -
                            apply(slipDraws3[[j]], 1, quantile, 0.025)),
                        max(apply(slipDraws4[[j]], 1, quantile, 0.975) -
                            apply(slipDraws4[[j]], 1, quantile, 0.025))))
  minWidthEvent = min(c(min(apply(slipDraws1[[j]], 1, quantile, 0.975) -
                            apply(slipDraws1[[j]], 1, quantile, 0.025)),
                        min(apply(slipDraws2[[j]], 1, quantile, 0.975) -
                            apply(slipDraws2[[j]], 1, quantile, 0.025)),
                        min(apply(slipDraws3[[j]], 1, quantile, 0.975) -
                            apply(slipDraws3[[j]], 1, quantile, 0.025)),
                        min(apply(slipDraws4[[j]], 1, quantile, 0.975) -
                            apply(slipDraws4[[j]], 1, quantile, 0.025))))
  
  if (minWidthEvent < slipWidthLimits[1]){
    slipWidthLimits[1] = minWidthEvent
  }
  if (maxWidthEvent > slipWidthLimits[2]){
    slipWidthLimits[2] = maxWidthEvent
  }
  
}

print(slipMeanLimits)  # # 1.17247 34.05697   => [1.0, 34.5]
print(slipSDLimits)    # 0.5715136 18.8986648 => [0.5, 19.0]
print(slipWidthLimits) # 2.173223 71.162589   => [2.0, 71.5]



#---------------------------Whole Fault Subsidences-----------------------------

# model 1
setwd("~/Uni/NTNU/Masters Project/CSZ/R/Code/MultiQuakeCSZ/Results")
load("Whole Subsidence.RData")
wholeSub1 = wholeSubsidence
rm(wholeSubsidence)

# model 2
setwd("~/Uni/NTNU/Masters Project/CSZ/R/Code/MultiQuakeSharedCSZ/Results")
load("Whole Subsidence.RData")
wholeSub2 = wholeSubsidence
rm(wholeSubsidence)

# model 3
setwd("~/Uni/NTNU/Masters Project/CSZ/R/Code/MultiQuakeAnisoCSZ/Results")
load("Whole Subsidence.RData")
wholeSub3 = wholeSubsidence
rm(wholeSubsidence)

# model 4
setwd("~/Uni/NTNU/Masters Project/CSZ/R/Code/MultiQuakeSharedAnisoCSZ/Results")
load("Whole Subsidence.RData")
wholeSub4 = wholeSubsidence
rm(wholeSubsidence)

subMeanLimits  = c(1000, -1000)
subSDLimits    = c(1000, -1000)

E = length(wholeSub1)

for (j in 1:E){
  # means
  maxMeanEvent = max(c(max(apply(wholeSub1[[j]], 1, mean)),
                       max(apply(wholeSub2[[j]], 1, mean)),
                       max(apply(wholeSub3[[j]], 1, mean)),
                       max(apply(wholeSub4[[j]], 1, mean))))
  minMeanEvent = min(c(min(apply(wholeSub1[[j]], 1, mean)),
                       min(apply(wholeSub2[[j]], 1, mean)),
                       min(apply(wholeSub3[[j]], 1, mean)),
                       min(apply(wholeSub4[[j]], 1, mean))))
  if (minMeanEvent < subMeanLimits[1]){
    subMeanLimits[1] = minMeanEvent
  }
  if (maxMeanEvent > subMeanLimits[2]){
    subMeanLimits[2] = maxMeanEvent
  }
  
  # SDs
  maxSDEvent = max(c(max(apply(wholeSub1[[j]], 1, sd)),
                     max(apply(wholeSub2[[j]], 1, sd)),
                     max(apply(wholeSub3[[j]], 1, sd)),
                     max(apply(wholeSub4[[j]], 1, sd))))
  minSDEvent = min(c(min(apply(wholeSub1[[j]], 1, sd)),
                     min(apply(wholeSub2[[j]], 1, sd)),
                     min(apply(wholeSub3[[j]], 1, sd)),
                     min(apply(wholeSub4[[j]], 1, sd))))
  if (minSDEvent < subSDLimits[1]){
    subSDLimits[1] = minSDEvent
  }
  if (maxSDEvent > subSDLimits[2]){
    subSDLimits[2] = maxSDEvent
  }
  
}

print(subMeanLimits)  # -11.056414   2.988366   => [-11.5, 3.0]
print(subSDLimits)    # 0.003942195 5.791293390 => [0, 6.0]




#---------------------------------Future Slips----------------------------------

# model 1
setwd("~/Uni/NTNU/Masters Project/CSZ/R/Code/MultiQuakeCSZ/Results")
load("Future Slips.RData")
slip1 = futureSlips
rm(futureSlips)

# model 2
setwd("~/Uni/NTNU/Masters Project/CSZ/R/Code/MultiQuakeSharedCSZ/Results")
load("Future Slips.RData")
slip2 = futureSlips
rm(futureSlips)

# model 3
setwd("~/Uni/NTNU/Masters Project/CSZ/R/Code/MultiQuakeAnisoCSZ/Results")
load("Future Slips.RData")
slip3 = futureSlips
rm(futureSlips)

# model 4
setwd("~/Uni/NTNU/Masters Project/CSZ/R/Code/MultiQuakeSharedAnisoCSZ/Results")
load("Future Slips.RData")
slip4 = futureSlips
rm(futureSlips)

# means
maxMeanEvent = max(c(max(apply(slip1, 1, mean)),
                     max(apply(slip2, 1, mean)),
                     max(apply(slip3, 1, mean)),
                     max(apply(slip4, 1, mean))))
minMeanEvent = min(c(min(apply(slip1, 1, mean)),
                     min(apply(slip2, 1, mean)),
                     min(apply(slip3, 1, mean)),
                     min(apply(slip4, 1, mean))))
slipMeanLimits = c(minMeanEvent, maxMeanEvent)
  
# SDs
maxSDEvent = max(c(max(apply(slip1, 1, sd)),
                   max(apply(slip2, 1, sd)),
                   max(apply(slip3, 1, sd)),
                   max(apply(slip4, 1, sd))))
minSDEvent = min(c(min(apply(slip1, 1, sd)),
                   min(apply(slip2, 1, sd)),
                   min(apply(slip3, 1, sd)),
                   min(apply(slip4, 1, sd))))
slipSDLimits = c(minSDEvent, maxSDEvent)

  # Widths
maxWidthEvent = max(c(max(apply(slip1, 1, quantile, 0.975) -
                          apply(slip1, 1, quantile, 0.025)),
                      max(apply(slip2, 1, quantile, 0.975) -
                          apply(slip2, 1, quantile, 0.025)),
                      max(apply(slip3, 1, quantile, 0.975) -
                          apply(slip3, 1, quantile, 0.025)),
                      max(apply(slip4, 1, quantile, 0.975) -
                          apply(slip4, 1, quantile, 0.025))))
minWidthEvent = min(c(min(apply(slip1, 1, quantile, 0.975) -
                          apply(slip1, 1, quantile, 0.025)),
                      min(apply(slip2, 1, quantile, 0.975) -
                          apply(slip2, 1, quantile, 0.025)),
                      min(apply(slip3, 1, quantile, 0.975) -
                          apply(slip3, 1, quantile, 0.025)),
                      min(apply(slip4, 1, quantile, 0.975) -
                          apply(slip4, 1, quantile, 0.025))))
slipWidthLimits = c(minWidthEvent, maxWidthEvent)

print(slipMeanLimits)  # 1.157111 28.258540   => [1.0, 28.5]
print(slipSDLimits)    # 0.6336757 18.6190194 => [0.5, 19.0]
print(slipWidthLimits) # 2.400427 69.591727   => [2.0, 70.0]




#-------------------------Future Whole Fault Subsidences------------------------

# model 1
setwd("~/Uni/NTNU/Masters Project/CSZ/R/Code/MultiQuakeCSZ/Results")
load("Future Subsidence.RData")
Sub1 = futureSubs
rm(futureSubs)

# model 2
setwd("~/Uni/NTNU/Masters Project/CSZ/R/Code/MultiQuakeSharedCSZ/Results")
load("Future Subsidence.RData")
Sub2 = futureSubs
rm(futureSubs)

# model 3
setwd("~/Uni/NTNU/Masters Project/CSZ/R/Code/MultiQuakeAnisoCSZ/Results")
load("Future Subsidence.RData")
Sub3 = futureSubs
rm(futureSubs)

# model 4
setwd("~/Uni/NTNU/Masters Project/CSZ/R/Code/MultiQuakeSharedAnisoCSZ/Results")
load("Future Subsidence.RData")
Sub4 = futureSubs
rm(futureSubs)

  # means
futureMeanLims = c(min(c(min(apply(Sub1, 1, mean)),
                        min(apply(Sub2, 1, mean)),
                        min(apply(Sub3, 1, mean)),
                        min(apply(Sub4, 1, mean)))),
                  max(c(max(apply(Sub1, 1, mean)),
                        max(apply(Sub2, 1, mean)),
                        max(apply(Sub3, 1, mean)),
                        max(apply(Sub4, 1, mean)))))
  
futureSDLims = c(min(c(min(apply(Sub1, 1, sd)),
                       min(apply(Sub2, 1, sd)),
                       min(apply(Sub3, 1, sd)),
                       min(apply(Sub4, 1, sd)))),
                 max(c(max(apply(Sub1, 1, sd)),
                       max(apply(Sub2, 1, sd)),
                       max(apply(Sub3, 1, sd)),
                       max(apply(Sub4, 1, sd)))))

print(futureMeanLims)  # -9.180193  2.487337     => [-9.5, 2.5]
print(futureSDLims)    # 0.004692782 5.406538414 => [0.0,  5.5]
