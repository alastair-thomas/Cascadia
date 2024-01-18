
setwd("C://users/alast/OneDrive/Documents/Uni/NTNU/Masters Project/CSZ/R")
load("Data/Subsidence/DR1.RData")
dr1$Site = sub("Port Al berni", "Port Alberni", dr1$Site)
  
Mode = function(x) {
  ux = unique(x)
  return(ux[which.max(tabulate(match(x, ux)))])
}


plotSubsidence = function(scale=2){
  # gets the base map of the CSZ
  g = plotBase(scale=scale)
  
  Sites = data.frame(Site=dr1$Site, Lat=dr1$Lat, Lon=dr1$Lon, Sub=dr1$subsidence, Unc=dr1$Uncertainty)
  Sites = aggregate(.~Site, data=Sites, mean)
  Sites2 = data.frame(Site=dr1$Site, Quality=dr1$quality)
  Sites2 = aggregate(.~Site, data=Sites2, Mode)
  Sites$Quality = factor(Sites2$Quality, levels=c(1,2,3), labels=c("High", "Medium", "Low"))
  
  g = g +
    geom_point(data=Sites, aes(x=Lon, y=Lat, colour=Unc, shape=Quality), size=scale*1.5) +
    scale_colour_gradient(low="red", high="green", name="Uncertainty (m)", trans = "reverse") +
    geom_text_repel(data=Sites, aes(x=Lon, y=Lat, label=Site),
                    size=scale*2, hjust=0,
                    xlim = c(-130, -127), force=5,
                    direction="y", box.padding=0.15) +
    coord_sf(xlim=-c(131, 120), ylim=c(40, 50)) +
    labs(shape='Quality') +
    theme(legend.position = "right", legend.key.height = unit(2, 'cm'))
  
  return(g)
}

# a function to plot the number of subsidence estimates for each earthquake
plotEvents = function(data){
  df = data.frame(data) # makes sure using a dataframe
  df = df[!grepl("a|b|R", df$event), ] # gets only full margin ruptures
  
  eventsOrder = c("T1", "T2", "T3", "T4", "T5", "T6", "T7", "T8", "T9", "T10", "T11", "T12")
  df$event = factor(df$event, levels=eventsOrder)
  
  ggplot(df, aes(x = event)) +
    geom_bar() +
    labs(title = "",
         x = "Seismic Event",
         y = "Number of Subsidence Estimates")
}

# plot the number of observations per site on the map
# size of spot is number of observations
# only for the full margin ruptures.

plotNumObersvations = function(scale=2){
  # gets the base map of the CSZ
  g = plotBase(scale=scale)
  
  Sites = data.frame(Site=dr1$Site, Lat=dr1$Lat, Lon=dr1$Lon)
  Sites = aggregate(.~Site, data=Sites, mean)
  Sites$count = as.numeric(data.frame(table(dr1$Site))$Freq)
  
  g = g +
    geom_point(data=Sites, aes(x=Lon, y=Lat, size=count, alpha=count, colour=count)) +
    scale_size_continuous(range = c(1, scale*8)) +
    scale_alpha_continuous(range = c(0.5, 1), trans="reverse") +
    scale_colour_gradient2(low="red", mid="yellow", high="green", midpoint=50) +
    coord_sf(xlim=-c(131, 120), ylim=c(40, 50)) + 
    guides(
      size = guide_legend(title = "Number of\nSubsidence Estimates"),
      alpha = guide_legend(title = "Number of\nSubsidence Estimates"),
      colour = guide_legend(title = "Number of\nSubsidence Estimates")
    )
  
  return(g)
}

