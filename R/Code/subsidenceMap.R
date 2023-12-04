
load("Data/Subsidence/DR1.RData")

Mode = function(x) {
  ux = unique(x)
  return(ux[which.max(tabulate(match(x, ux)))])
}


plotSubsidence = function(scale=2){
  # gets the base map in given projection
  g = plotBase(scale=scale)
  
  Sites <- data.frame(Site=dr1$Site, Lat=dr1$Lat, Lon=dr1$Lon, Sub=dr1$subsidence, Unc=dr1$Uncertainty)
  Sites <- aggregate(.~Site, data=Sites, mean)
  Sites2 <- data.frame(Site=dr1$Site, Quality=dr1$quality)
  Sites2 <- aggregate(.~Site, data=Sites2, Mode)
  Sites$Quality <- factor(Sites2$Quality, levels=c(1,2,3), labels=c("High","Medium", "Low"))
  
  g = g +
    geom_point(data=Sites, aes(x=Lon, y=Lat, colour=Unc, shape=Quality), size=scale*1.5) +
    scale_colour_gradient(low="green", high="red", name="Uncertainty:") +
    #scale_color_manual(values = c("High" = "green", "Medium" = "orange", "Low" = "red"), name="Quality:") +
    geom_text_repel(data=Sites, aes(x=Lon, y=Lat, label=Site),
                    size=scale*2, hjust=0,
                    xlim = c(-130, -127), force=5,
                    direction="y", box.padding=0.15) +
    coord_sf(xlim=-c(131, 120), ylim=c(40, 50)) +
    theme(legend.position="top") +
    labs(shape='Quality:') 
  
  return(g)
}

#ggsave(file="Subsidence.png", path="Plots/", height=6, width=6)

