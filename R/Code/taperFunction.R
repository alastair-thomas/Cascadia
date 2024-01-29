require(plotly)
# The taper function used in https://arxiv.org/pdf/1605.02863.pdf
taperPaper = function(lambda, depth){
  dmax = 22500
  return(1 - exp(-lambda*(depth - dmax)/dmax))
}

# A new taper with fewer parameters
# Not sure it can give the correct shape
taperNew = function(lambda, depth){
  return(1 - exp(lambda*(depth - 30000)))
}

# The one we actually use
taperSimple = function(lambda, depth){
  return(exp(-lambda*depth))
}

depths = c(0:30000)/1000 # depths in km to calculate taper over

# create tapers
l1 = 20
y1 = taperPaper(l1, depths)

l2 = 1/2000
y2 = taperNew(l2, depths)




l3 = exp(-2.2)
y3 = taperSimple(l3, depths)

data = data.frame(depth=depths, y3=y3)

# plot using plotly

fig = ggplot(data, aes(x=depth, y=y3)) +
  geom_line(color="blue") +
  labs(x = "Depth (km)", y = "Taper")
plot(fig)

