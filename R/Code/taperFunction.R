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

taperSimple = function(lambda, depth){
  return(exp(-lambda*depth))
}

depths = c(0:30000) # depths in meters to calculate taper over

# create tapers
l1 = 20
y1 = taperPaper(l1, depths)

l2 = 1/2000
y2 = taperNew(l2, depths)

l3 = 1/10000
y3 = taperSimple(l3, depths)

data = data.frame(depth=depths,
                  y1=y1, y2=y2, y3=y3)

# plot using plotly
fig = plot_ly(data, x = ~depth) 
#fig = fig %>% add_trace(y = ~y1, name = 'Old', type="scatter", mode = 'lines')
fig = fig %>% add_trace(y = ~y2, name = '1 - exp(lambda*(depth - 30000))', type="scatter", mode = 'lines')
fig = fig %>% add_trace(y = ~y3, name = 'exp(-lambda*depth)', type="scatter", mode = 'lines')
fig = fig %>% layout(title = "Example Taper Functions",
                      xaxis = list(title = "Depth (m)"),
                      yaxis = list (title = "Taper"))
fig

