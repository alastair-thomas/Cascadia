
h0 = seq(0.01,10,0.1)
h1 = seq(0.01,10,0.1)
h = expand.grid(x=h0,y=h1)

angles = c()
stretchs = c()
for (i in 1:dim(h)[1]){
  h0 = h[i,1]
  h1 = h[i,2]
  
  H = matrix(c(h0,      log(h1),
               log(h1), (1+log(h1)^2)/h0),
             nrow=2)
  
  eigen_result = eigen(H)
  eigenvalues = eigen_result$values
  eigenvectors = eigen_result$vectors
  
  angles = c(angles, rad2deg(atan(eigenvectors[1,1]/eigenvectors[2,1])))
  stretchs = c(stretchs, sqrt(eigenvalues[1]))
}

plotDF = data.frame(h0=h[,1],
                    h1=h[,2],
                    h01 = h[,1]/h[,2],
                    stretch=stretchs,
                    angle=angles)

ggplot(plotDF)+
  geom_point(aes(x=h1,y=angle,colour=h0))+
  scale_colour_viridis()

ggplot(plotDF)+
  geom_point(aes(x=h0,y=stretch,colour=h1))+
  scale_colour_viridis()

ggplot(plotDF[plotDF$h01 < 10,])+
  geom_point(aes(x=h01, y=stretch))+
  scale_colour_viridis()

plotDF[(plotDF$angle < -10) & (plotDF$angle > -14),]


h0 = 5
h1 = 3
H = matrix(c(h0,      log(h1),
             log(h1), (1+log(h1)^2)/h0),
           nrow=2)
# Compute eigenvalues and eigenvectors
eigen_result = eigen(H)
eigenvalues = eigen_result$values
eigenvectors = eigen_result$vectors
print(rad2deg(atan(eigenvectors[1,1]/eigenvectors[2,1])))
print(sqrt(eigenvalues))

# plot ellipse
major_length <- sqrt(eigenvalues[1])
minor_length <- sqrt(eigenvalues[2])
major_direction <- eigenvectors[,1]
minor_direction <- eigenvectors[,2]
angles <- seq(0, 2*pi, length.out = 100)
x <- major_length * cos(angles) * major_direction[1] + minor_length * sin(angles) * minor_direction[1]
y <- major_length * cos(angles) * major_direction[2] + minor_length * sin(angles) * minor_direction[2]
plot(x, y, type = "l", asp = 1, xlab = "X", ylab = "Y", main = "Coviance Structure")


# finding min and max strike angles
varName = "strike"
grid = readGrid(varName) # read data
grid$Z[grid$Z > 100] = grid$Z[grid$Z > 100] - 360 # center around zero
tempGrid = readGrid("depth") # get depths
depths = -tempGrid$Z # make them positive
depthMask = which(depths <= 30) # find points shallower than 30km
grid2 = grid[depthMask,] # get those strike points
print(min(grid2$Z))
print(max(grid2$Z))
