

# Precision Q = t(K) %*% solve(K) %*% K
# K = kappa^2 * C + G
# THEREFORE: Q = kappa^4 * C + 2*kappa^2 * G + G %*% solve(C) %*% G
# Re-parameterized as:   
# M0 = C
# M1 = G
# M2 = G %*% solve(C) %*% G

library(INLA)
library(abind)

n_s = 25
set.seed( 1 )
loc_mesh = matrix( runif(n_s*2), ncol=2)

  mesh_stations = inla.mesh.create( loc_mesh, plot.delay=NULL, extend=list(n=8,offset=-0.15), refine=F )  # loc_samp
  #mesh_stations = inla.mesh.create( loc_mesh, plot.delay=NULL, extend=list(n=8,offset=-0.15), refine=list(min.angle=26,max.edge.data=0.08,max.edge.extra=0.2) )  # loc_samp
  n_i = nrow(mesh_stations$loc)

  # R-INLA constructor
  spde_stations = inla.spde2.matern(mesh_stations, alpha=2)  # Given 2D field, alpha=2 <=> Matern Nu=1 (nu = alpha - D/2: Lindgren and Rue 2013, between Eq. 1&2)

  # My re-coding in R
    TV = mesh_stations$graph$tv       # Triangle to vertex indexing
    V0 = mesh_stations$loc[TV[,1],]   # V = vertices for each triangle
    V1 = mesh_stations$loc[TV[,2],]
    V2 = mesh_stations$loc[TV[,3],]
    E0 = V2 - V1                      # E = edge for each triangle
    E1 = V0 - V2
    E2 = V1 - V0
    EdgeTF = cbind( !(TV[,2]%in%mesh_stations$idx$loc|TV[,3]%in%mesh_stations$idx$loc), !(TV[,1]%in%mesh_stations$idx$loc|TV[,3]%in%mesh_stations$idx$loc), !(TV[,1]%in%mesh_stations$idx$loc|TV[,2]%in%mesh_stations$idx$loc) ) # T: triangle I has edge pos2=E0; pos2=E1, pos3=E2; F=not edge
    H = diag( rep(1,3) )
    # Calculate G0
    TmpFn = function(Vec1,Vec2) abs(det( rbind(rep(1,3),Vec1,Vec2) ))
    T = rep(NA, nrow(E0))
    for(i in 1:length(T)) T[i] = TmpFn( E0[i,],E1[i,] )/2   # T = area of each triangle
    M0 = diag( tapply( rep(T/3,3), INDEX=TV, FUN=sum) )  # Expectation: spde_stations$param.inla$M0
    # Calculate G1
    Gtmp = array(NA, dim=c(nrow(TV),3,3))
    for(i in 1:dim(Gtmp)[1]) Gtmp[i,,] = rbind( E0[i,], E1[i,], E2[i,]) %*% H %*% cbind( E0[i,], E1[i,], E2[i,]) / (4*T[i])
    # Calculate boundary effects (NOT SURE HOW TO USE THIS!)
    Btmp = array(NA, dim=c(nrow(TV),3,3))
    for(i in 1:dim(Btmp)[1]) Btmp[i,,] = t(rbind(cbind(0,E0[i,],E0[i,]),cbind(E1[i,],0,E1[i,]),cbind(E2[i,],E2[i,],0))) %*% rbind(diag(rep(EdgeTF[i,1],3)),diag(rep(EdgeTF[i,2],3)),diag(rep(EdgeTF[i,3],3))) %*% cbind(E0[i,],E1[i,],E2[i,]) / (-4*T[i])
    # Sum across triangles
    Gtmp2 = array(0, dim=c(n_i,n_i))
    for(k in 1:nrow(TV)){
      Gtmp2[TV[k,2],TV[k,1]] = Gtmp2[TV[k,1],TV[k,2]] = Gtmp2[TV[k,1],TV[k,2]] + (Gtmp[k,1,2])  
      Gtmp2[TV[k,3],TV[k,2]] = Gtmp2[TV[k,2],TV[k,3]] = Gtmp2[TV[k,2],TV[k,3]] + (Gtmp[k,2,3])
      Gtmp2[TV[k,3],TV[k,1]] = Gtmp2[TV[k,1],TV[k,3]] = Gtmp2[TV[k,1],TV[k,3]] + (Gtmp[k,1,3])
      Gtmp2[TV[k,1],TV[k,1]] = Gtmp2[TV[k,1],TV[k,1]] + (Gtmp[k,1,1])
      Gtmp2[TV[k,2],TV[k,2]] = Gtmp2[TV[k,2],TV[k,2]] + (Gtmp[k,2,2])
      Gtmp2[TV[k,3],TV[k,3]] = Gtmp2[TV[k,3],TV[k,3]] + (Gtmp[k,3,3])
    }
    M1 = Gtmp2                         # Expectation: M1==spde_stations$param.inla$M1
    # Check correspondance: plot( as.vector(Gtmp2), as.vector(spde_stations$param.inla$M1)); cor( as.vector(Gtmp2), as.vector(spde_stations$param.inla$M1))
    M2 = M1 %*% diag(1/diag(M0)) %*% M1 # Expectation: spde_stations$param.inla$M2
    #spde_stations = list( "param.inla"=list("M0"=as(M0,"dgTMatrix"), "M1"=as(M1,"dgTMatrix"), "M2"=as(M2,"dgTMatrix")), "n.spde"=n_i )

    # Plot vertex numbers
    plot(mesh_stations); text( x=mesh_stations$loc[,1], y=mesh_stations$loc[,2], labels=1:n_i)
    Loc = apply(abind(V0[,1:2],V1[,1:2],V2[,1:2],along=3),MARGIN=1:2,FUN=mean)
    text( x=Loc[,1], y=Loc[,2], labels=1:nrow(Loc), col="red")

summary( as.vector(M0 - spde_stations$param.inla$M0) )
summary( as.vector(M1 - spde_stations$param.inla$M1) )
summary( as.vector(M2 - spde_stations$param.inla$M2) )
  