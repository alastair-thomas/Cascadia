
## install missing packages
pkgs <- c('data.table', 'ggplot2', 'RColorBrewer', 'RandomFields',
          'raster', 'TMB', 'viridis')
new.packages <- pkgs[!(pkgs %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
if(!('INLA' %in% installed.packages()[, 'Package'])){
  ## INLA is not on CRAN
  install.packages("INLA",
                   repos=c(getOption("repos"),
                           INLA="https://inla.r-inla-download.org/R/stable"),
                   dep=TRUE)
}
# load packages
invisible(lapply(c(pkgs, 'INLA'), library, character.only = TRUE))

## setup continuous domain
set.seed(413206)
x <- seq(0, 10, length = 200)
grid.pts <- expand.grid(x, x)

## set up matern params, also set param priors to be used in modeling
sp.alpha <- 2
sp.kappa <- 0.5
sp.var   <- 0.5
gp.int  <- -2
# prior on spde parameters: c(a, b, c, d), where
# P(sp.range < a) = b
# P(sp.sigma > c) = d
matern.pri <- c(10, .95, 1., .05) ## a, b, c, d
# mean and sd for normal prior on fixed effects (alpha and betas)
alpha.pri <- c(0, 3) ## N(mean, sd)

## sample from matern RF on our grid
model <- RMmatern(nu    = sp.alpha - 1, ## from INLA book
                  scale = sqrt(2 * (sp.alpha - 1)) / sp.kappa,
                  var   = 1)
true.gp <- RFsimulate(model, x = x, y = x, n =1, spConform = FALSE)

## insert into a raster
gp.rast <- raster(nrows=length(x), ncols=length(x),
                  xmn=0, xmx=10, ymn=0, ymx=10,
                  vals=(true.gp + gp.int))

## define cluster locations and sample size at each
n.clust <- 500
clust.mean.ss <- 35
dat <- data.table(x = runif(n.clust, min = min(x), max = max(x)),
                  y = runif(n.clust, min = min(x), max = max(x)),
                  n = rpois(n.clust, clust.mean.ss)
)

## extract value of raster at cluster locs and logit transform
##   to binom probs
dat[, latent.truth := raster::extract(x = gp.rast, y = cbind(x, y))]
dat[, p.truth := plogis(latent.truth)]

## sample binomial data
dat[, obs := rbinom(n = .N, size = n, p = p.truth)]

## make SPDE triangulation mesh over our domain
mesh.s <- inla.mesh.2d(loc.domain = grid.pts,
                       max.e = c(0.25, 5))
## check number of vertices
mesh.s[['n']]

## plot true latent field, the observed/empirical binom probs at
##   cluster locs, and the mesh
par(mfrow = c(3, 1))
plot(gp.rast, maxpixels = length(x) ^ 2,
     xlim = range(x), ylim = range(x), main = 'latent truth')
fields::quilt.plot(dat[, x], dat[, y], dat[, obs] / dat[, n],
                   main = 'empirical binom probs')
plot(mesh.s)
polygon(x = c(0, 0, 10, 10, 0), y = c(0, 10, 10, 0, 0),
        col = NA, border = 2, lwd = 5)

## make the SPDE objects (including prec components)
spde <- inla.spde2.pcmatern(mesh = mesh.s, alpha = 2,
                            prior.range = matern.pri[1:2],
                            prior.sigma = matern.pri[3:4])

## make projector matrices to:
## 1) project data to mesh
## 2) project mesh to raster grid
A.proj <- inla.spde.make.A(mesh = mesh.s,
                           loc = dat[, as.matrix(x, y)])
A.pred <- inla.spde.make.A(mesh = mesh.s,
                           loc = as.matrix(grid.pts),
                           group = 1)
Continuous GP modeling with SPDE in R-INLA
## prep inputs for INLA
design_matrix <- data.frame(int = rep(1, nrow(dat)))
stack.obs <- inla.stack(tag='est',
                        data=list(Y = dat$obs, ## response
                                  N = dat$n), ## binom trials
                        A=list(A.proj, ## A.proj for space
                               1),     ## 1 for design.mat
                        effects=list(
                          space = 1:mesh.s[['n']],
                          design_matrix))

## define the INLA model
formula <- formula(Y ~ -1 + int + f(space, model = spde))

## run INLA
i.fit <- inla(formula,
              data = inla.stack.data(stack.obs),
              control.predictor = list(A = inla.stack.A(stack.obs),
                                       compute = FALSE),
              control.fixed = list(expand.factor.strategy = 'inla',
                                   prec = list(default = 1 / alpha.pri[2] ^ 2)),
              control.inla = list(strategy = 'simplified.laplace',
                                  int.strategy = 'ccd'),
              control.compute=list(config = TRUE),
              family = 'binomial',
              Ntrials = N,
              verbose = FALSE,
              keep = FALSE)

## take draws from inla
i.draws <- inla.posterior.sample(n = 500, i.fit,
                                 use.improved.mean = TRUE,
                                 skew.corr = TRUE)

## summarize the draws
par_names <- rownames(i.draws[[1]][['latent']])
s_idx <- grep('^space.*', par_names)
a_idx <- which(!c(1:length(par_names)) %in%
                 grep('^space.*|Predictor|clust.id', par_names))

# project from mesh to raster, add intercept
pred_s <- sapply(i.draws, function (x) x[['latent']][s_idx])
pred_inla <- as.matrix(A.pred %*% pred_s)
alpha_inla_draws <- sapply(i.draws, function (x) x[['latent']][a_idx])
pred_inla <- sweep(pred_inla, 2, alpha_inla_draws, '+')


## find the median and sd across draws, as well as 90% intervals
summ_inla <- cbind(median = (apply(pred_inla, 1, median)),
                   sd     = (apply(pred_inla, 1, sd)),
                   lower = (apply(pred_inla, 1, quantile, .05)),
                   upper = (apply(pred_inla, 1, quantile, .95)))

## make summary rasters
ras_med_inla <- ras_sdv_inla <- ras_lower_inla <-
  ras_upper_inla <- ras_inInt_inla <- gp.rast
values(ras_med_inla)   <- summ_inla[, 1]
values(ras_sdv_inla)   <- summ_inla[, 2]
values(ras_lower_inla) <- summ_inla[, 3]
values(ras_upper_inla) <- summ_inla[, 4]
values(ras_inInt_inla) <- 0
ras_inInt_inla[gp.rast < ras_lower_inla | ras_upper_inla < gp.rast] <- 1

## plot truth, pixels falling within/without the 90% interval,
##  post. median, and post sd

# set the range for the truth and median
rast.zrange <- range(c(values(gp.rast), values(ras_med_inla)), na.rm = T)

# plot
par(mfrow = c(2, 2))
plot(gp.rast, main = 'Truth', zlim = rast.zrange, col = (viridis(100)))
points(dat[, .(x, y)])
plot(ras_inInt_inla, main = 'Pixels where 90% CIs did not cover Truth')
points(dat[, .(x, y)])
plot(ras_med_inla, main = 'INLA Posterior Median',
     zlim = rast.zrange, col = (viridis(100)))
points(dat[, .(x, y)])
plot(ras_sdv_inla, main = 'INLA Posterior Standard Deviation')
points(dat[, .(x, y)])

tmb_spde <-
  "// include libraries
#include <TMB.hpp>
#include <Eigen/Sparse>
#include <vector>
using namespace density;
using Eigen::SparseMatrix;

// helper function for detecting NAs in the data supplied from R
template<class Type>
bool isNA(Type x){
  return R_IsNA(asDouble(x));
}

// helper function to make sparse SPDE precision matrix
// Inputs:
//    logkappa: log(kappa) parameter value
//    logtau: log(tau) parameter value
//    M0, M1, M2: these sparse matrices are output from:
//     R::INLA::inla.spde2.matern()$param.inla$M*
template<class Type>
SparseMatrix<Type> spde_Q(Type logkappa, Type logtau, SparseMatrix<Type> M0,
                          SparseMatrix<Type> M1, SparseMatrix<Type> M2) {
    SparseMatrix<Type> Q;
    Type kappa2 = exp(2. * logkappa);
    Type kappa4 = kappa2*kappa2;
    Q = pow(exp(logtau), 2.)  * (kappa4*M0 + Type(2.0)*kappa2*M1 + M2);
    return Q;
}

// helper function to use the same penalized complexity prior on
//  matern params that is used in INLA

template<class Type>
Type dPCPriSPDE(Type logtau, Type logkappa,
                Type matern_par_a, Type matern_par_b,
                Type matern_par_c, Type matern_par_d,
                //vector<Type> matern_pri(4),
                int give_log=0)
{

  // matern_pri = c(a, b, c, d): P(range < a) = b; P(sigma > c) = d

  Type penalty; // prior contribution to jnll

  Type d = 2.;  // dimension
  Type lambda1 = -log(matern_par_b) * pow(matern_par_a, d/2.);
  Type lambda2 = -log(matern_par_d) / matern_par_c;
  Type range   = sqrt(8.0) / exp(logkappa);
  Type sigma   = 1.0 / sqrt(4.0 * 3.14159265359 * exp(2.0 * logtau) *
                            exp(2.0 * logkappa));

  penalty = (-d/2. - 1.) * log(range) - lambda1 * pow(range, -d/2.) -
               lambda2 * sigma;
  // Note: (rho, sigma) --> (x=log kappa, y=log tau) -->
  //  transforms: rho = sqrt(8)/e^x & sigma = 1/(sqrt(4pi)*e^x*e^y)
  //  --> Jacobian: |J| propto e^(-y -2x)
  Type jacobian = - logtau - 2.0*logkappa;
  penalty += jacobian;

  if(give_log)return penalty; else return exp(penalty);
}

///////////////////////////
// the main function     //
// to calculate the jnll //
///////////////////////////
template<class Type>
Type objective_function<Type>::operator() ()
{

  // ~~~~~~~~~------------------------------------------------------~~
  // FIRST, we define params/values/data that will be passed in from R
  // ~~~~~~~~~~~------------------------------------------------------

  // normalization flag - used for speed-up
  DATA_INTEGER( flag ); // flag == 0 => no data contribution added to jnll

  // Indices
  DATA_INTEGER( num_i );   // Number of data points in space
  DATA_INTEGER( num_s );   // Number of mesh points in space mesh

  // Data (all except for X_ij is a vector of length num_i)
  DATA_VECTOR( y_i );   // obs per binomial experiment at point i (clust)
  DATA_VECTOR( n_i );   // Trials per cluster
  DATA_MATRIX( X_alpha );  // 'design matrix' for just int

  // SPDE objects
  DATA_SPARSE_MATRIX( M0 );
  DATA_SPARSE_MATRIX( M1 );
  DATA_SPARSE_MATRIX( M2 );
  DATA_SPARSE_MATRIX( Aproj );

  // Options
  DATA_VECTOR( options );
  // options[0] == 1 : use normalization trick
  // options[1] == 1 : adreport transformed params

  // Prior specifications
  DATA_VECTOR( alpha_pri );
  DATA_VECTOR( matern_pri);
  // matern_pri = c(a, b, c, d): P(range < a) = b; P(sigma > c) = d
  Type matern_par_a = matern_pri[0]; // range limit:    rho0
  Type matern_par_b = matern_pri[1]; // range prob:     alpha_rho
  Type matern_par_c = matern_pri[2]; // field sd limit: sigma0
  Type matern_par_d = matern_pri[3]; // field sd prob:  alpha_sigma

  // Fixed effects
  PARAMETER( alpha ); // Intercept
  // Log of INLA tau param (precision of space covariance matrix)
  PARAMETER( log_tau );
  // Log of INLA kappa (related to spatial correlation and range)
  PARAMETER( log_kappa );

  // Random effects for each spatial mesh vertex
  PARAMETER_VECTOR( Epsilon_s );

  // ~~~~~~~~~------------------------------------------------~~
  // SECOND, we define all other objects that we need internally
  // ~~~~~~~~~------------------------------------------------~~

  // objective function -- joint negative log-likelihood
  Type jnll = 0;

  // Make spatial precision matrix
  SparseMatrix<Type> Q_ss = spde_Q(log_kappa, log_tau, M0, M1, M2);

  // Transform some of our parameters
  Type sp_range = sqrt(8.0) / exp(log_kappa);
  Type sp_sigma = 1.0 / sqrt(4.0 * 3.14159265359 *
                  exp(2.0 * log_tau) * exp(2.0 * log_kappa));

  // Define objects for derived values
  vector<Type> fe_i(num_i); // main effect: alpha
  // Logit estimated prob for each cluster i
  vector<Type> latent_field_i(num_i);
  // value of gmrf at data points
  vector<Type> projepsilon_i(num_i);

  // fixed effects is just alpha in this example
  fe_i = X_alpha * Type(alpha); // initialize

  // Project GP approx from mesh points to data points
  projepsilon_i = Aproj * Epsilon_s.matrix();

  // ~~~~~~~~~------------------------------------------------~~-
  // THIRD, we calculate the contribution to the likelihood from:
  // 1) priors
  // 2) GP field
  // 3) data
  // ~~~~~~~~~------------------------------------------------~~-

  /////////
  // (1) //
  /////////
  // the random effects. we do this first so to do the
  //   normalization outside of every optimization step
  // NOTE: likelihoods from namespace 'density' already return NEGATIVE
  //       log-liks so we add other likelihoods return positive log-liks
  if(options[0] == 1){
    // then we are not calculating the normalizing constant in the inner opt
    // that norm constant means taking an expensive determinant of Q_ss
    jnll += GMRF(Q_ss, false)(Epsilon_s);
    // return without data ll contrib to avoid unneccesary log(det(Q)) calcs
    if (flag == 0) return jnll;
  }else{
    jnll += GMRF(Q_ss)(Epsilon_s);
  }

  /////////
  // (2) //
  /////////
  // Prior contributions to joint likelihood (if options[1]==1)

  // add in priors for spde gp
  jnll -= dPCPriSPDE(log_tau, log_kappa,
                     matern_par_a, matern_par_b, matern_par_c, matern_par_d,
                     true);

  // prior for intercept
  jnll -= dnorm(alpha, alpha_pri[0], alpha_pri[1], true); // N(mean, sd)

  /////////
  // (3) //
  /////////
  // jnll contribution from each datapoint i

  for (int i = 0; i < num_i; i++){

    // latent field estimate at each obs
    latent_field_i(i) = fe_i(i) + projepsilon_i(i);

    // and add data contribution to jnll
    if(!isNA(y_i(i))){

     // Uses the dbinom_robust function, which takes the logit probability
      	jnll -= dbinom_robust( y_i(i), n_i(i), latent_field_i(i), true );

    } // !isNA

  } // for( i )


  // ~~~~~~~~~~~
  // ADREPORT: used to return estimates and cov for transforms?
  // ~~~~~~~~~~~
  if(options[1]==1){
    ADREPORT(sp_range);
    ADREPORT(sp_sigma);
  }

  return jnll;

}"
  
  ## write model to file, compile, and load it into R
  dir.create('TMB_spde_example')
  write(tmb_spde,file="TMB_spde_example/tmb_spde.cpp")
  compile( "TMB_spde_example/tmb_spde.cpp")
  dyn.load( dynlib("TMB_spde_example/tmb_spde") )
  
  ## prep inputs for TMB
  data_full <- list(num_i = nrow(dat),  # Total number of observations
                    num_s = mesh.s[['n']], # num. of vertices in SPDE mesh
                    y_i   = dat[, obs],# num. of pos. obs in the cluster
                    n_i   = dat[, n],  # num. of exposures in the cluster
                    X_alpha  = matrix(1, nrow = nrow(dat), ncol = 1),# des.mat
                    M0    = spde[['param.inla']][['M0']], # SPDE sparse matrix
                    M1    = spde[['param.inla']][['M1']], # SPDE sparse matrix
                    M2    = spde[['param.inla']][['M2']], # SPDE sparse matrix
                    Aproj = A.proj,             # Projection matrix
                    options = c(1, ## if 1, use normalization trick
                                1), ## if 1, run adreport
                    # normalization flag.
                    flag = 1,
                    alpha_pri = alpha.pri, ## normal
                    matern_pri = matern.pri
  )
  
  ## Specify starting values for TMB params
  tmb_params <- list(alpha = 0.0, # intercept
                     log_tau = 0, # Log inverse of tau (Epsilon)
                     log_kappa = 0, # Matern range parameter
                     Epsilon_s = rep(0, mesh.s[['n']]) # RE on mesh vertices
  )
  
  ## make a list of things that are random effects
  rand_effs <- c('Epsilon_s')
  
  ## make the autodiff generated liklihood func & gradient
  obj <- MakeADFun(data=data_full,
                   parameters=tmb_params,
                   random=rand_effs,
                   hessian=TRUE,
                   DLL='tmb_spde')
  
  ## we can normalize the GMRF outside of the nested optimization,
  ## avoiding unnecessary and expensive cholesky operations.
  obj <- normalize(obj, flag="flag", value = 0)
  
  ## run TMB
  opt0 <- nlminb(start       =    obj[['par']],
                 objective   =    obj[['fn']],
                 gradient    =    obj[['gr']],
                 lower = rep(-10, length(obj[['par']])),
                 upper = rep( 10, length(obj[['par']])),
                 control     =    list(trace=1))
  
  ## Get standard errors
  SD0 <- TMB::sdreport(obj, getJointPrecision=TRUE,
                       bias.correct = TRUE,
                       bias.correct.control = list(sd = TRUE))
  ## summary(SD0, 'report')
  
  ## take samples from fitted model
  mu <- c(SD0$par.fixed,SD0$par.random)
  
  ## simulate draws
  rmvnorm_prec <- function(mu, chol_prec, n.sims) {
    z <- matrix(rnorm(length(mu) * n.sims), ncol=n.sims)
    L <- chol_prec #Cholesky(prec, super=TRUE)
    z <- Matrix::solve(L, z, system = "Lt") ## z = Lt^-1 %*% z
    z <- Matrix::solve(L, z, system = "Pt") ## z = Pt    %*% z
    z <- as.matrix(z)
    mu + z
  }
  
  L <- Cholesky(SD0[['jointPrecision']], super = T)
  t.draws <- rmvnorm_prec(mu = mu , chol_prec = L, n.sims = 500)
  
  ## summarize the draws
  parnames <- c(names(SD0[['par.fixed']]), names(SD0[['par.random']]))
  epsilon_tmb_draws  <- t.draws[parnames == 'Epsilon_s',]
  alpha_tmb_draws    <- matrix(t.draws[parnames == 'alpha',], nrow = 1)
  
  # project from mesh to raster, add intercept
  pred_tmb <- as.matrix(A.pred %*% epsilon_tmb_draws)
  pred_tmb <- sweep(pred_tmb, 2, alpha_tmb_draws, '+')
  
  ## find the median and sd across draws, as well as 90% intervals
  summ_tmb <- cbind(median = (apply(pred_tmb, 1, median)),
                    sd     = (apply(pred_tmb, 1, sd)),
                    lower = (apply(pred_tmb, 1, quantile, .05)),
                    upper = (apply(pred_tmb, 1, quantile, .95)))
  
  ## make summary rasters
  ras_med_tmb <- ras_sdv_tmb <- ras_lower_tmb <-
    ras_upper_tmb <- ras_inInt_tmb <- gp.rast
  values(ras_med_tmb)   <- summ_tmb[, 1]
  values(ras_sdv_tmb)   <- summ_tmb[, 2]
  values(ras_lower_tmb) <- summ_tmb[, 3]
  values(ras_upper_tmb) <- summ_tmb[, 4]
  values(ras_inInt_tmb) <- 0
  ras_inInt_tmb[gp.rast < ras_lower_tmb | ras_upper_tmb < gp.rast] <- 1
  
  ## plot truth, pixels falling within/without the 90% interval,
  ##  post. median, and post sd
  
  # set the range for the truth and median
  rast.zrange <- range(c(values(gp.rast), values(ras_med_tmb)), na.rm = T)
  
  # plot tmb
  png(file='figures/example_tmb.png', width=9, height=9, units='in', res=300)
  par(mfrow = c(2, 2))
  plot(gp.rast, main = 'Truth', zlim = rast.zrange, col = (viridis(100)))
  points(dat[, .(x, y)])
  plot(ras_inInt_tmb, main = 'Pixels where 90% CIs did not cover Truth')
  points(dat[, .(x, y)])
  plot(ras_med_tmb, main = 'TMB Posterior Median',
       zlim = rast.zrange, col = (viridis(100)))
  points(dat[, .(x, y)])
  plot(ras_sdv_tmb, main='TMB Posterior Standard Deviation')
  points(dat[, .(x, y)])
  dev.off()
  
  
  ## compare INLA and TMB meds and stdevs
  med.zrange <- range(c(values(ras_med_tmb), values(ras_med_inla)), na.rm = T)
  sdv.zrange <- range(c(values(ras_sdv_tmb), values(ras_sdv_inla)), na.rm = T)
  
  png(file='figures/example_tmb_v_inla.png', width=9, height=9, units='in', res=300)
  par(mfrow = c(2, 2))
  plot(ras_med_inla, main = 'INLA Posterior Median',
       zlim = med.zrange, col = (viridis(100)))
  points(dat[, .(x, y)])
  plot(ras_sdv_inla, main = 'INLA Posterior Standard Deviation',
       zlim = sdv.zrange)
  points(dat[, .(x, y)])
  plot(ras_med_tmb, main = 'TMB Posterior Median',
       zlim = med.zrange, col = (viridis(100)))
  points(dat[, .(x, y)])
  plot(ras_sdv_tmb, main = 'TMB Posterior Standard Deviation',
       zlim = sdv.zrange)
  points(dat[, .(x, y)])
  dev.off()