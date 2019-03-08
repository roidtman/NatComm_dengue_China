
count = 1


##### set up libraries
if(!require(BayesianTools)){install.packages('BayesianTools'); library(BayesianTools)}
if(!require(mgcv)){install.packages('mgcv'); library(mgcv)}
if(!require(coda)){install.packages('coda'); library(coda)}
if(!require(mvtnorm)){install.packages('mvtnorm'); library(mvtnorm)}
if(!require(fda)){install.packages('fda'); library(fda)}
if(!require(lubridate)){install.packages('lubridate'); library(lubridate)}
if(!require(VGAM)){install.packages('VGAM'); library(VGAM)}

# =============================================================== # 
###### load & set up data
source('runFunctions_bayes.R')
source('runSerialinterval.R')

mosq = read.csv('../data/guangzhou_vector_05-15.csv')
GZ.data = read.csv('../data/guangzhou_case_weather_05-15.csv')

# load mosquito MLE 
load('../output/mosqOnly.optim.RData')
mosq.median = out.list[[which(out.LL == min(out.LL))]]
generate_spline = function(y.in, dat = GZ.data, bs.in, browse = F){
  if(browse) browser()
  fd.params = y.in
  fd.bs = fd(
    coef = fd.params,
    basisobj = bs.in
  )
  spline.tmp = predict(fd.bs, 1:nrow(dat))
  return(spline.tmp)
}
mosq_basis = create.bspline.basis(
  rangeval = c(1, nrow(GZ.data)),
  nbasis = 3 * 11,
  norder = 3
)
mosq.spline = exp(generate_spline(mosq.median[1:33], bs.in = mosq_basis))

# for MOI replace NA values w/ means from other years
ind.NA = which(is.na(mosq$MOI_pnas))
for(ii in 1:length(ind.NA)){
  mosq$MOI_pnas[ind.NA[ii]] = mean(mosq$MOI_pnas[which(mosq$Month==mosq$Month[ind.NA[ii]])], na.rm = T)
}
# set up MOI object
MOI = mosq$MOI_pnas

# for BI replace NA values w/ means from other years
ind.NA = which(is.na(mosq$BI))
for(ii in 1:length(ind.NA)){
  mosq$BI[ind.NA[ii]] = mean(mosq$BI[which(mosq$Month==mosq$Month[ind.NA[ii]])], na.rm = T)
}
# set up BI object
BI = mosq$BI

# load prior distribution
load('../output/prior_beta.RData')
sigma.mean = sigma.optim
#####

# =============================================================== # 
##### set values
# Guangzhou population as of 2016
GZ.pop = 14.04e6

# general inputs
years = 2005:2015

# mosquito spline input
k = 3
set.order = 3

# biv spline input
k.biv = 3
set.order.biv = 3

# beta0 spline basis 
b0_basis = create.bspline.basis(
  rangeval = c(1, nrow(GZ.data)),
  nbasis = 3 * length(years),
  norder = set.order
)

# inputs for bivariate temp spline
# temp.range = seq(floor(min(GZ.data$T_mean)), ceiling(max(GZ.data$T_mean)), by = 0.1)
temp.range = seq(4, 36, by = 0.1)
temp.range.round = signif(temp.range, 3)
l = 1:49

align_index_T = function(dat = GZ.data, browse = F){
  if(browse) browser()
  temp.ind = matrix(NA, nrow = nrow(dat), ncol = length(l))
  
  for(ii in (length(l) + 1):nrow(dat)){
    t = GZ.data$T_mean[ii - l]
    tmp = unlist(sapply(1:length(l), function(kk) which(temp.range.round == t[kk])))
    if(length(tmp) == length(l)){
      temp.ind[ii, ] = tmp
    }else{
      temp.ind[ii, ] = c(tmp, rep(NA, (length(l) - length(tmp))))
    }
  }
  
  return(temp.ind)
}
temp.ind = align_index_T(browse = F)

# inputs for bivariate mosq spline
mosq.range = seq(0, 3, by = 0.01)
mosq.range.round = round(mosq.range, 2)
align_index_m = function(dat = GZ.data, browse = F){
  if(browse) browser()
  mosq.ind = matrix(NA, nrow = nrow(dat), ncol = length(l))
  mosq.spline = round(mosq.spline, 2)
  
  for(ii in (length(l) + 1):nrow(dat)){
    m = mosq.spline[ii - l]
    tmp = unlist(sapply(1:length(l), function(kk) which(mosq.range.round == m[kk])))
    if(length(tmp) == length(l)){
      mosq.ind[ii, ] = tmp
    }else{
      mosq.ind[ii, ] = c(tmp, rep(NA, (length(l) - length(tmp))))
    }
  }
  return(mosq.ind)
}
mosq.ind = align_index_m(browse = F)

# days which has non-zero cases
ind.cases = which(GZ.data$I.old != 0)
I.new = GZ.data$local.cases
I.old = GZ.data$I.old

# inputs for MCMC - initial values
# ==== # mosquito spline & weather initial values
y.in = mosq.median[1:(k* length(years))]
initial.biv = param.mean[1:18]
initial.b0 = rep(0, (length(years) * 3))
initial.param = c(initial.biv, initial.b0)

# =============================================================== # 
##### set up functions
# create daily spline for mosquito data
# input: y.in and dat
# y.in takes mean MOI
# output: mosq.spline w/ set.order number of knots
generate_spline = function(y.in, dat = GZ.data, bs.in, browse = F){
  if(browse) browser()
  fd.params = y.in
  fd.bs = fd(
    coef = fd.params,
    basisobj = bs.in
  )
  spline.tmp = predict(fd.bs, 1:nrow(dat))
  return(spline.tmp)
}

# fxn to generate the bivariate spline surface
# inputs: temp.dat = temp range vector, lag.dat = lag range vector, 
# coef.prop = proposed paramters for temp and lag, n.knots = # of knots, 
# poly.order = order of spline, r = # of repitions in the coef mat
# outputs: fd.biv = covariance matrix corresponding to temp.range and lag.range
generate_biv_spline = function(coef.prop, dat.range, lag.dat = l, 
                               n.knots = k.biv, poly.order = set.order.biv, r = 2, browse = F){
  if(browse) browser()
  
  basis.dat = create.bspline.basis(rangeval = c(min(dat.range), max(dat.range)),
                                   nbasis = n.knots, 
                                   norder = poly.order)
  basis.lag = create.bspline.basis(rangeval = c(min(lag.dat), max(lag.dat)),
                                   nbasis = n.knots,
                                   norder = poly.order)
  
  coefs.dummy = array(0, c(n.knots, r, 2))
  fd.dat = fd(coef = coefs.dummy[,,1],
              basisobj = basis.dat)
  fd.lag = fd(coef = coefs.dummy[,,2], 
              basisobj = basis.lag)
  fd.var = var.fd(fdobj1 = fd.dat, 
                  fdobj2 = fd.lag)
  
  fd.var$coefs = matrix(coef.prop, n.knots, n.knots)
  fd.biv = eval.bifd(sevalarg = dat.range,
                     tevalarg = lag.dat,
                     bifd = fd.var)
  return(fd.biv)
}

# sim epidemic function
# simulate epidemic
# input: GZ.data = GZ.data
# n.reps = # of simulations, beta_mat = matrix of betas produced from write_beta
# output: matrix with nrows corresponding to n.reps and ncols as time series
sim_epidemic = function(dat = GZ.data, n.reps = 100, beta.vec, browse = F){
  if(browse) browser()
  
  I.new = matrix(0, nrow = n.reps, ncol = nrow(dat))
  I.old = matrix(0, nrow = n.reps, ncol = 1)
  imp.cases = matrix(rep(GZ.data$imported.cases, n.reps), nrow = n.reps, ncol = nrow(dat), byrow = T)
  SI = matrix(weights.serialinterval)
  
  for(tt in (length(weights.serialinterval) + 1):nrow(dat)){
    A = I.new[, seq(tt-1, tt-length(weights.serialinterval), -1)]
    A.SI = A %*% SI
    
    Imp = imp.cases[, seq(tt-1, tt-length(weights.serialinterval), -1)]
    Imp.SI = Imp %*% SI
    
    I.old = A.SI + Imp.SI
    
    ind.sim = which(I.old > 0)
    if(length(ind.sim) != 0){
      I.new[ind.sim, tt] = rnbinom(length(ind.sim), size = I.old[ind.sim], mu = I.old[ind.sim] * beta.vec[tt])
      # print(tt)
    }
    # if(sum(I.new[, tt] > 1e4) > 0){
    #   browser()
    # }
  }
  return(I.new)
}

# likelihood function
lik_weather = function(sim = I.new.sim, n.reps, dat = GZ.data, browse = F){
  if(browse) browser()
  
  sim.mean = colMeans(sim)
  
  LL.calc = 
    tryCatch({
      sum(dbetabinom.ab(dat$local.cases, GZ.pop, (1 + sim.mean), (1 + GZ.pop - sim.mean), log = T))
    }, error = function(err){
      if(sum(grep('missing value where TRUE/FALSE needed', err$message) > 0)){
        # print(sprintf("The error is:: %s",err))
        LL.calc = NA
      }else{
        stop(err)
      }
    })
  
  # print(paste0('lik_weather=',LL.calc))
  return(LL.calc)
}

# proposal draw from multivariate normal for weather and mosquito 
# parameters w/ hard cutoff if 20 < beta
# input: param.old = current weather parameters and mosquito paramers,
# sigma.in = covariance matrix
# output: randomly generated values for each parameter
calc_prop_beta = function(tt, biv.temp, biv.mosq, b0.spline.tmp, browse = F){
  if(browse) browser()
  
  T.coefs = unlist(sapply(1:length(l), function(ii) biv.temp[temp.ind[tt, ii], l[ii]]))
  
  mosq.coefs = unlist(sapply(1:length(l), function(ii) biv.mosq[mosq.ind[tt, ii], l[ii]]))
  
  param.prop.new = c(mosq.coefs, T.coefs, b0.spline.tmp[tt])
  sum(param.prop.new, na.rm = T)
}

# =============================================================== # 
##### MCMC

ll = function(params, browse = F){
  if(browse) browser()
  
  b0.spline = generate_spline(y.in = params[19:length(params)], bs.in = b0_basis)
  
  mosq.biv = generate_biv_spline(coef.prop = params[1:9], dat.range = mosq.range)
  temp.biv = generate_biv_spline(coef.prop = params[10:18], dat.range = temp.range)
  
  beta.ts = c(rep(0, length(weights.serialinterval)), exp(sapply((length(weights.serialinterval) + 1):nrow(GZ.data), function(tt)
    calc_prop_beta(tt, biv.temp = temp.biv, biv.mosq = mosq.biv, browse = F, b0.spline.tmp = b0.spline))))
  
  I.new.sim = sim_epidemic(n.reps = 100, beta.vec = beta.ts)
  
  LL = lik_weather(sim = I.new.sim)
  if(is.na(LL)){
    LL = -(.Machine$double.xmax)/10
  }
  return(0.001 * (LL))
}

# prior function
p = function(params){
  p1 = dmvnorm(params[1:18], param.mean, makePositiveDefinite(sigma.mean), log = TRUE)
  p2 = dnorm(mean(params[19:51]), 0, 2, log = TRUE)
  return(0.001 * (sum(p1) + sum(p2)))
}
# sampler function
s = function(n1=1, n2=33){
  # browser()
  s1 = rmvnorm(n1, param.mean, makePositiveDefinite(sigma.mean))
  s2 = rnorm(n2, 0, 0.5)
  return(c(s1, s2))
}
priors = createPrior(density = p, best = initial.param, sampler = s)

# bayesian setup for 
bayesianSetup = createBayesianSetup(likelihood = ll, prior = priors)
settings = list(iterations = 10, nrChains = 1, initialParticles = 10000)
out = runMCMC(bayesianSetup = bayesianSetup, settings = settings, sampler = 'SMC')
f = paste0('../output/posterior_SMC_V3_', count, '.RData')
save(out, file = f)
