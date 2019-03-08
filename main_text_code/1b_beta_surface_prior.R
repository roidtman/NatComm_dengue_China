
##### load libraries 
require(fda)
require(parallel)

# =============================================================== # 
##### loading and setting up data 
allR0 = read.csv('../data/siraj_et_al_allRo_with_uncertainty.csv')

# biv spline input
k.biv = 3
reps.biv = 2
set.order.biv = 3
# =============================================================== # 
###### setting up functions 

# fxn to generate the bivariate spline surface
# inputs: temp.dat = temp range vector, lag.dat = lag range vector, 
# coef.prop = proposed paramters for temp and lag, n.knots = # of knots, 
# poly.order = order of spline, r = # of repitions in the coef mat
# outputs: fd.biv = covariance matrix corresponding to temp.range and lag.range
generate_biv_spline = function(coef.prop, dat.range, lag.dat, 
                               n.knots = k.biv, poly.order = set.order.biv, r = reps.biv, browse = F){
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

# optim function
# inputs: param = 
# outputs: sum of squared error
optim_surface = function(param, mosq.in = mosq.range, lag.mosq.in = l, temp.in = temp.range, lag.temp.in = l, browse = F){
  if(browse) browser()
  
  mosq.param = c(rep(param[1:3], 3))
  temp.param = c(rep(param[4:6], 3))
  
  mosq.biv = generate_biv_spline(coef.prop = mosq.param,
                                 dat.range = mosq.in, 
                                 lag.dat = lag.mosq.in)
  temp.biv = generate_biv_spline(coef.prop = temp.param,
                                 dat.range = temp.in,
                                 lag.dat = lag.temp.in)
  
  beta.tmp = matrix(rep(1, (nrow(mosq.biv) * nrow(temp.biv))), nrow = nrow(mosq.biv), ncol = nrow(temp.biv))
  
  for(ii in 1:nrow(beta.tmp)){ # loop across mosq indices (nrow(beta.tmp))
    lag.eff.mosq = sum(mosq.biv[ii, ])
    
    for(tt in 1:ncol(beta.tmp)){ # loop across temp indices (ncol(beta.tmp))
      lag.eff.temp = sum(temp.biv[tt, ])
      beta.tmp[ii, tt] = exp(sum(lag.eff.temp, lag.eff.mosq))
    }
  }
  
  SSq = sum((beta.idealized - beta.tmp) ^ 2)
  return(SSq)
}

# 
# function to create random inputs for optim_surface
# input: n = numCores
# output: initList = list of randomly generated initial starting parameters
rand_init = function(n){
  initList = list()
  length(initList) = n
  for(ii in 1:n){
    initList[[ii]] = rnorm(6, mean = 0, sd = 0.15)
  }
  return(initList)
}

# function over our optimize fxn in order to parallelize
# inputs: inits = randomly generated inputs
optim_wrap = function(inits){
  print(inits)
  optim(par = inits, fn = optim_surface, method = 'Nelder-Mead',
        hessian = TRUE, control = list(maxit = 10000))
}

# =============================================================== # 
###### setting up idealized surface
l = 1:49
allR0 = as.matrix(allR0)
R0.arr = array(NA, c(nrow(allR0)/1000, ncol(allR0), 1000))
t = -49
for(ii in 1:1000){
  for(kk in 1:50){
    if(kk == 1){
      t = t + 50
      R0.arr[kk,, ii] = allR0[t, ]
    }else{
      R0.arr[kk,, ii] = allR0[(t + kk - 1), ]
    }
  }
}; rm(ii,kk)

mosq.range = seq(0, 3, by = 0.1)
mosq.range.round = round(mosq.range, 2)
temp.range = seq(4, 36, by = 0.1)
temp.range.round = round(temp.range, 2)

# want to use  Amir's temperature relationship @ max m
R0.mat = matrix(NA, nrow = 1000, ncol = length(temp.range))
for(ii in 1:1000){
  R0.mat[ii, ] = R0.arr[50,, ii]
}

# determine mean max R0
max.R0 = sapply(1:nrow(R0.mat), function(kk) max(R0.mat[kk,]))
mean.max.R0 = mean(max.R0)

# we are following estimates from Reiner et al. 2015 and setting maximum R0 = ~2.5
mult = rep(NA, 1000)
for(ii in 1:1000){
  max.R0.idealized = abs(rnorm(1, mean = 2.5, sd = 1))
  mult[ii] = (max.R0.idealized/3) * (1 / mean.max.R0)
}

R0.arr.new = array(NA, c(length(mosq.range), length(temp.range), 1000))
for(ii in 1:1000){
  R0.reduced = mult[ii] * R0.mat[ii, ]
  R0.arr.new[,,ii] = t(sapply(1:length(mosq.range), function(kk) R0.reduced * mosq.range[kk]))
}

###### setting up optimization to run in parallel 
optim.out.full = list()
length(optim.out.full) = dim(R0.arr.new)[3]
for(ii in 1:dim(R0.arr)[3]){
  beta.idealized = R0.arr.new[,,ii]
  numCores = 3
  initial.param = rand_init(numCores)
  optim.out = mclapply(initial.param, optim_wrap, mc.cores = numCores)
  SSE = rep(NA, numCores)
  for(kk in 1:numCores){
    SSE[kk] = optim.out[[kk]]$value
  }
  ind.min = which(SSE == min(SSE))
  optim.out.full[[ii]] = list('par' = optim.out[[ind.min]]$par, 'SSE' = optim.out[[ind.min]]$value)
  save(optim.out.full, file = '../output/optim.surface.RData')
}


# =============================================================== #
# determining median 
temp.optimum = median(sapply(1:nrow(R0.mat), function(i) temp.range[which(R0.mat[i, ] == max(R0.mat[i, ]))]))

# determining mean and 95% CI
temp.optimum.range = rep(NA, nrow(R0.mat))
for(ii in 1:nrow(R0.mat)){
  temp.optimum.range[ii] = temp.range[which(R0.mat[ii,] == max(R0.mat[ii, ]))]
}



# =============================================================== #
####### creating prior distributions
if(!require(mvtnorm)){install.packages('mvtnorm'); library(mvtnorm)}
source('runFunctions_bayes.R')

reps = 1000
param.optim = matrix(NA, nrow = reps, ncol = 6)

for(ii in 1:reps){
  param.optim[ii, ] = optim.out.full[[ii]]$par
}
sigma.optim = cov(param.optim)

param.mean = colMeans(param.optim)
param.mean = c(rep(param.mean[1:3], 3), rep(param.mean[4:6], 3))

mat = matrix(NA, nrow = length(param.mean[1:18]), ncol = length(param.mean[1:18]))
colnames(mat) = c('m1', 'm2', 'm3', 'm4', 'm5', 'm6', 'm7', 'm8', 'm9', 't1', 't2', 't3',
                  't4', 't5', 't6', 't7', 't8', 't9')
rownames(mat) = c('m1', 'm2', 'm3', 'm4', 'm5', 'm6', 'm7', 'm8', 'm9', 't1', 't2', 't3',
                  't4', 't5', 't6', 't7', 't8', 't9')

m_m = sigma.optim[1:3, 1:3]
m_t = sigma.optim[4:6, 1:3]
t_m = sigma.optim[1:3, 4:6]
t_t = sigma.optim[4:6, 4:6]

mat[1:3, 1:3] = m_m
mat[4:6, 4:6] = m_m
mat[7:9, 7:9] = m_m
mat[10:12, 10:12] = t_t
mat[13:15, 13:15] = t_t
mat[16:18, 16:18] = t_t

mat[10:12, 1:3] = m_t
mat[13:15, 4:6] = m_t
mat[16:18, 7:9] = m_t

mat[1:3, 10:12] = t_m
mat[4:6, 13:15] = t_m
mat[7:9, 16:18] = t_m

mat[which(is.na(mat))]=0

sigma.optim = mat

save(param.mean, sigma.optim, file = '../output/prior_beta.RData')
