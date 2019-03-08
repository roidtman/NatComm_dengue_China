rm(list = ls())

setwd("~/Dropbox/DENV_China_Drivers/code")
# setwd("~/Private/DENV_China_Drivers/code")


##### set up libraries
library(mgcv)
library(coda)
library(mvtnorm)
library(fda)
library(lubridate)

# =============================================================== # 
###### load & set up data
source('runSerialinterval.R')

mosq = read.csv('../data/guangzhou_vector_05-15.csv')
GZ.data = read.csv('../data/guangzhou_case_weather_05-15.csv')

# load mosquito mcmc output
for(load.count in 1:3){
  eval(parse(text = paste0(
    'load("../output/mosquito_mcmc_', load.count, '.RData")'
  )))
  name.chain = paste('param.chain_', load.count, sep = '')
  assign(name.chain, param.chain)
}
rm(param.chain)
mosq.chain_noTrans = rbind(param.chain_1, param.chain_2, param.chain_3)
mosq.median = sapply(1:ncol(mosq.chain_noTrans), function(mm) median(mosq.chain_noTrans[,mm]))

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
#####

# =============================================================== # 
##### set values
# general inputs
years = 2005:2015

# mosquito spline input
k = 3
set.order = 3

# biv spline input
k.biv = 3
set.order.biv = 3

# mosquito spline basis
mosq_basis = create.bspline.basis(
  rangeval = c(1, nrow(GZ.data)),
  nbasis = k * length(years),
  norder = set.order
)

# inputs for mosquito likelihood fxn
ind = list()
for(ii in 1:nrow(mosq)){
  ind[[ii]] = which(GZ.data$year == mosq$Year[ii] & month(GZ.data$date) == mosq$Month[ii])
}

initial.mosq = rep(0, 36)
#####

# =============================================================== # 
##### set up functions
generate_spline = function(y.in, dat = GZ.data, bs.in){
  fd.params = y.in
  fd.bs = fd(
    coef = fd.params,
    basisobj = bs.in
  )
  spline.tmp = predict(fd.bs, 1:nrow(dat))
  return(spline.tmp)
}

nll = function(y.in, dat = GZ.data, mosq.dat = mosq, bs.in = mosq_basis, browse = F){
  if(browse) browser()
  
  mosq.in = y.in[1:33]
  sigma.moi.in = y.in[34]
  sigma.bi.in = y.in[35]
  mu.bi.in = y.in[36]
  
  fd.bs = fd(
    coef = mosq.in,
    basisobj = bs.in
  )
  spline.tmp = predict(fd.bs, 1:nrow(dat))
  mean.mosq.spline.in = sapply(1:nrow(mosq.dat), function(ii) mean(spline.tmp[ind[[ii]]]))
  
  LL.moi = dnorm(MOI, mean = exp(mean.mosq.spline.in), sd = exp(sigma.moi.in), log = TRUE)
  LL.bi = dnorm(BI, mean = exp(mu.bi.in) * exp(mean.mosq.spline.in), sd = exp(sigma.bi.in), log = TRUE)
  LL = sum(LL.moi) + sum(LL.bi)
  
  return(-LL)
}

# =============================================================== # 
##### optim

ss = sample(1:nrow(mosq.chain_noTrans), 1000, replace = F)
init = mosq.chain_noTrans[ss, ]

f = '../output/mosqOnly.optim.RData'
if(file.exists(f)){
  load(f)
}else{
  out.list = list()
  length(out.list) = nrow(init)
  out.LL = rep(NA, nrow(init))
  for(ii in 1:nrow(init)){
    par.tmp = init[ii, ]
    
    out = optim(par = par.tmp,
                fn = nll,
                method = 'Nelder-Mead',
                control = list(maxit = 10000))
    out.list[[ii]] = out$par
    out.LL[ii] = out$value
  }
  save(out.LL, out.list, file = f)
}

MLE = out.list[[which(out.LL == min(out.LL))]]





