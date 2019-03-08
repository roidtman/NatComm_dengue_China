###### Main text figures ######
rm(list = ls())

setwd('~/Dropbox/DENV_China_Drivers/code/')
require(lubridate)
require(fda)
require(RColorBrewer)
require(BayesianTools)
# =============================================================== # 
###### load & set up data

source('runSerialinterval.R')

# weather data, imported cases data, mosquito data 
GZ.data = read.csv('../data/Dengue_case_weather_data_2005-2015_20170203.csv')
mosq = read.csv('../data/Guangzhou Vector 2005-2015 - update.csv')

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
# posterior_SMC_V3_precip_20180811_1.RData
for(ii in 1:3){
  f = paste0('../output/posterior_SMC_V3_precip_20180811_', ii, '.RData')
  load(f)
  if(ii == 1){
    param.chain = getSample(out, numSamples = 334)
  }else{
    param.chain = rbind(param.chain, getSample(out, numSamples = 333))
  }
  rm(out)
}

mosq.biv.chain = param.chain[, 1:9]
weather.biv.chain = param.chain[, 10:18]
precip.biv.chain = param.chain[, 19:27]
int.chain = param.chain[, 28:ncol(param.chain)]

# polynomial order of spline
k = 3
set.order = 3
years = 2005:2015

b0_k = 3

# mosquito spline basis
mosq_basis = create.bspline.basis(
  rangeval = c(1, nrow(GZ.data)),
  nbasis = k * length(years),
  norder = set.order
)
# beta0 spline basis 
b0_basis = create.bspline.basis(
  rangeval = c(1, nrow(GZ.data)),
  nbasis = b0_k * length(years),
  norder = set.order
)

# set up MOI & BI objects
MOI = mosq$MOI_pnas
BI = mosq$BI

# biv spline input
k.biv = 3
set.order.biv = 3

# inputs for bivariate temp spline
temp.range = seq(4, 36, by = 0.1)
temp.range.round = round(temp.range, 2)
l = 1:49
align_index = function(dat = GZ.data, browse = F){
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
temp.ind = align_index(browse = F)

# inputs for bivariate mosq spline
generate_spline = function(y.in, dat = GZ.data, bs.in){
  fd.params = y.in
  fd.bs = fd(
    coef = fd.params,
    basisobj = bs.in
  )
  spline.tmp = predict(fd.bs, 1:nrow(dat))
  return(spline.tmp)
}

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


# inputs for bivariate precipiration spline 
rain = GZ.data$avgRainfall
rain[which(GZ.data$avgRainfall == 70)] = NA
rain[which(GZ.data$avgRainfall == 319.667285499247)] = NA

precip.range = seq(0, max(round(rain), na.rm = T), by = 0.1)
precip.range.round = signif(precip.range, 4)
l = 1:49

align_index_precip = function(dat = GZ.data, browse = F){
  if(browse) browser()
  precip.ind = matrix(NA, nrow = nrow(dat), ncol = length(l))
  
  for(ii in (length(l) + 1): nrow(dat)){
    p = dat$avgRainfall[ii - l]
    tmp = unlist(sapply(1:length(l), function(kk) which(precip.range.round == p[kk])))
    if(length(tmp) == length(l)){
      precip.ind[ii, ] = tmp
    }else{
      precip.ind[ii, ] = c(tmp, rep(NA, (length(l) - length(tmp))))
    }
  }
  return(precip.ind)
}
precip.ind = align_index_precip()


# =============================================================== # 
###### set up functions

# create daily spline for mosquito data
# input: y.in and dat
# y.in takes mean MOI
# output: mosq.spline w/ set.order number of knots
generate_spline = function(y.in, dat = GZ.data, bs.in){
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

# write 95% CI
# input: param_posterior=matrix of posterior dist of param of interest
# output: matrix w/ 3 rows - 1st row=upper limit; 2nd row=lower limit; 3rd row=median
write_CI = function(param_posterior, lwr = 0.025, upr = 0.975){
  
  mat = matrix(0, ncol = ncol(param_posterior), nrow = 3)
  
  for(ii in 1:ncol(param_posterior)){
    mat[1, ii] = quantile(param_posterior[, ii], c(lwr, upr), na.rm = T)[1]
    mat[2, ii] = quantile(param_posterior[, ii], c(lwr, upr), na.rm = T)[2]  
    mat[3, ii] = median(param_posterior[, ii], na.rm = T)
  }
  
  return(mat)
  
}

# generate an array of coefficients to input into beta
# each row represents a day; each column represents the temp effect on a given day,
# and each matrix is created with one iteration in the chain
generate_coef = function(m.biv.chain, t.biv.chain, p.biv.chain, dat = GZ.data, browse = F){
  if(browse) browser()
  
  # generate array
  mosq.coefs = array(NA, dim = c(nrow(dat), length(l), nrow(m.biv.chain)))
  temp.coefs = array(NA, dim = c(nrow(dat), length(l), nrow(t.biv.chain)))
  precip.coefs = array(NA, dim = c(nrow(dat), length(l), nrow(p.biv.chain)))
  
  for(ll in 1:nrow(m.biv.chain)){
    biv.mosq = generate_biv_spline(coef.prop = m.biv.chain[ll, ], dat.range = mosq.range, lag.dat = l)
    biv.temp = generate_biv_spline(coef.prop = t.biv.chain[ll, ], dat.range = temp.range, lag.dat = l)
    biv.precip = generate_biv_spline(coef.prop = p.biv.chain[ll, ], dat.range = precip.range, lag.dat = l)
    
    for(tt in (length(l) + 1):nrow(dat)){
      
      temp.tmp = unlist(sapply(1:length(l), function(ii) biv.temp[temp.ind[tt, ii], l[ii]]))
      temp.coefs[tt, , ll] = temp.tmp
      
      mosq.tmp = unlist(sapply(1:length(l), function(ii) biv.mosq[mosq.ind[tt, ii], l[ii]]))
      mosq.coefs[tt, , ll] = mosq.tmp
      
      precip.tmp = unlist(sapply(1:length(l), function(ii) biv.precip[temp.ind[tt, ii], l[ii]]))
      precip.coefs[tt, , ll] = precip.tmp
    } 
  }
  
  list.coef = list(mosq.coefs, temp.coefs, precip.coefs)
  names(list.coef) = c('mosq.coefs', 'temp.coefs', 'precip.coefs')
  return(list.coef)
}

# generate a matrix of betas evaluated on every day in the time period
# each row corresponds to a different sample from the posterior
# inputs: chain = weather.chain, spline.in = spline generated from 
# generate_spline command, dat = GZ.data
# outputs: matrix of beta estimates over the time series where each
# row represents another draw from the posterior
write_beta = function(int.spline.mat, biv.arr, dat = GZ.data, browse = F, extra = F){
  if(browse) browser()
  
  if(!extra){
    mosq.arr = biv.arr$mosq.coefs
    temp.arr = biv.arr$temp.coefs
    precip.arr = biv.arr$precip.coefs
    
    mat = matrix(0, nrow = nrow(int.spline.mat), ncol = nrow(dat))
    lagEnd = length(l) + 1
    newdata.in = matrix(0, nrow = length(lagEnd:nrow(dat)), ncol = (1 + length(l) * 3))
    
    for(ll in 1:nrow(int.spline.mat)){
      newdata.in[, 1] = int.spline.mat[ll, lagEnd:nrow(dat)]
      newdata.in[, 2:50] = mosq.arr[lagEnd:nrow(dat), , ll]
      newdata.in[, 51:99] = temp.arr[lagEnd:nrow(dat), , ll]
      newdata.in[, 100:148] = precip.arr[lagEnd:nrow(dat), , ll]
      tmp = sapply(1:nrow(newdata.in), function(ii) sum(newdata.in[ii, ], na.rm = T))
      mat[ll, ] = c(rep(0, lagEnd - 1), exp(tmp))
    }
    return(mat)
  }else{
    
    mosq.arr = biv.arr$mosq.coefs
    temp.arr = biv.arr$temp.coefs
    precip.arr = biv.arr$precip.coefs
    
    mat = matrix(0, nrow = nrow(int.spline.mat), ncol = nrow(dat))
    lagEnd = length(l) + 1
    newdata.in_keep = array(0, c(length(lagEnd:nrow(dat)), (1+length(l)*3), nrow(int.spline.mat)/10))
    newdata.in_dontKeep = matrix(0, nrow = length(lagEnd:nrow(dat)), ncol = (1 + length(l) * 3))
    
    for(ll in 1:nrow(int.spline.mat)){
      print(ll)
      if(ll %% 10 == 0){
        newdata.in_keep[, 1, ll/10] = int.spline.mat[ll, lagEnd:nrow(dat)]
        newdata.in_keep[, 2:50, ll/10] = mosq.arr[lagEnd:nrow(dat), , ll]
        newdata.in_keep[, 51:99, ll/10] = temp.arr[lagEnd:nrow(dat), , ll]
        newdata.in_keep[, 100:148, ll/10] = precip.arr[lagEnd:nrow(dat), , ll]
        tmp = sapply(1:nrow(newdata.in_keep), function(ii) sum(newdata.in_keep[ii,, ll/10], na.rm = T))
        mat[ll, ] = c(rep(0, lagEnd - 1), exp(tmp))
      }else{
        newdata.in_dontKeep[, 1] = int.spline.mat[ll, lagEnd:nrow(dat)]
        newdata.in_dontKeep[, 2:50] = mosq.arr[lagEnd:nrow(dat), , ll]
        newdata.in_dontKeep[, 51:99] = temp.arr[lagEnd:nrow(dat), , ll]
        newdata.in_dontKeep[, 100:148] = precip.arr[lagEnd:nrow(dat), , ll]
        tmp = sapply(1:nrow(newdata.in_dontKeep), function(ii) sum(newdata.in_dontKeep[ii, ], na.rm = T))
        mat[ll, ] = c(rep(0, lagEnd - 1), exp(tmp))
      }
    }
    return(list(mat, newdata.in_keep))
  }
}

# simulate epidemic
# input: GZ.data = GZ.data, alpha = alpha.value for specific mcmc chain, 
# n.reps = # of simulations, beta_mat = matrix of betas produced from write_beta
# output: matrix with nrows corresponding to n.reps and ncols as time series
sim_epidemic = function(dat = GZ.data, alpha = 1, n.reps, beta.vec, browse = F){
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
      I.new[ind.sim, tt] = rnbinom(length(ind.sim), size = I.old[ind.sim], mu = (I.old[ind.sim] ^ alpha) * beta.vec[tt])
    }
    # print(tt)
  }
  return(I.new)
}

# one-step ahead simulation
# input: GZ.data = GZ.data, alpha = alpha.value for specific mcmc chain, 
# n.reps = # of simulations, beta_mat = matrix of betas produced from write_beta
# output: matrix with nrows corresponding to n.reps and ncols as time series
sim_oneStep = function(dat = GZ.data, alpha = 1, n.reps, beta.vec, browse = F){
  if(browse) browser()
  
  I.new = matrix(0, nrow = n.reps, ncol = nrow(dat))
  I.old = matrix(0, nrow = n.reps, ncol = 1)
  GZ.data_I.old = matrix(GZ.data$I.old, nrow = 1, ncol = nrow(GZ.data))
  imp.cases = matrix(rep(GZ.data$imported.cases, n.reps), nrow = n.reps, ncol = nrow(dat), byrow = T)
  SI = matrix(weights.serialinterval)
  
  for(tt in (length(weights.serialinterval) + 1):nrow(dat)){
    A = GZ.data_I.old[, seq(tt-1, tt-length(weights.serialinterval), -1)]
    A.SI = A %*% SI
    
    Imp = imp.cases[, seq(tt-1, tt-length(weights.serialinterval), -1)]
    Imp.SI = Imp %*% SI
    
    I.old = A.SI + Imp.SI
    
    ind.sim = which(I.old > 0)
    if(length(ind.sim) != 0){
      I.new[ind.sim, tt] = rnbinom(length(ind.sim), size = I.old[ind.sim], mu = (I.old[ind.sim] ^ alpha) * 
                                     min(20, beta.vec[tt]))
    }
    # print(tt)
  }
  return(I.new)
}


# aggregate GZ.data local cases weekly
# inputs
# outputs
casesByWeek = function(dat = GZ.data){
  
  # aggregating weekly
  dat.week = sapply(1:nrow(dat), function(ww) week(GZ.data$date[ww]))
  GZ.dat.tmp = cbind(GZ.data, dat.week)
  cases.weekly = aggregate(local.cases ~ year + dat.week, data = GZ.dat.tmp, sum)
  
  # order by year and week
  cases.weekly = cases.weekly[order(cases.weekly$year, cases.weekly$dat.week), ]
  
  return(cases.weekly$local.cases)
}

# aggregate simulations weekly
# inputs
# outputs
simByWeek = function(sim = I.new.sim, dat = GZ.data){
  
  dat.week = sapply(1:nrow(dat), function(ww) week(GZ.data$date[ww]))
  GZ.dat.tmp = cbind(GZ.data, dat.week)
  
  sim.agg = t(sim)
  sim.aggWeekly = list()
  
  for(ii in 1:ncol(sim.agg)){
    sim.aggWeekly[[ii]] = 
      aggregate(sim.agg[, ii], by = list(GZ.dat.tmp$year, GZ.dat.tmp$dat.week), FUN = sum)
  }
  
  sim.weekly = list()
  for(ii in 1:length(sim.aggWeekly)){
    sim.weekly[[ii]] = sim.aggWeekly[[ii]][order(sim.aggWeekly[[ii]]$Group.1, 
                                                 sim.aggWeekly[[ii]]$Group.2), ]
  }
  
  mat = matrix(0, nrow = length(sim.weekly), ncol = nrow(sim.weekly[[1]]))
  for(ii in 1:length(sim.weekly)){
    mat[ii, ] = sim.weekly[[ii]]$x
  }
  
  return(mat)
  
}

# rank simulations yearly
# inputs 
# outputs
simRank = function(I.new.sim = I.new.sim){
  
  # subset simulations by year
  sim.2005 = I.new.sim[, which(GZ.data$year == 2005)]
  sim.2006 = I.new.sim[, which(GZ.data$year == 2006)]
  sim.2007 = I.new.sim[, which(GZ.data$year == 2007)]
  sim.2008 = I.new.sim[, which(GZ.data$year == 2008)]
  sim.2009 = I.new.sim[, which(GZ.data$year == 2009)]
  sim.2010 = I.new.sim[, which(GZ.data$year == 2010)]
  sim.2011 = I.new.sim[, which(GZ.data$year == 2011)]
  sim.2012 = I.new.sim[, which(GZ.data$year == 2012)]
  sim.2013 = I.new.sim[, which(GZ.data$year == 2013)]
  sim.2014 = I.new.sim[, which(GZ.data$year == 2014)]
  sim.2015 = I.new.sim[, which(GZ.data$year == 2015)]
  
  sim.order = matrix(nrow = nrow(I.new.sim), ncol = 11)
  
  for(ii in 1:nrow(I.new.sim)){
    sim.order[ii, ] = 
      order(c(sum(sim.2005[ii, ]), sum(sim.2006[ii, ]), sum(sim.2007[ii, ]), sum(sim.2008[ii, ]),
              sum(sim.2009[ii, ]), sum(sim.2010[ii, ]), sum(sim.2011[ii, ]), sum(sim.2012[ii, ]),
              sum(sim.2013[ii, ]), sum(sim.2014[ii, ]), sum(sim.2015[ii, ])), decreasing = TRUE)
  }
  
  mat.rank = matrix(0, nrow = 11, ncol = 11)
  
  # gives year as rows and rank as columns
  for(ii in 1:11){
    for(jj in 1:11){
      mat.rank[jj, ii] = length(which(sim.order[, ii] == jj))
    }
  }
  return(mat.rank)  
}

# =============================================================== # 
###### calculate values for figures
setDate = 'rainfall'

# creating int spline matrix
int.mat = matrix(0, ncol = nrow(GZ.data), nrow = nrow(int.chain))
for(ii in 1:nrow(int.mat)){
  int.mat[ii, ] = generate_spline(as.numeric(int.chain[ii, ]), bs.in = b0_basis)
}

# calculate the coefficient inputs for beta
f = paste0('../output/BIVcoef_', setDate, '.RData')
if(!file.exists(f)){
  biv.coef.arr = generate_coef(m.biv.chain = mosq.biv.chain, t.biv.chain = weather.biv.chain, 
                               p.biv.chain = precip.biv.chain, browse = F)
  save(biv.coef.arr, file = f)
}else{
  load(f)
}

# calculate beta and mosquito matrix
# generate beta
f = paste0('../output/betaMat_', setDate, '.RData')
if(!file.exists(f)){
  beta_mat = write_beta(int.spline.mat = int.mat, biv.arr = biv.coef.arr, browse = F)
  save(beta_mat, file = f)
}else{
  load(f)
}

beta_CI = write_CI(beta_mat)


# creating simulation
f = paste0('../output/sim_', setDate, '.RData')
# f = '../output/sim_V1_1_V1_1_20170918_1.RData'
if(!file.exists(f)){
  reps = 2e3
  I.new.sim = matrix(NA, nrow = reps, ncol = nrow(GZ.data))
  for(ii in 1:reps){
    pp = sample(1:nrow(beta_mat), 1)
    beta.samp = beta_mat[pp, ]
    I.new.sim[ii, ] = sim_epidemic(n.reps = 1, beta.vec = beta.samp, browse = F)
  }
  # correlation between median simulation and data
  median.sim = write_CI(param_posterior = I.new.sim)[3,]
  median.corr = cor(GZ.data$local.cases, median.sim)
  save(I.new.sim, median.corr, file = f)
}else{
  load(f)
}

# =============================================================== # 
###### Figure 5. Dengue cases aggregated weekly for 2005-2015

f = paste0('../output/fig_modelVal_rainfall.pdf')
pdf(f, height = 4.5, width = 6)
layout(matrix(1:2))
par(oma = c(2.5,2,1,1), mar = c(0.5, 1.5, 0,0))

cases_agg = log10(casesByWeek())
cases_agg[is.infinite(cases_agg)] = 0

# regular simulation
sim_agg = log10(simByWeek(sim = I.new.sim))
sim_CI = write_CI(param_posterior = sim_agg)
sim_CI[is.infinite(sim_CI)] = 0
time = 1:ncol(sim_agg)
plot(-1000, -1000, xlim = c(0, ncol(sim_agg)), ylim = c(0, 5), xaxs = 'i', xlab = '', ylab = '',
     xaxt = 'n', las = 1)
polygon(c(time, rev(time)), c(sim_CI[1, ], rev(sim_CI[2, ])), density = NA, 
        col = rgb(0,0,1,0.3), border = NA)
lines(cases_agg, lwd = 2, col = rgb(0,0,0,0.75))
# cor(write_CI(I.new.sim)[3,], GZ.data$local.cases)
text(x = length(cases_agg) * 0.92, y = 4.8, labels = expression(rho* ' = 0.271'))
# paste0('corr = ', round(median.corr, 3))

text(x=15, y = 5 * 0.92, labels = 'A', cex = 1.5)

plot(-1000, -1000, xlim = c(0, ncol(sim_agg)), ylim = c(0, 5), xaxs = 'i', xlab = '', ylab = '',
     xaxt = 'n', las = 1)
polygon(c(time, rev(time)), c(sim_CI[1, ], rev(sim_CI[2, ])), density = NA, 
        col = rgb(0,0,1,0.3), border = NA)
set.seed(1)
ss = sample(1:2000, 4)
samples = sim_agg[ss,]
samples[is.infinite(samples)] = 0
lines(samples[1,], lwd = 2, col = rgb(1,0,0,0.75))
lines(samples[2,], lwd = 2, col = rgb(1,0.5,0,0.6))
lines(samples[3,], lwd = 2, col = rgb(0,1,1,0.5))
lines(samples[4,], lwd = 2, col = rgb(0,1,0,0.6))
text(x=15, y = 5 * 0.92, labels = 'B', cex = 1.5)

# axis
tick = seq(from = 0, to = 583, by = 54)
lbl = 2005:2015
axis(side = 1, at = tick, labels = FALSE)
axis(side = 1, at = tick + 25, tick = FALSE, labels = lbl)
mtext(expression(log[10]*'(incidence)'), side=2, line=0.5, outer = T)
mtext('Year', side=1, line=2.1)

dev.off()

# =============================================================== # 
###### Figure 6. Ranking of years by total annual dengue incidence simulated by the fitted model
###### given weather, mosquito, and imported case data from that year 

f = paste0('../output/fig_modelVal_red-blue_rainfall.pdf')
pdf(f, height = 6.5, width = 8)
layout(matrix(c(1,2,2,2,2), 5, 1, byrow = T))
par(oma = c(4,3.5,1.5,1), mar = c(0,2,0,2))

# first pane
# incidence
yearly.sum = sapply(2005:2015, function(yy)sum(GZ.data$local.cases[GZ.data$year==yy]))
# yearly.sum[1] = 1
yearly.sum = as.numeric(as.character(yearly.sum))
bb = barplot(log10(yearly.sum), space = rep(0, 11), xaxs = 'i', ylim = c(0, 1.4*max(log10(yearly.sum))),
             col = rgb(0,0,0,0.5), yaxs = 'i', cex.axis = 1.2, las = 1)
# bb = barplot(yearly.sum, space = rep(0, 11), xaxs = 'i', ylim = c(0, 1.17*max(yearly.sum)),
#              col = rgb(0,0,0,0.5), yaxs = 'i', cex.axis = 1.2, las = 1)
text(x =bb, y = log10(yearly.sum), labels = yearly.sum, pos = 3, cex =1.3)
text(x = 0.5, y = 0, labels = '0', pos = 3, cex = 1.2)
mtext(expression(log[10]*'(incidence)'), side=2, line=3.3, cex = 1.2)
box()

# second pane
# barplot
mat.rank = simRank(I.new.sim = I.new.sim)
barplot(t(mat.rank), col=brewer.pal(11, 'RdBu'), space = rep(0, 11), xaxs = 'i', yaxt = 'n')
tick = c(0, nrow(I.new.sim) * 0.25 , nrow(I.new.sim) * .5, nrow(I.new.sim) * .75)
lbl = c(0, 0.25, 0.5, 0.75)
axis(side = 2, at = tick, labels = lbl, cex.axis = 1.2, las = 1)
axis(side = 1, at = (1:11) - 0.5, labels = 2005:2015, tick = F, lty = 0, line = -0.5)

mtext('Year', side=1, line=2.7, cex = 1.2)
mtext('Posterior probability of a given ranking', side=2, line=3.5, cex = 1.2)

dev.off()


# =============================================================== # 
###### Figure 4. Contribution to beta plot

f = '../output/contribution_to_beta_rainfall.RData'
if(!file.exists(f)){
  beta_extra = write_beta(int.spline.mat = int.mat, biv.arr = biv.coef.arr, 
                          extra = T, browse = F)
  save(beta_extra, file = f)
}else{
  load(f)
}
beta_b0 = beta_extra[[2]][,1,]
beta_mosq = beta_extra[[2]][,2:50,]
beta_temp = beta_extra[[2]][,51:99,]
beta_precip = beta_extra[[2]][,100:148,]

f = paste0('../output/contribution_to_beta_rainfall.pdf')
pdf(f, height = 6, width = 8)
layout(matrix(1:5),5,1)
# layout.show()
par(oma = c(3.1,4,1.3,1), mar = c(0,1,0,0))

col.vec = rainbow(5)

# beta
time = 1:nrow(GZ.data)
plot(-1000, -1000, xlim = c(0, nrow(GZ.data)), ylim = c(0, 20), xaxt = 'n',
     xlab = '', ylab = '', xaxs = 'i', las = 1)
polygon(c(time, rev(time)), c(beta_CI[1,], rev(beta_CI[2,])), density = NA,
        col = rgb(0,0,1,0.3), border = NA)
ii.lines = seq(0, 50, by = 10)
ii.col = 1
for(ii in ii.lines){
  lines(beta_mat[ii, ], col = adjustcolor(col = col.vec[ii.col], alpha.f = 0.8), lwd = 1.5)
  ii.col = ii.col + 1
}
lines(beta_CI[3,], lwd = 2, col = rgb(0,0,0,0.7))
abline(h = 1, col = rgb(1,0,0,0.75))
abline(v = seq(from = 1, to = 4017, by = 366), lwd = 1, col = rgb(0,0,0,0.3))
mtext(text = 'Transmission', side = 2, line = 3.4, cex = 0.7)
mtext(text = expression('coefficient ('*beta*')'), side = 2, line = 2, cex = 0.7)
text(x=65, y = 20*0.9, labels = 'A', cex = 2)

# b0
b0_CI = write_CI(int.mat)
plot(-100,-100, xlab = '', ylab = '', xlim = c(0,4017), 
     # ylim = c(min(exp(beta_b0)), max(exp(beta_b0))), 
     ylim = c(min(beta_b0), max(beta_b0)),
     xaxt = 'n', las = 1,
     xaxs = 'i')
polygon(c(time, rev(time)), c(b0_CI[1,], rev(b0_CI[2,])), density = NA, 
        col = rgb(0,0,1,0.3),
        border = NA)
for(ii in 1:5){
  # lines(c(rep(0,49),exp(beta_b0[, ii])), col = rgb(0,0,0,0.2))
  lines(c(rep(0,49),beta_b0[, ii]), col = adjustcolor(col = col.vec[ii], alpha.f = 0.8), lwd = 1.5)
}
lines(b0_CI[3,], lwd = 2, col = rgb(0,0,0,0.7))
abline(v = seq(from = 1, to = 4017, by = 366), lwd = 1, col = rgb(0,0,0,0.3))
text(x=65, y = max(beta_b0)*0.9, labels = 'B', cex = 2)
mtext(text = expression(''*beta[0]*'(t)'), side = 2, line = 3.2, cex = 0.7)

# calculate the lines for mosq and temp and precip contributions
mosq.cont = matrix(NA, nrow = 100, ncol = nrow(GZ.data))
for(ii in 1:dim(beta_mosq)[3]){
  tmp = sapply(1:nrow(beta_mosq), function(kk) sum(beta_mosq[kk,,ii]))
  mosq.cont[ii, ] = c(rep(0, 49), (tmp - mean(tmp)))
}
mosq.cont_CI = write_CI(mosq.cont)

temp.cont = matrix(NA, nrow = 100, ncol = nrow(GZ.data))
for(ii in 1:dim(beta_temp)[3]){
  tmp = sapply(1:nrow(beta_temp), function(kk) sum(beta_temp[kk,,ii]))
  temp.cont[ii, ] = c(rep(0, 49), (tmp - mean(tmp, na.rm = T)))
}
temp.cont_CI = write_CI(temp.cont)

precip.cont = matrix(NA, nrow = 100, ncol = nrow(GZ.data))
for(ii in 1:dim(beta_precip)[3]){
  tmp = sapply(1:nrow(beta_precip), function(kk) sum(beta_precip[kk,,ii]))
  precip.cont[ii, ] = c(rep(0, 49), (tmp - mean(tmp, na.rm = T)))
}
precip.cont_CI = write_CI(precip.cont)


# mosq
plot(-100,-100, xlab = '', ylab = '', xlim = c(0,4017), ylim = c(min(mosq.cont, na.rm = T), 
                                                                 max(mosq.cont, na.rm = T)), 
     xaxt = 'n', xaxs = 'i', las = 1)
polygon(c(time, rev(time)), c(mosq.cont_CI[1,], rev(mosq.cont_CI[2,])), density = NA, 
        col = rgb(0,0,1,0.3),
        border = NA)
for(ii in 1:5){
  lines(mosq.cont[ii,],col = adjustcolor(col = col.vec[ii], alpha.f = 0.8), lwd = 1.5)
}
lines(mosq.cont_CI[3,], lwd = 2, col = rgb(0,0,0,0.7))
abline(v = seq(from = 1, to = 4017, by = 366), lwd = 1, col = rgb(0,0,0,0.3))
text(x=65, y = max(mosq.cont, na.rm = T) * 0.87, labels = 'C', cex = 2)
mtext(text = expression(sum(s[m](m(t-tau), tau), tau = 1, 49)), side =2 , line = 1.95, cex = 0.7)

# temp
ind = 1:3987
plot(-100,-100, xlab = '', ylab = '', xlim = c(0,4017), ylim = c(min(temp.cont, na.rm = T), 
                                                                 max(temp.cont, na.rm = T)), 
     xaxt = 'n', xaxs = 'i', las = 1)
polygon(c(time[ind], rev(time[ind])), 
        c(temp.cont_CI[1, ind], rev(temp.cont_CI[2,ind])), density = NA, 
        col = rgb(0,0,1,0.3),
        border = NA)
for(ii in 1:5){
  lines(temp.cont[ii,], col = adjustcolor(col = col.vec[ii], alpha.f = 0.8), lwd = 1.5)
}
lines(temp.cont_CI[3,], lwd = 2, col = rgb(0,0,0,0.7))
abline(v = seq(from = 1, to = 4017, by = 366), lwd = 1, col = rgb(0,0,0,0.3))
mtext(text = expression(sum(s[T](T(t-tau), tau), tau = 1, 49)), side =2 , line = 1.95, cex = 0.7)
text(x=65, y = max(temp.cont, na.rm = T) * 0.85, labels = 'D', cex = 2)

# temp
ind = 1:3987
plot(-100,-100, xlab = '', ylab = '', xlim = c(0,4017), ylim = c(min(precip.cont, na.rm = T), 
                                                                 max(precip.cont, na.rm = T)), 
     xaxt = 'n', xaxs = 'i', las = 1)
polygon(c(time[ind], rev(time[ind])), 
        c(precip.cont_CI[1, ind], rev(precip.cont_CI[2,ind])), density = NA, 
        col = rgb(0,0,1,0.3),
        border = NA)
for(ii in 1:5){
  lines(precip.cont[ii,], col = adjustcolor(col = col.vec[ii], alpha.f = 0.8), lwd = 1.5)
}
lines(precip.cont_CI[3,], lwd = 2, col = rgb(0,0,0,0.7))
abline(v = seq(from = 1, to = 4017, by = 366), lwd = 1, col = rgb(0,0,0,0.3))
mtext(text = expression(sum(s[P](P(t-tau), tau), tau = 1, 49)), side =2 , line = 1.95, cex = 0.7)
text(x=65, y = max(precip.cont, na.rm = T) * 0.87, labels = 'E', cex = 2)



# adding axes
tick = seq(from = 1, to = 4017, by = 366)
axis(side = 1, at = tick, labels = FALSE)
tmp = rep(0, 11)
names(tmp) = 2005:2015
par(new=TRUE)
barplot(tmp, ylim = c(0,2000), axes = F, yaxs = 'i', xaxs = 'i', space = rep(0,11))
mtext(side = 1, outer = T, text = "Year", line = 2.05, cex = 0.7)

dev.off()

# =============================================================== # 
## supplementary figures


require(lubridate)
# =============================================================== # 
###### set up functions

# aggregate GZ.data local cases weekly
# inputs
# outputs
casesByWeek = function(dat = GZ.data){
  
  # aggregating weekly
  dat.week = sapply(1:nrow(dat), function(ww) week(GZ.data$date[ww]))
  GZ.dat.tmp = cbind(GZ.data, dat.week)
  cases.weekly = aggregate(local.cases ~ year + dat.week, data = GZ.dat.tmp, sum)
  
  # order by year and week
  cases.weekly = cases.weekly[order(cases.weekly$year, cases.weekly$dat.week), ]
  
  return(cases.weekly$local.cases)
}

# aggregate simulations weekly
# inputs
# outputs
simByWeek = function(sim = I.new.sim, dat = GZ.data){
  
  dat.week = sapply(1:nrow(dat), function(ww) week(GZ.data$date[ww]))
  GZ.dat.tmp = cbind(GZ.data, dat.week)
  
  sim.agg = t(sim)
  sim.aggWeekly = list()
  
  for(ii in 1:ncol(sim.agg)){
    sim.aggWeekly[[ii]] = 
      aggregate(sim.agg[, ii], by = list(GZ.dat.tmp$year, GZ.dat.tmp$dat.week), FUN = sum)
  }
  
  sim.weekly = list()
  for(ii in 1:length(sim.aggWeekly)){
    sim.weekly[[ii]] = sim.aggWeekly[[ii]][order(sim.aggWeekly[[ii]]$Group.1, 
                                                 sim.aggWeekly[[ii]]$Group.2), ]
  }
  
  mat = matrix(0, nrow = length(sim.weekly), ncol = nrow(sim.weekly[[1]]))
  for(ii in 1:length(sim.weekly)){
    mat[ii, ] = sim.weekly[[ii]]$x
  }
  
  return(mat)
  
}

# aggregate case data weekly
date_week = sapply(1:nrow(GZ.data), function(ww) week(GZ.data$date[ww]))
GZ.data = cbind(GZ.data, date_week)
casesWeekly = aggregate(local.cases ~ year + date_week, data = GZ.data, sum)
# order by year and week
cases_agg = casesWeekly[order(casesWeekly$year, casesWeekly$date_week), ]

sim_agg = simByWeek(sim = I.new.sim)

# subset simulations by year
sim_agg_yearly = function(y){
  return(sim = sim_agg[, which(cases_agg$year == y)])
}
years = 2005:2015
for(yy in years){
  eval(parse(text = paste0(
    'sim.', yy, ' = sim_agg_yearly(y = ', yy, ')'
  )))
}

# write 95% CI
# input: param_posterior=matrix of posterior dist of param of interest
# output: matrix w/ 3 rows - 1st row=upper limit; 2nd row=lower limit; 3rd row=median
write_CI = function(param_posterior, lwr = 0.025, upr = 0.975){
  
  mat = matrix(0, ncol = ncol(param_posterior), nrow = 3)
  
  for(ii in 1:ncol(param_posterior)){
    mat[1, ii] = quantile(param_posterior[, ii], c(lwr, upr), na.rm = T)[1]
    mat[2, ii] = quantile(param_posterior[, ii], c(lwr, upr), na.rm = T)[2]  
    mat[3, ii] = median(param_posterior[, ii], na.rm = T)
  }
  
  return(mat)
  
}

# =============================================================== # 
###### density plots
####### Figure S1. Posterior predictive distribution of total dengue incidence for
####### each year based on 2,000 simulations from the fitted transmission model. 

plot_S1 = function(sim, y, browse = F){
  if(browse) browser()
  
  dens = log10(rowSums(sim))
  # if(sum(is.infinite(dens))!=0){
  #   dens = dens[-which(is.infinite(dens))]
  # }
  # if(sum(is.na(dens))!=0){
  #   dens[which(is.na(dens))] = 100
  # }
  
  # print(range(dens))
  
  if(y == 2005){
    v = 0
  }else{
    v = log10(sum(cases_agg$local.cases[cases_agg$year == y]))
  }
  quant = sum(dens <= v)
  txt = quant / length(dens)
  
  # print(txt)
  print(y)
  print(paste0('median = ', median(dens, na.rm = T)))
  print(paste0('range = ', range(dens, na.rm = T)))
  print(txt)
}

for(yy in years){
  eval(parse(text = paste0(
    'plot_S1(sim = sim.', yy, ', y = ', yy, ')'
  )))
}


####### Figure S2. Posterior predictive distribution of peak weekly incidence for each year 
####### based on 2,000 simulations from the fitted transmission model.  
max_S2 = function(sim){
  return(sapply(1:nrow(sim), function(ii) max(sim[ii, ])))
}
for(yy in years){
  eval(parse(text = paste0(
    'max.', yy, ' = max_S2(sim.', yy, ')'
  )))
}

plot_S2 = function(max.in, y){
  dens = log10(max.in)
  
  # if(sum(is.infinite(dens))!=0){
  #   dens = dens[-which(is.infinite(dens))]
  # }
  # if(sum(is.na(dens))!=0){
  #   dens[which(is.na(dens))] = 100
  # }
  
  if(y == 2005){
    v = 0
  }else{
    v = log10(max(cases_agg$local.cases[cases_agg$year == y]))
  }
  quant = sum(dens <= v)
  txt = quant / length(dens)

  print(y)
  print(paste0('median = ', median(dens)))
  print(paste0('range = ', range(dens)))
  print(txt)
}

for(yy in years){
  eval(parse(text = paste0(
    'plot_S2(max.in = max.', yy, ', y = ', yy, ')'
  )))
}

####### Figure S3. Posterior predictive distribution of the number of weeks with non-zero
####### dengue incidence in each year based on 2,000 simulations from the fitted transmission model.
out_S3 = function(sim){
  return(sapply(1:nrow(sim), function(ii) length(which(sim[ii, ] != 0))))
}
for(yy in years){
  eval(parse(text = paste0(
    'out.', yy, ' = out_S3(sim.', yy, ')'
  )))
}

plot_S3 = function(out, y){
  dens = out
  v = length(which(cases_agg$local.cases[cases_agg$year == y] != 0))
  quant = sum(dens <= v)
  txt = quant / length(dens)
  
  # print(txt)
  print(y)
  print(paste0('median = ', median(dens)))
  print(paste0('range = ', range(dens)))
}

for(yy in years){
  eval(parse(text = paste0(
    'plot_S3(out = out.', yy, ', y = ', yy, ')'
  )))
}

####### Figure S4. Posterior predictive distribution of length of dengue seasons
####### in each year based on 2,000 simulations from the fitted transmission model.

diff_S4 = function(sim, browse = F){
  if(browse) browser()
  seas = sapply(1:nrow(sim), function(ii) which(sim[ii, ] != 0))
  seas.0 = sapply(1:length(seas), function(ii) which(length(seas[[ii]]) == 0))
  ind = which(seas.0 == 1)
  if(length(ind) != 0){
    for(ii in ind){
      seas[[ii]] = 0
    }
  }
  diff = sapply(1:length(seas), function(ii) seas[[ii]][length(seas[[ii]])] - seas[[ii]][1])
  return(unlist(diff))
}
for(yy in years){
  eval(parse(text = paste0(
    'diff.', yy, ' = diff_S4(sim.', yy, ')'
  )))
}

plot_S4 = function(dif, y){
  dens = dif
  seas = which(cases_agg$local.cases[cases_agg$year == y] != 0)
  if(yy == 2005){
    v = 0
  }else{
    v =  seas[length(seas)] - seas[1]
  }
  quant = sum(dens <= v)
  txt = quant / length(dens)
  
  # print(txt)
  print(y)
  print(paste0('median = ', median(dens)))
  print(paste0('range = ', range(dens)))
}

for(yy in years){
  eval(parse(text = paste0(
    'plot_S4(dif = diff.', yy, ', y = ', yy, ')'
  )))
}
