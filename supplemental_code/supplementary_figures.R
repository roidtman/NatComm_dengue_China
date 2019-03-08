###### Supplementary figures ######

# =============================================================== # 
###### load & set up data
mosq = read.csv('../data/guangzhou_vector_05-15.csv')
GZ.data = read.csv('../data/guangzhou_case_weather_05-15.csv')
# load('../output/sim.RData')
# load simulation output that is generated in main_text_figures.R

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

plot_S1 = function(sim, y){
  dens = log10(rowSums(sim))
  if(y == 2005){
    v = 0
  }else{
    v = log10(sum(cases_agg$local.cases[cases_agg$year == y]))
  }
  quant = sum(dens <= v)
  txt = quant / length(dens)
  
  plot(density(dens), col = rgb(0,0,1,0.2),
       xlim = c(0,6), xlab = '', ylab = '', main = '', las = 1)
  polygon(density(dens), col = rgb(0,0,1,0.3), border = NA)
  abline(v = v, col = 'red', lwd = 1.5)
  mtext(y, side = 3, cex = 0.8, line = 0.5)
  
  z = max(density(dens)$y)
  text(x = 5.5, y = z - (.1*z), labels = signif(txt, 3))
}
pdf('../output/figS1.pdf', height = 6, width = 8)
layout(matrix(1:12, 3,4, byrow = T))
par(oma=c(2,2,2,2), mar=rep(2, 4))

for(yy in years){
  eval(parse(text = paste0(
    'plot_S1(sim = sim.', yy, ', y = ', yy, ')'
  )))
}
mtext('Posterior probability', outer = T, side = 2, line = 0.5)
mtext(expression(log[10] * '(incidence)'), side = 1, line = 0.5, outer = T)
mtext('Annual local incidence', side = 3, line = 0.5, outer = T)

dev.off()

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
  if(y == 2005){
    v = 0
  }else{
    v = log10(max(cases_agg$local.cases[cases_agg$year == y]))
  }
  quant = sum(dens <= v)
  txt = quant / length(dens)
  
  plot(density(dens), col = rgb(0,0,1,0.3), xlim = c(0, 6), xlab = '', 
       ylab = '', main = '', las = 1)
  polygon(density(dens), col = rgb(0,0,1,0.3), border = NA)
  abline(v = v, col = 'red', lwd = 1.5)
  mtext(y, side = 3, cex = 0.8, line = 0.5)
  
  z = max(density(dens)$y)
  text(x = 5, y = z - (.1*z), labels = signif(txt, 3))
}

pdf('../output/figS2.pdf', height = 6, width = 8)
layout(matrix(1:12, 3,4, byrow = T))
par(oma=c(2,2,2,2), mar=rep(2, 4))

for(yy in years){
  eval(parse(text = paste0(
    'plot_S2(max.in = max.', yy, ', y = ', yy, ')'
  )))
}
mtext('Posterior probability', outer = T, side = 2, line = 0.5)
mtext(expression(log[10] * '(incidence)'), side = 1, line = 0.5, outer = T)
mtext('Peak weekly incidence', side = 3, line = 0.5, outer = T)

dev.off()

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
  
  plot(density(dens), col = rgb(0,0,1,0.3), xlim = c(0, 53), main = '', xlab = '', 
       ylab = '', las = 1)
  mtext(y, side = 3, cex = 0.8, line = 0.5)
  polygon(density(dens), col = rgb(0,0,1,0.3), border = NA)
  abline(v = v, col = 'red', lwd = 1.5)
  
  z = max(density(dens)$y)
  text(x = 45, y = z - (.1*z), labels = signif(txt, 3))
}
pdf('../output/figS3.pdf', heigh = 6, width = 8)
layout(matrix(1:12, 3,4, byrow = T))
par(oma=c(2,2,2,2), mar=rep(2, 4))
for(yy in years){
  eval(parse(text = paste0(
    'plot_S3(out = out.', yy, ', y = ', yy, ')'
  )))
}
mtext('Posterior probability', outer = T, side = 2, line = 0.75)
mtext('Number of weeks', side = 1, line = 0.5, outer = T)
mtext('Number of weeks with nonzero transmission', side = 3, line = 0.5, outer = T)

dev.off()

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
  
  plot(density(dens), col = rgb(0,0,1,0.3), xlim = c(0,53),
       main = '', xlab = '', ylab = '', las = 1)
  polygon(density(dens), col = rgb(0,0,1,0.3), border = NA)
  abline(v = v, col = 'red', lwd = 1.5)
  mtext(y, side = 3, cex = 0.8, line = 0.5)
  
  z = max(density(dens)$y)
  text(x = 45, y = z - (.1*z), labels = signif(txt, 3))
}

pdf('../output/figS4.pdf', height = 6, width = 8)
layout(matrix(1:12, 3,4, byrow = T))
par(oma=c(2,2,2,2), mar=rep(2, 4))
for(yy in years){
  eval(parse(text = paste0(
    'plot_S4(dif = diff.', yy, ', y = ', yy, ')'
  )))
}
mtext('Posterior probability', outer = T, side = 2, line = 0.75)
mtext('Number of weeks', side = 1, line = 0.5, outer = T)
mtext('Length of transmission season', side = 3, line = 0.5, outer = T)

dev.off()

####### Figure S7. Density plots
load('../output/posterior_SMC_V3_20171215_1.RData')
param.chain_1 = getSample(out, numSamples = 1000)
mosq.biv.chain_1 = param.chain_1[, 1:9]
weather.biv.chain_1 = param.chain_1[, 10:18]
int.chain_1 = param.chain_1[, 19:ncol(param.chain_1)]

load('../output/posterior_SMC_V3_20171215_2.RData')
param.chain_2 = getSample(out, numSamples = 1000)
mosq.biv.chain_2 = param.chain_2[, 1:9]
weather.biv.chain_2 = param.chain_2[, 10:18]
int.chain_2 = param.chain_2[, 19:ncol(param.chain_2)]

load('../output/posterior_SMC_V3_20171215_3.RData')
param.chain_3 = getSample(out, numSamples = 1000)
mosq.biv.chain_3 = param.chain_3[, 1:9]
weather.biv.chain_3 = param.chain_3[, 10:18]
int.chain_3 = param.chain_3[, 19:ncol(param.chain_3)]

# function 

densityPlots = function(chain1, chain2, chain3){
  name.vec = c('m1', 'm2', 'm3', 'm4', 'm5', 'm6', 'm7', 'm8', 'm9',
               't1', 't2', 't3', 't4', 't5', 't6', 't7', 't8', 't9',
               'b0.1', 'b0.2', 'b0.3', 'b0.4', 'b0.5', 'b0.6', 'b0.7', 
               'b0.8', 'b0.9', 'b0.10', 'b0.11', 'b0.12', 'b0.13', 'b0.14', 
               'b0.15', 'b0.16', 'b0.17', 'b0.18', 'b0.19', 'b0.20', 'b0.21', 
               'b0.22', 'b0.23', 'b0.24', 'b0.25', 'b0.26', 'b0.27', 'b0.28',
               'b0.29', 'b0.30', 'b0.31', 'b0.32', 'b0.33')
  
  for(ii in 1:ncol(chain1)){
    ymax = max(max(density(chain1[,ii])$y), max(density(chain2[,ii])$y),
               max(density(param.chain_3[,ii])$y))
    xmax = max(max(density(chain1[,ii])$x), max(density(chain2[,ii])$y),
               max(density(param.chain_3[,ii])$x))
    xmin = min(min(density(chain1[,ii])$x), min(density(chain2[,ii])$y),
               min(density(chain3[,ii])$x))
    
    plot(density(chain1[,ii]), main = '', col = rgb(0,0,1,0.1), ylim = c(0, ymax), 
         xlim = c(xmin, xmax))
    polygon(density(chain1[,ii]), border = NA, col = rgb(0,0,1,0.2))
    mtext(text = name.vec[ii], side = 3, line = 0.5, cex = 0.75)
    par(new = T)
    plot(density(chain2[,ii]), main = '', col = rgb(0,0,1,0.1), ylim = c(0, ymax),
         xlim = c(xmin, xmax), xaxt = 'n', yaxt = 'n')
    polygon(density(chain2[,ii]), border = NA, col = rgb(0,0,1,0.2))
    par(new = T)
    plot(density(chain3[,ii]), main = '', col = rgb(0,0,1,0.1), ylim = c(0, ymax),
         xlim = c(xmin, xmax), yaxt = 'n', xaxt = 'n')
    polygon(density(chain3[,ii]), border = NA, col = rgb(0,0,1,0.2))
  }
}

pdf('../output/figS5.pdf', height = 9.5, width = 6)
layout(matrix(1:54, nrow = 9, ncol = 6, byrow = T))
par(mar = c(1.5, 1.5, 2.5, 1), oma = c(1,1,1,0))
densityPlots(param.chain_1, param.chain_2, param.chain_3)
dev.off()

####### Table S2. Gelman-Rubin stats 
require(coda)
require(BayesianTools)
for(load.count in 1:3){
  eval(parse(text = paste0(
    'load("../output/posterior_SMC_V3_20171215_', load.count, '.RData")'
  )))
  param.chain = getSample(out, numSamples = 1000)
  
  posterior.mosq = as.mcmc(as.matrix(data.frame(
    m1 = param.chain[,1],
    m2 = param.chain[,2],
    m3 = param.chain[,3],
    m4 = param.chain[,4],
    m5 = param.chain[,5],
    m6 = param.chain[,6],
    m7 = param.chain[,7],
    m8 = param.chain[,8],
    m9 = param.chain[,9],
    t1 = param.chain[,10],
    t2 = param.chain[,11],
    t3 = param.chain[,12],
    t4 = param.chain[,13],
    t5 = param.chain[,14],
    t6 = param.chain[,15],
    t7 = param.chain[,16],
    t8 = param.chain[,17],
    t9 = param.chain[,18],
    b0.1 = param.chain[,19],
    b0.2 = param.chain[,20],
    b0.3 = param.chain[,21],
    b0.4 = param.chain[,22],
    b0.5 = param.chain[,23],
    b0.6 = param.chain[,24],
    b0.7 = param.chain[,25],
    b0.8 = param.chain[,26],
    b0.9 = param.chain[,27],
    b0.10 = param.chain[,28],
    b0.11 = param.chain[,29],
    b0.12 = param.chain[,30],
    b0.13 = param.chain[,31],
    b0.14 = param.chain[,32],
    b0.15 = param.chain[,33],
    b0.16 = param.chain[,34],
    b0.17 = param.chain[,35],
    b0.18 = param.chain[,36],
    b0.19 = param.chain[,37],
    b0.20 = param.chain[,38],
    b0.21 = param.chain[,39],
    b0.22 = param.chain[,40],
    b0.23 = param.chain[,41],
    b0.24 = param.chain[,42],
    b0.25 = param.chain[,43],
    b0.26 = param.chain[,44],
    b0.27 = param.chain[,45],
    b0.28 = param.chain[,46],
    b0.29 = param.chain[,47],
    b0.30 = param.chain[,48],
    b0.31 = param.chain[,49],
    b0.32 = param.chain[,50],
    b0.33 = param.chain[,51]
  )))
  name.post = paste('posterior.mosq_', load.count, sep = '')
  assign(name.post, posterior.mosq)
  
  name.chain = paste('param.chain_', load.count, sep = '')
  assign(name.chain, param.chain)
}

posterior.mosq = mcmc.list(posterior.mosq_1, posterior.mosq_2, posterior.mosq_3)
GR.diag = gelman.diag(posterior.mosq)
print(GR.diag)
gelman.plot(posterior.mosq)

####### Figure S8. Covariance matrix
if(!require(corrplot)){install.packages('corrplot'); library(corrplot)}
load('../output/posterior_SMC_V3_20171215_1.RData')
param.chain = getSample(out, numSamples = 333)
load('../output/posterior_SMC_V3_20171215_2.RData')
param.chain = rbind(param.chain, getSample(out, numSamples = 333))
load('../output/posterior_SMC_V3_20171215_3.RData')
param.chain = rbind(param.chain, getSample(out, numSamples = 334))

pdf('../output/figS6.pdf', height = 7, width = 8)
layout(1)
par(mar = c(5,4,4,4))
corrplot(cor(param.chain), method = 'square', tl.col = 'white', tl.cex = 0.2)
dev.off()

####### Figure S9. Factorial experiment epidemic trajectory
load('../output/interchange.RData')

pdf('../output/figS5.pdf', height = 8, width = 9)
layout(matrix(1:11, nrow = 11, ncol = 1))
par(oma = c(3.1,3.5,1.3,1), mar = c(0,0,0,0))
years = 2005:2015

for(ii in 1:11){
  plot(-100, -100, xlim = c(0, nrow(GZ.data)), ylim = c(0, 5), xaxt = 'n', xaxs = 'i')
  CI = log10(write_CI(t(Arr.ts[ii,,])))
  CI[which(is.infinite(CI))] = 0
  time = 1:nrow(GZ.data)
  polygon(c(time, rev(time)), c(CI[1,], rev(CI[2,])), density = NA,
          col = rgb(0,0,1,0.3), border = NA)
  abline(v = seq(from = 1, to = 4017, by = 366), lwd = 1, col = rgb(0,0,0,1))
  tmp = log10(GZ.data$local.cases)
  tmp[which(is.infinite(tmp))] = 0
  lines(tmp, col = rgb(0,0,0,0.7), lwd = 2)
  ll = sample(x = 1:dim(Arr)[3], 1)
  tmp = log10(Arr.ts[ii,,ll])
  tmp[which(is.infinite(tmp))] = 0
  lines(tmp, lwd = 1.5, col = rgb(1,0,0,0.4))
  ll = sample(x = 1:dim(Arr)[3], 1)
  tmp = log10(Arr.ts[ii,,ll])
  tmp[which(is.infinite(tmp))] = 0
  lines(tmp, lwd = 1.5, col = adjustcolor(col = 'forestgreen', alpha.f = 0.4))
  ll = sample(x = 1:dim(Arr)[3], 1)
  tmp = log10(Arr.ts[ii,,ll])
  tmp[which(is.infinite(tmp))] = 0
  lines(tmp, lwd = 1.5, col = adjustcolor(col = 'navy', alpha.f = 0.4))
  mtext(text = years[ii], side = 2, cex = 0.75, line = 2)
}
# maybe add a rug or ceiling rug w/ importation timing
# like I wonder what was going on in 2006 in terms of importation 

# adding axes
tick = seq(from = 1, to = 4017, by = 366)
axis(side = 1, at = tick, labels = FALSE)
tmp = rep(0, 11)
names(tmp) = 2005:2015
par(new=TRUE)
barplot(tmp, ylim = c(0,2000), axes = F, yaxs = 'i', xaxs = 'i', space = rep(0,11), cex.names = 1.1)
# mtext(side = 1, outer = T, text = "Year", line = 2.05, cex = 0.9)

dev.off()

