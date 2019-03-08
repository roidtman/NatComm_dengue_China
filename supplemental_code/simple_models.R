

##### set up libraries
require(mgcv)
require(fda)
require(lubridate)
require(scam)

# =============================================================== # 
###### load & set up data
source('runFunctions_bayes.R')
source('runSerialinterval.R')

mosq = read.csv('../data/guangzhou_vector_05-15.csv')
GZ.data = read.csv('../data/guangzhou_case_weather_05-15.csv')

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


# aggregate GZ.data local cases weekly
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


# =============================================================== # 
###### fit poisson spline

to_keep = which(GZ.data$local.cases != 0 & GZ.data$I.old != 0)

df_poisson = data.frame(I_new = GZ.data$local.cases[to_keep], 
                        I_old = GZ.data$I.old[to_keep] + GZ.data$imported.cases[to_keep],
                        T_mean = GZ.data$T_mean[to_keep],
                        mosq_spline = mosq.spline[,1][to_keep]
                        ,
                        time = (1:nrow(GZ.data))[to_keep]
                        )


model_gam = scam(I_new ~ I_old + s(T_mean,bs='cv') 
                 + s(mosq_spline,bs='mpi')
                 + s(time)
                 , 
                 family = 'poisson', data = df_poisson)

AIC(model_gam)

# need data.frame names to be the same as in df_poisson
df_predict = data.frame(I_new = GZ.data$local.cases,
                        I_old = GZ.data$I.old + GZ.data$imported.cases,
                        T_mean = GZ.data$T_mean,
                        mosq_spline = mosq.spline[,1]
                        ,
                        time = 1:nrow(GZ.data)
                        )

# use model to predict local cases
pred_gam = predict(model_gam, newdata = df_predict, se.fit = T)

# calculate correlation between model fit and observed local case data
corr = cor(exp(pred_gam$fit), GZ.data$local.cases)


# =============================================================== # 
## log figure 

exp_fit_upper = exp(pred_gam$fit + pred_gam$se.fit)
exp_fit_lower = exp(pred_gam$fit - pred_gam$se.fit)
log_fit_upper = log10(exp_fit_upper)
log_fit_lower = log10(exp_fit_lower)
log_cases = log10(GZ.data$local.cases)

pdf('../output/fig_A1.pdf',
    height = 3, width = 6)

layout(1)
par(oma = c(3,2.5,1,1), mar = c(0.5, 1.5, 0,0))

plot(-100, -100, ylim = c(min(log_fit_upper, log_fit_lower), max(log_fit_upper, log_fit_lower)),
     xlim = c(0, nrow(GZ.data)), xaxt = 'n', las = 1)

lines(log_cases, col = rgb(0,0,0,0.5), lwd = 1)

polygon(c(1:nrow(df_predict), rev(1:nrow(df_predict))),
        c(log_fit_lower, rev(log_fit_upper)),
        border = rgb(0,0,1,0.3),
        col = rgb(0,0,1,0.5))

text(x = nrow(df_predict) * 0.965, y = max(log_fit_upper, log_fit_lower)*0.95, 
     labels = expression(rho* ' = 0.269'), cex = 0.8)

tick = seq(from = 1, to = 4017, by = 366)
axis(side = 1, at = tick, labels = FALSE)
axis(side = 1, at = tick + (366/2), tick = FALSE, labels = 2005:2015)
mtext(expression(log[10]*'(incidence)'), side=2, line=0.9, outer = T)
mtext('Year', side=1, line=2.1)

dev.off()
