# estimates of the serial interval distribution were derived from...
#
# LOCAL AND REGIONAL SPREAD OF CHIKUNGUNYA FEVER IN THE AMERICAS
# S Cauchemez, M Ledrans, C Poletto, P Quenel, H de Valk, V Colizza, P Y Boelle
# Eurosurveillance 2014
#
# ...who report that the distribution has a mean of 23, standard devition of 6,
# and in a plot they show it looks somewhat like a gamma

# by the method of moments, the parameters of the gamma are
param.1 = (23/6)^2
param.2 = 23/(6^2)

# as one can see, most of the mass falls in weeks 3-5, with the amount in each
weights.serialinterval =
  sapply(1:49,function(ww)
    integrate(function(t){
      pgamma(t+1,param.1,param.2)-pgamma(t,param.1,param.2)},
      (ww-1)*1,ww*1)$value/1)
weights.serialinterval = weights.serialinterval / sum(weights.serialinterval)

# get rid of days that are below 1% probability
weights.serialinterval[weights.serialinterval < .01] = 0
weights.serialinterval = weights.serialinterval / sum(weights.serialinterval)
