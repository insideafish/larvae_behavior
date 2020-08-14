# Dose-response functions. This is the SECOND file to run.

# Intercept-only, for comparison
Int = function(p,Y,C){
  g = p[1]
  s = exp(p[2])
  -sum(dnorm(x=Y,mean=g,sd=s,log=TRUE),na.rm=TRUE)}
parnames(Int) = c('g','s')

# Linear regression:
Linear = function(p,X,Y){
  b0 = p[1]
  b1 = p[2]
  s = exp(p[3])
  P = b0 + b1*X # equation of the line
  -sum(dnorm(x=Y,mean=P,sd=s,log=TRUE),na.rm=TRUE)}
parnames(Linear) = c('b0','b1','s')

# Quadratic:
Quadratic = function(p,X,Y){
  b0 = p[1] # maximum (or minimum)
  b1 = p[2] # shape parameter
  b2 = p[3] # X-value at the max/min
  s = exp(p[4])
  P = b0 + b1*(X-b2)^2 # equation of the line
  -sum(dnorm(x=Y,mean=P,sd=s,log=TRUE),na.rm=TRUE)}
parnames(Quadratic) = c('b0','b1','b2','s')

# Simpler Unimodal Dose-response: (5-parameter)
uniDR2 = function(p,X,Y){
  
  g = p[1] # max or min
  h = p[2] # direction (+ or -) and intensity
  #a1 = p[3] # exponential intercept term
  b1 = abs(p[3]) # exponential rate term (rate of increase before the hump)
  a2 = p[4] # exponential intercept term
  b2 = (p[5]) # exponential rate term (rate of decrease after the hump)
  s = exp(p[6]) # standard deviation of response (use exp to avoid negative values)
  P = g+h*(1 + b1*X)/(1+exp(-(a2 + b2*X))) # the equation for the curve
  # Assume a normal distribution for the likelihood (probably OK if data are log transformed), this is defining the actual function that you are fitting to the data
  
  -sum(dnorm(x=Y,mean=P,sd=s,log=TRUE),na.rm=TRUE)}
parnames(uniDR2) = c('g','h','b1','a2','b2','s')

plot_uniDR2 = function(p,X,Y){
  g = p[1] # max or min
  h = p[2] # direction (+ or -) and intensity
  #a1 = p[3] # exponential intercept term
  b1 = abs(p[3]) # exponential rate term (rate of increase before the hump)
  a2 = p[4] # exponential intercept term
  b2 = (p[5]) # exponential rate term (rate of decrease after the hump)
  s = exp(p[6]) # standard deviation of response (use exp to avoid negative values)
  Xtmp = seq(-1.5,2,0.1)
  P = g+h*(1 + b1*Xtmp)/(1+exp(-(a2 + b2*Xtmp))) # the equation for the curve
  plot(X,Y)
  lines(Xtmp,P)
}
parnames(plot_uniDR2) = c('g','h','b1','a2','b2','s') 


# Simpler Unimodal Dose-response: (5-parameter)
uniDR3 = function(p,X,Y){
  
  g = p[1] # max or min
  h = p[2] # direction (+ or -) and intensity
  #a1 = p[3] # exponential intercept term
  b1 = abs(p[3]) # exponential rate term (rate of increase before the hump)
  #a2 = p[4] # exponential intercept term
  b2 = (p[4]) # exponential rate term (rate of decrease after the hump)
  s = exp(p[5]) # standard deviation of response (use exp to avoid negative values)
  P = g+h*(1 + b1*X)/(1+b2*X) # the equation for the curve
  # Assume a normal distribution for the likelihood (probably OK if data are log transformed), this is defining the actual function that you are fitting to the data
  
  -sum(dnorm(x=Y,mean=P,sd=s,log=TRUE),na.rm=TRUE)}
parnames(uniDR3) = c('g','h','b1','b2','s')

# Unimodal Dose-response: (6-parameter)
uniDR = function(p,X,Y){
  
  g = p[1] # maximum
  h = p[2] # minimum
  a1 = p[3] # exponential intercept term
  b1 = abs(p[4]) # exponential rate term (rate of increase before the hump)
  a2 = p[5] # exponential intercept term
  b2 = -abs(p[6]) # exponential rate term (rate of decrease after the hump)
  s = exp(p[7]) # standard deviation of response (use exp to avoid negative values)
  P = g+h*(1 + exp(-(a1+b1*X)))/(1+exp(-(a2 + b2*X))) # the equation for the curve
  # Assume a normal distribution for the likelihood (probably OK if data are log transformed), this is defining the actual function that you are fitting to the data
  
  -sum(dnorm(x=Y,mean=P,sd=s,log=TRUE),na.rm=TRUE)}

parnames(uniDR) = c('g','h','a1','b1','a2','b2','s')

# Regular sigmoidal Dose-response:
regDR = function(p,X,Y){
  g = p[1] # maximum
  h = p[2] # minimum
  a = p[3] # exponential intercept term
  b = p[4] # exponential rate term
  s = exp(p[5]) # standard deviation of response; use exp to avoid negative values
  P = (h-g)/(1+exp(a + b*X))+g # the equation for the curve
  # Assume a normal distribution for the likelihood (probably OK if data are log transformed), this is defining the actual function that you are fitting to the data
  
  -sum(dnorm(x=Y,mean=P,sd=s,log=TRUE),na.rm=TRUE)}

parnames(regDR) = c('g','h','a','b','s')

