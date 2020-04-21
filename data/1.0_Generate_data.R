###############################################
# 1.0 - Generate Data
# 
# Function called gen.data() that generates data from a 
# multivariate normal distribution. 
#
# The function takes the following arguments:
#   n - the sample size
#   beta.star - beta coefficients for the selection model
#   a.star - beta coefficients for the causal model
#   delta - true treatment effect size (this is what we are ultimately trying to estimate)
#
# The function returns the following:
#   x - the resulting X matrix
#   y - the resulting Y vector (the response variable)
#   t - indicator for whether subject was assigned to treatment (1) or control (0)
#   true.psocre - true propensity score (the probability that a unit was assigned to the treatment group)
#   Sigma - covariance matrix (note that covariance matrix is not an option for the user, I have pre-specified it)
#
###############################################

library(MASS)

gen.data = function(n, beta.star=c(0.25, 0.75, 0, 0.5), a.star=c(1, 2, 1, 0), delta) {
  
  #extract p (number of parameters, which includes intercept)
  p = length(beta.star) 
  
  # generate Sigma, the covariance matrix
  Sigma = matrix(0, nrow=p, ncol=p)
  Sigma[2,1] = 0.1; Sigma[1,2] = 0.1
  Sigma[3,1] = -0.1; Sigma[1,3] = -0.1
  Sigma[3,2] = 0.7; Sigma[2,3] = 0.7
  diag(Sigma) = 1

  # generate design matrix
  x = cbind(mvrnorm(n, mu=rep(0,(p)), Sigma))
  
  # estimate the true propensity score (pr(treatment))
  # additivity and linearity (this is the simplest case - the causal and selection model is linear)
  est = x %*% beta.star
  true.pscore = 1 / (1 + exp(-1*(est -2))) # results in about 13% of units being treated
  t = 1*(runif(n) < true.pscore)

  # generate y | t with normal error
  y = x%*%a.star + t*delta + rnorm(n)
  
  # values to be returned
  return(list(x=x, y=y, t=t, true.pscore=true.pscore, Sigma=Sigma))
  
}