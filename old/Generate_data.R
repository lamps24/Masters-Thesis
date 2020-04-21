# function that generates data
# n : sample size
# beta : vector of betas that influence treatment assignment (directly leads to pscores)
# a : vector of covariates that influence the outcome
library(MASS)

gen.data = function(n, beta.star=c(0.25, 0.75, 0, 0.5), a.star=c(1, 2, 1, 0), delta) {
  
  #extract p (number of parameters, which includes intercept)
  p = length(beta.star) 
  
  # generate Sigma
  Sigma = matrix(0, nrow=p, ncol=p)
  Sigma[2,1] = 0.1; Sigma[1,2] = 0.1
  Sigma[3,1] = -0.1; Sigma[1,3] = -0.1
  Sigma[3,2] = 0.7; Sigma[2,3] = 0.7
  diag(Sigma) = 1

  # design matrix
  x = cbind(mvrnorm(n, mu=rep(0,(p)), Sigma))
  
  # estimate the true propensity score (pr(treatment))
  # additivity and linearity 
  est = x %*% beta.star
  true.pscore = 1 / (1 + exp(-1*(est + log(0.1/0.9)))) # results in about 13% of units being treated
  t = 1*(runif(n) < true.pscore)

  # y | t
  y = x%*%a.star + t*delta + rnorm(n)
  
  return(list(x=x, y=y, t=t, true.pscore=true.pscore, Sigma=Sigma))
  
}





