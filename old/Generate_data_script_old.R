# function that generates data
# n : sample size
# beta : vector of betas that influence treatment assignment (directly leads to pscores)
# a : vector of covariates that influence the outcome
library(MASS)

gen.data = function(n, beta.star=c(0.8,-0.25,0.6,-0.4,-0.8,-0.5,0.7,0,0,0), 
                    a.star=c(0.3,-0.36,-0.73,-0.2,0,0,0,0.71,-0.19,0.26), 
                    delta) {
  
  #extract p (number of parameters, which includes intercept)
  p = length(beta.star) 
  
  # generate Sigma
  Sigma = matrix(0, nrow=p, ncol=p)
  Sigma[5,1] = 0.2; Sigma[1,5] = 0.2
  Sigma[6,2] = 0.9; Sigma[2,6] = 0.9
  Sigma[8,3] = 0.2; Sigma[3,8] = 0.2
  Sigma[9,4] = 0.9; Sigma[4,9] = 0.9
  diag(Sigma) = 1

  # design matrix
  x = cbind(mvrnorm(n, mu=rep(0,(p)), Sigma))
  
  # estimate the true propensity score (pr(treatment))
  # additivity and linearity 
  est = x %*% beta.star
  true.pscore = 1 / (1 + exp(-1*(est + log(0.1/0.9)))) # results in about 17% of units being treated
  t = 1*(runif(n) < true.pscore)

  # y | t
  y = x%*%a.star + t*delta + rnorm(n)
  
  return(list(x=x, y=y, t=t, true.pscore=true.pscore))
  
}


qr.solve( t(xmat) %*% wmat %*% xmat ) %*% t(xmat) %*% wmat %*% ymat



