# calculates euclidean distance of each treated unit to every control unit
# a treated unit's estimate is the weighted average of the control units, 
#   where the weights are inversely proportional to the maha distance to the treated unit
#
# the treatment effect is the average 
#
# STEPS
#
# (1) compute mahalanobis distance of every control unit to every treatment
# (2) compute mahalanobis distance of every treatment unit to every control
# (3) estimate counterfactual for every treatment unit
# (4) estimate counterfactual for every control unit
# (5) compute ATE 
##################################################################################

#data = gen.data(200, delta=5)

maha = function(x, y, t, Sigma, alpha=1, weighted=T, cutoff=0.5) {
  
  y = as.vector(y)
  
  # estimate betas, [remove intercept]
  beta = as.vector(coef(lm(y ~ x)))[-1]
  
  # multiply each column of X by its beta
  newX = t(beta*t(x))
  
  # inverse of maha distances
  if (weighted==T)
  {
    # compute mahalanobis distance matrix (distance of each point to each point, taking into account Sigma)
    maha.mat = matrix(NA, nrow=nrow(x), ncol=nrow(x))
    for (i in 1: nrow(x))
    {
      maha.mat[i,] = mahalanobis(x, x[i,], Sigma)
    }
    wgts = 1 / (maha.mat^alpha)  
  }
  if (weighted==F)
  {
    # compute mahalanobis distance matrix (distance of each point to each point, taking into account Sigma)
    maha.mat = matrix(NA, nrow=nrow(x), ncol=ncol(x))
    for (i in 1: nrow(x))
    {
      maha.mat[i] = mahalanobis(x, x[i,], Sigma)
    }
    wgts = 1 / (maha.mat^alpha)  
  }
  
  # distance to yourself is zero
  diag(wgts) = 0
  
  # remove observations that are too far away (small wgts)
  cut = as.numeric(quantile(wgts, cutoff))
  wgts = ifelse(wgts < cut, 0, wgts)
  
  
  # flip the 1s and 0s 
  t.prime = as.vector(1 * (t < 0.5))
  t = as.vector(t)
  
  # weight matrices for counterfactual components (cf)
  wgts.trt.cf = t.prime * wgts
  wgts.ctrl.cf = t * wgts
  
  # normalize the weights so that they add to 1
  normalized.wgts.trt.cf = wgts.trt.cf / colSums(wgts.trt.cf)[col(wgts.trt.cf)]
  normalized.wgts.ctrl.cf = wgts.ctrl.cf / colSums(wgts.ctrl.cf)[col(wgts.ctrl.cf)]
  
  # calculate counterfactuals (sum of (wgt  * y))
  y.trt.cf = t * (y - colSums(y * normalized.wgts.trt.cf)) #[trt effect on treated - cf]
  y.ctrl.cf = t.prime * (colSums(y * normalized.wgts.ctrl.cf) - y) #[cf - trt effect on non-treated]
  
  # set NaNs to zero, if they exist
  y.trt.cf = ifelse(y.trt.cf=="NaN", 0, y.trt.cf)
  y.ctrl.cf = ifelse(y.ctrl.cf=="NaN", 0, y.ctrl.cf)
  
  # combine the vectors, calculate an average treatmet effect (ATE)
  ate = mean(y.trt.cf + y.ctrl.cf)
  
  return(list(maha.delta=ate))
  
}




# maha simulation
# which alpha does well?
alphas = c(1, 5, 50, 100)
ns = c(200)
deltas = c(0.5, 5)  
ks = c(1)
reps = 1000
which.betas = c(1,2,3,4) # controls which of the X1 - X10 covariates to keep in the estimations

# initialize results matrix
rows = length(ns) * length(deltas) * length(ks) * length(alphas)
results = matrix(NA, nrow=rows, ncol=5)
row.index = 1
set.seed(1234)

# simulation
for (n in ns) {
  for (delta in deltas) {
    for (k in ks) {
      for (a in alphas) {
        
        # initialize temp results matrix
        tmp.results = matrix(NA, nrow=reps, ncol=1)
        
        for (i in 1:reps) {
          
          # generate data
          data = gen.data(n=n, delta=delta)
          
          # maha method
          maha.dist = maha(data$x, data$y, data$t, alpha=a, Sigma=data$Sigma)
          
          # put results in temp matrix
          tmp.results[i,1] = maha.dist$maha.delta
          
          cat("n=", n, "delta=", delta, "k=", k, "a=", a, "i=", i, "\n")
          
        }
        
        # calculate simulation means and SEs
        mean.maha = mean(tmp.results[,1])
        se.maha = sqrt(var(tmp.results[,1]))
        
        # put results into the matrix
        results[row.index,] = c(n, delta, a, mean.maha, se.maha)
        row.index = row.index + 1
        
      }
    }
  }
}
