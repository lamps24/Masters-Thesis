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

#data = gen.data(500, delta=1)

euclid = function(x, y, t, alpha=1, weighted=T) {
  
  y = as.vector(y)
  
  # estimate betas, [remove intercept]
  beta = as.vector(coef(lm(y ~ x)))[-1]
  
  # multiply each column of X by its beta
  newX = t(beta*t(x))
  
  # inverse of euclidean distances
  if (weighted==T)
  {
    wgts = 1 / ((as.matrix(dist(newX)))^alpha)  
  }
  if (weighted==F)
  {
    wgts = 1 / ((as.matrix(dist(x)))^alpha)
  }
  diag(wgts) = 0
  
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
  y.trt.cf = t*((t*y) - colSums(y * normalized.wgts.trt.cf)) #[trt effect on treated - cf]
  y.ctrl.cf = t.prime*(colSums(y * normalized.wgts.ctrl.cf) - (t.prime*y)) #[cf - trt effect on non-treated]
  
  # combine the vectors, calculate an average treatmet effect (ATE)
  ate = mean(y.trt.cf + y.ctrl.cf)
  
  return(list(euclid.delta=ate))

}




# euclid simulation
# which alpha does well?
alphas = c(1, 2)
ns = c(200, 500)
deltas = c(0.1, 0.5, 5)  
ks = c(1)
reps = 1000
which.betas = c(1,2,3,4,5,6,7,8,9,10) # controls which of the X1 - X10 covariates to keep in the estimations

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
        tmp.results = matrix(NA, nrow=reps, ncol=7)
        
        for (i in 1:reps) {
          
          # generate data
          data = gen.data(n=n, delta=delta)
          
          # euclid method
          euclid.dist = euclid(data$x, data$y, data$t, alpha=a)
          
          # put results in temp matrix
          tmp.results[i,1] = euclid.dist$euclid.delta
          
          cat("n=", n, "delta=", delta, "k=", k, "a=", a, "i=", i, "\n")

        }
        
        # calculate simulation means and SEs
        mean.euclid = mean(tmp.results[,1])
        se.euclid = sqrt(var(tmp.results[,1]))
        
        # put results into the matrix
        results[row.index,] = c(n, delta, a, mean.euclid, se.euclid)
        row.index = row.index + 1
      
      }
    }
  }
}
