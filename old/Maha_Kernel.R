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

maha = function(x, y, t, Sigma, h) {
  
  y = as.vector(y)
  
  # compute mahalanobis distance matrix (distance of each point to each point, taking into account Sigma)
  maha.mat = matrix(NA, nrow=nrow(x), ncol=nrow(x))
  for (i in 1: nrow(x))
  {
    maha.mat[i,] = (mahalanobis(x, x[i,], Sigma)) / h
  }
  
  # apply kernel to maha distances
  wgts = (3/4) * (1 - maha.mat^2)
  wgts = ifelse(wgts < 0, 0, wgts)

  # distance to yourself is zero
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
ns = c(200)
deltas = c(0.5)  
hs = c(0.1)
reps = 1000

# initialize results matrix
rows = length(ns) * length(deltas) * length(hs)
results = matrix(NA, nrow=rows, ncol=5)
row.index = 1
set.seed(1234)

# simulation
for (n in ns) {
  for (delta in deltas) {
    for (h in hs) {

      # initialize temp results matrix
      tmp.results = matrix(NA, nrow=reps, ncol=1)
      
      for (i in 1:reps) {
        
        # generate data
        data = gen.data(n=n, delta=delta)
        
        # maha method
        maha.dist = maha(data$x, data$y, data$t, Sigma=data$Sigma, h=h)
        
        # put results in temp matrix
        tmp.results[i,1] = maha.dist$maha.delta
        
        cat("n=", n, "delta=", delta, "i=", i, "\n")
        
      }
      
      # calculate simulation means and SEs
      mean.maha = mean(tmp.results[,1])
      se.maha = sqrt(var(tmp.results[,1]))
      
      # put results into the matrix
      results[row.index,] = c(n, delta, h, mean.maha, se.maha)
      row.index = row.index + 1
    
    }
  }
}
