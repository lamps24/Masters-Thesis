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

data = gen.data(10, delta=5)

euclid = function(x, y, t) {
  
  # useful
  y = as.vector(data$y)
  
  # inverse of euclidean distances
  wgts = 1 / (as.matrix(dist(data$x)))
  diag(wgts) = 0
  
  # flip the 1s and 0s 
  t.prime = as.vector(1 * (data$t < 0.5))
  t = as.vector(data$t)
  
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

}


data = matrix(c(1,2,3,2,6,8,9,6,4), nrow=3, ncol=3)
as.matrix(dist(data))

euc = matrix(NA, nrow=3, ncol=3)
for (i in 1:3) 
{
  for (j in 1:3)
  {
    euc[i,j] = sqrt(sum((data[i,]-data[j,])^2))
  }
}



