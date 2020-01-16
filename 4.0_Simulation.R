###############################################
# 4.0 - The Simulation Study
# 
# Simulation study. "reps" iterations. Loops through ns and deltas "reps" times.
# 
# Input values are as follows:
#   ns - sample sizes 
#   deltas - treatment effect sizes
#   reps - number of repetitions, or iterations, to perform
#   alpha - parameter for the failed "euclidean/mahalonobis" method
#   which.betas - controls which variables are known to the simulation study
#
###############################################

# setup
ns = c(200, 500, 2000)
deltas = c(0.1, 0.5, 5)  
alpha = 0.05
reps = 1000
which.betas = c(1,2,3) # controls which of the X1 - X10 covariates to keep in the estimations

# initialize results matrix
rows = length(ns) * length(deltas)
results = matrix(NA, nrow=rows, ncol=32)
row.index = 1
set.seed(1234)

# simulation
for (n in ns) {
  for (delta in deltas) {

    # initialize temp results matrix
    tmp.results = matrix(NA, nrow=reps, ncol=14)
     
    for (i in 1:reps) {
      
      # generate data
      data = gen.data(n=n, delta=delta)
      
      # generate matches
      matches = psm(data$x, data$y, data$t, which.betas, data$true.pscore)
      
      # euclid method - did not work very well. Turned out to be very similar to "kernel matching"
      #euclid.dist = euclid(data$x, data$y, data$t, alpha=5, weighted=T)
      
      # put results in temp matrix
      tmp.results[i,1] = matches$unmatched.delta
      tmp.results[i,2] = matches$reg.delta
      tmp.results[i,3] = matches$reg.delta.nom
      tmp.results[i,4] = matches$ridge.delta
      tmp.results[i,5] = matches$matched.delta
      tmp.results[i,6] = matches$matched.delta.nom
      tmp.results[i,7] = matches$wgtd.delta
      tmp.results[i,8] = matches$wgtd.delta.nom
      tmp.results[i,9] = matches$wgtd.delta.cap
      tmp.results[i,10] = matches$wgtd.delta.cap.nom
      tmp.results[i,11] = matches$true.wgtd.delta
      tmp.results[i,12] = matches$true.wgtd.delta.nom
      tmp.results[i,13] = matches$true.wgtd.delta.cap
      tmp.results[i,14] = matches$true.wgtd.delta.cap.nom
      
      # print results so that you can track progress of simulation (can take upwards of 4 hours due to cross validation)
      cat("n=", n, "delta=", delta, "i=", i, "\n")
      
    }
    
    # calculate simulation means and SEs
    mean.unmatched = mean(tmp.results[,1])
    se.unmatched = sqrt(var(tmp.results[,1]))
    rmse.unmatched = abs(delta - mean.unmatched) + se.unmatched
    
    mean.reg = mean(tmp.results[,2])
    se.reg = sqrt(var(tmp.results[,2]))
    nom.se.reg = mean(tmp.results[,3])
    rmse.reg = abs(delta - mean.reg) + se.reg
    
    mean.ridge = mean(tmp.results[,4])
    se.ridge = sqrt(var(tmp.results[,4]))
    rmse.ridge = abs(delta - mean.ridge) + se.ridge
    
    mean.matched = mean(tmp.results[,5])
    se.matched = sqrt(var(tmp.results[,5]))
    nom.se.matched = mean(tmp.results[,6])
    rmse.matched = abs(delta - mean.matched) + se.matched
    
    mean.wgtd = mean(tmp.results[,7])
    se.wgtd = sqrt(var(tmp.results[,7]))
    nom.se.wgtd = mean(tmp.results[,8])
    rmse.wgtd = abs(delta - mean.wgtd) + se.wgtd
    
    mean.wgtd.cap = mean(tmp.results[,9])
    se.wgtd.cap = sqrt(var(tmp.results[,9]))
    nom.se.wgtd.cap = mean(tmp.results[,10])
    rmse.wgtd.cap = abs(delta - mean.wgtd.cap) + se.wgtd.cap
    
    mean.true.wgtd = mean(tmp.results[,11])
    se.true.wgtd = sqrt(var(tmp.results[,11]))
    nom.se.true.wgtd = mean(tmp.results[,12])
    rmse.true.wgtd = abs(delta - mean.true.wgtd) + se.true.wgtd
    
    mean.true.wgtd.cap = mean(tmp.results[,13])
    se.true.wgtd.cap = sqrt(var(tmp.results[,13]))
    nom.se.true.wgtd.cap = mean(tmp.results[,14])
    rmse.true.wgtd.cap = abs(delta - mean.true.wgtd.cap) + se.true.wgtd.cap
  
    
    # put results into the matrix
    results[row.index,] = c(n, delta, 
                            mean.unmatched, se.unmatched, rmse.unmatched, mean.reg, se.reg, nom.se.reg, rmse.reg, 
                            mean.ridge, se.ridge, rmse.ridge, mean.matched, se.matched, nom.se.matched, rmse.matched, 
                            mean.wgtd, se.wgtd, nom.se.wgtd, rmse.wgtd, 
                            mean.wgtd.cap, se.wgtd.cap, nom.se.wgtd.cap, rmse.wgtd.cap, 
                            mean.true.wgtd, se.true.wgtd, nom.se.true.wgtd, rmse.true.wgtd,  
                            mean.true.wgtd.cap, se.true.wgtd.cap, nom.se.true.wgtd.cap, rmse.true.wgtd.cap)
    row.index = row.index + 1
          
  }
}
  


