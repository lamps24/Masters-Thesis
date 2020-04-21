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

# call user written functions
source("C:/Users/lamps/Documents/Research/Code/NL_1.0_Generate_data.R")
source("C:/Users/lamps/Documents/Research/Code/2.0_PS_weight_cap.R")
source("C:/Users/lamps/Documents/Research/Code/3.0_Functions_robust_SE.R")

# setup
ns = c(2000)
deltas = c(5)  
alpha = 0.05
reps = 50
which.betas = c(1,2,3,4) # controls which of the X1 - X4 covariates to keep in the estimations

# initialize results matrix
rows = length(ns) * length(deltas)
results = matrix(NA, nrow=rows, ncol=36)
row.index = 1
set.seed(12345)

# simulation
for (n in ns) {
  for (delta in deltas) {
    
    # initialize temp results matrix
    tmp.results = matrix(NA, nrow=reps, ncol=16)
    
    for (i in 1:reps) {
      
      # generate data
      data = gen.data.nl(n=n, delta=delta)
      
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
      tmp.results[i,15] = matches$nn.delta
      tmp.results[i,16] = matches$nn.delta.nom
      
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
    
    mean.nn = mean(tmp.results[,15])
    se.nn = sqrt(var(tmp.results[,15]))
    nom.se.nn = mean(tmp.results[,16])
    rmse.nn = abs(delta - mean.nn) + se.nn
    
    
    # put results into the matrix
    results[row.index,] = c(n, delta, 
                            mean.unmatched, se.unmatched, rmse.unmatched, mean.reg, se.reg, nom.se.reg, rmse.reg, 
                            mean.ridge, se.ridge, rmse.ridge, mean.matched, se.matched, nom.se.matched, rmse.matched, 
                            mean.wgtd, se.wgtd, nom.se.wgtd, rmse.wgtd, 
                            mean.wgtd.cap, se.wgtd.cap, nom.se.wgtd.cap, rmse.wgtd.cap, 
                            mean.true.wgtd, se.true.wgtd, nom.se.true.wgtd, rmse.true.wgtd,  
                            mean.true.wgtd.cap, se.true.wgtd.cap, nom.se.true.wgtd.cap, rmse.true.wgtd.cap,
                            mean.nn, se.nn, nom.se.nn, rmse.nn)
    row.index = row.index + 1
    
    # create boxplot
    # subset the tmp.results dataset to create side by side boxplots
    tmp.results2 = tmp.results[, c(2, 4, 5, 7, 9, 11, 13, 15)]
    
    # establish bounds of boxplot
    lbound = min(tmp.results2) - 0.01*delta
    ubound = max(tmp.results2) + 0.01*delta
    
    # name boxplot file
    outfile = paste("C:/Users/lamps/Documents/Research/Output/Plots/nl_boxplot", n, delta , ".pdf", sep="_")
    
    # open a pdf file
    pdf(outfile) 
    
    # create the plot
    par(mar=c(4,11,2,1))
    boxplot(tmp.results2, main="Performance of Estimators", ylim=c(lbound, ubound), xlab="Estimate (Red line is truth)",
            names=c("OLS", "Ridge", "Matched", "Weighted", "Weighted (capped)", "True Weighted", "True Weighted (capped)", "NN"),
            horizontal=T, las=1)
    abline(v=delta, col="Red")
    
    # close the pdf file
    dev.off() 
    
  }
}

# name the columns of the results matrix
colnames(results) = c("n", "delta", "Unmatched Est", "Unmatched SE", "Unmatched RMSE", "OLS Est", "OLS SE", "OLS nom SE", 
                      "OLS RMSE", "Ridge Est", "Ridge SE", "Ridge RMSE", "Matched Est", "Matched SE", "Matched nom SE", 
                      "Matched RMSE", "PSW Est", "PSW SE", "PSW nom SE", "PSW RMSE", "PSWC est", "PSWC SE", "PSWC nom SE",
                      "PSWC RMSE", "PSWT Est", "PSWT SE", "PSWT nom SE", "PSWT RMSE", "PSWCT Est", "PSWCT SE", "PSWCT nom SE", "PSCT RMSE",
                      "NN Est", "NN SE", "NN nom SE", "NN RMSE")

#write.csv(results, "C:/Users/lamps/Documents/Research/Output/NL_Results_021020.csv")

