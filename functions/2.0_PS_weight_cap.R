###############################################
# 2.0 - PS Weight Cap
# 
# Function called psw.cap.selection() that uses cross-validation to 
# select a maximum value for the weights that minimizes out of sample
# MSE. This is necessary because the propensity score estimation 
# leads to a handful of weights that are extremely large and leads to 
# unecessary variation in the weighted estimates. 
#
# The function takes the following arguments:
#   data - data object from the gen.data() function
#   weights - vector of weights from the gen.data() function
#   fmla.reg - syntax for causal model 
#   cap.vec - vector of percentiles to try for the caps
#   K - number of cross-validation folds
#
# The function returns the following:
#   min.cap - the cap that minimized the out of sample MSE. Given as percentile.
#
###############################################

psw.cap.selection = function(data, weights, fmla.reg, cap.vec=seq(0.1, 1, by=0.05), K=5) {
  
  # initialize results matrix, iteration count
  results = matrix(NA, length(cap.vec), 2)
  i = 1
  
  # loop through the candidate caps
  for (cap in cap.vec)
  {
    # initialize error value
    n = nrow(data)
    total.pred.err=0
    
    # "ind" holds the indices for cross validation observations (sample function randomly samples w/o replacement)
    ind = sample(n)
    
    # cap the weights at the specified quantile
    cap.value = quantile(weights, probs=cap)
    wgt.cap = ifelse(weights > cap.value, cap.value, weights)
    data$wgt.cap = wgt.cap

    # cross validation
    for(k in 1:K)
    {
      # indices for validation data
      leave.out = ind[(1 + floor((k-1)*n/K)):floor(k*n/K)]
      
      # training data
      train = data[-leave.out,]

      # validation data
      x.val = data[leave.out,]
      y.val = as.matrix(data[leave.out, 2])

      # evaluate [need to confirm whether to weight the prediction error]
      fit = lm(fmla.reg, data=train, weights=wgt.cap)
      y.val.hat = predict(fit, x.val)
      total.pred.err = total.pred.err + sum((y.val - y.val.hat)^2)
    }
    
    # enter the cap percentile alongside the prediction error  
    results[i,] = c(cap.vec[i], total.pred.err)  
    i = i + 1
    
  }

  # capture the index of the row that holds the minimum prediction error
  min.row.index = which.min(results[,2])
  min.cap = results[min.row.index, 1] # extract the cap percentile from row #min.row.index
  return(min.cap)
  
}