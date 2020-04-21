# cross validation for PS weight cap selection
# function returns the cap that minimizes out of sample MSE
psw.cap.selection = function(data, weights, fmla.reg, cap.vec=seq(0.1, 1, by=0.05), K=5) {

  results = matrix(NA, length(cap.vec), 2)
  i = 1
  
  for (cap in cap.vec)
  {
    n = nrow(data)
    total.pred.err=0
    ind = sample(n)
    
    # cap the weights at the specified quantile
    cap.value = quantile(weights, probs=cap)
    wgt.cap = ifelse(weights > cap.value, cap.value, weights)
    data$wgt.cap = wgt.cap

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
  
    results[i,] = c(cap.vec[i], total.pred.err)  
    i = i + 1
    
  }

  min.row.index = which.min(results[,2])
  min.cap = results[min.row.index, 1]
  return(min.cap)
  
}




# test it out 
# psw.cap.selection(data=data, weights=w.est)