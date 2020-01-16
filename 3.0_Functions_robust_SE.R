###############################################
# 3.0 - Functions (updated for robust SE calculation)
# 
# Function called psm() that estimates the treatment effect from the data 
# generated in 1.0_Generate_data using various methods
#
# The function takes the following arguments:
#   X - the design matrix
#   Y - the vector of the response variable
#   t - the vector of treatment/control indicators
#   which.betas - which variables are known? For ex. c(1,2,4) means that we are trying to estimate 
#       the treatment effect without the third variable, X3
#   true.pscores - the vector of true probabilities of selection into the treatment group
#   robust - variance estimator. Needs to be modified to compute robust variance estimates due to the 
#       fact that the propensity scores are estimated and therefore induce added variability
#
# The function returns the estimates of treatment effect using each method, along with the nominal
# standard errors, which are used to compare to the simulation based estimate of standard error
#
###############################################
library(MatchIt)
library(glmnet)
library(jtools)

psm = function(X, Y, t, which.betas=c(1,2,3,4), true.pscores, robust="HC3") {
  
  # assemble dataframe 
  data = as.data.frame(cbind(t, Y, X))
  
  # formula needs to be implemented in this messy way in order to accept different model syntax
  betas = which.betas + 2 # X matrix begins in column 3 (X matrix doesn't include column of 1s)
  xnam = paste0("V", betas)
  (fmla.ps = as.formula(paste("V1 ~ ", paste(xnam, collapse= "+"))))
  (fmla.reg = as.formula(paste("V2 ~ V1 +", paste(xnam, collapse= "+"))))

  #########################
  ### PS MATCHING + OLS ###
  #########################
  
  # unmatched treatment effect (this is the estimate of treatment effect without doing any modeling)
  unmatched.delta = mean(Y[t==1]) - mean(Y[t==0]) 
  
  # test whether unmatched treatment effect significant
  unmatched.delta.p = t.test(Y[t==1], Y[t==0], conf.level=0.95)$p.value
  
  # match 
  out = matchit(fmla.ps, data=data, ratio=1)
  
  # add match flag to data
  data$w.m = out$weights
  
  # regression adjustment with matches
  lm = lm(fmla.reg, data=data, weights=w.m)
  delta = as.numeric(coef(lm)[2])
  delta.nom = as.numeric(coef(summary(lm))[2,2])
  
  ##########################
  ### PS WEIGHTING + OLS ###
  ##########################
  
  # compile data into a more usable format
  psw.true = cbind(Y, t)
  data$w.true = ifelse(psw.true[,2]==1, 1/true.pscores, 1/(1-true.pscores)) # true pscores
  
  # estimate weighted delta (true weights) - the estimate of treatment effect use ATE weights derived from pscores
  lm = lm(fmla.reg, data=data, weights=w.true)
  results = summ(lm, robust=robust)
  true.wgtd.delta = results$coeftable[2,1]
  true.wgtd.delta.nom = results$coeftable[2,2]

  # estimate weighted delta (true weights), capped at cross-validated best percentile
  best.cap = psw.cap.selection(data=data, weights=data$w.true, fmla.reg)
  cap = quantile(data$w.true, probs=(best.cap))
  data$w.true.cap = ifelse(data$w.true > cap, cap, data$w.true)
  lm = lm(fmla.reg, data=data, weights=w.true.cap)
  results = summ(lm, robust=robust)
  true.wgtd.delta.cap = results$coeftable[2,1]
  true.wgtd.delta.cap.nom = results$coeftable[2,2]
  
  
  # estimate pscores using logistic regression
  out.psw = glm(fmla.ps, data=data, family=binomial())
  est.pscores = predict(out.psw, type = "response")
  psw = cbind(t, est.pscores)
  
  # ATE weights:
  w.est = ifelse(psw[,1]==1, 1/psw[,2], 1/(1-psw[,2])) # estimated pscores
  
  # regression adjustment with weights
  lm = lm(fmla.reg, data=data, weights=w.est)
  results = summ(lm, robust=robust)
  wgtd.delta = results$coeftable[2,1]
  wgtd.delta.nom = results$coeftable[2,2]
  
  # regression adjustment with capped weights
  best.cap = psw.cap.selection(data=data, weights=w.est, fmla.reg)
  cap = quantile(w.est, probs=(best.cap))
  w.est.cap = ifelse(w.est > cap, cap, w.est)
  lm = lm(fmla.reg, data=data, weights=w.est.cap)
  results = summ(lm, robust=robust)
  wgtd.delta.cap = results$coeftable[2,1]
  wgtd.delta.cap.nom = results$coeftable[2,2]
  
  ################
  ### OLS only ###
  ################
  
  # regression without any data pre-processing
  lm = lm(fmla.reg, data=data)
  reg.delta = as.numeric(coef(lm)[2])
  reg.delta.nom = as.numeric(coef(summary(lm))[2,2])
  
  # ridge penalized estimate
  X = as.matrix(cbind(t, X))
  X = X[,which.betas]
  Y = as.matrix(Y)
  lambdas = 10^seq(3, -2, by = -.1)
  cv = cv.glmnet(X, Y, alpha = 0, lambda = lambdas)  
  m = glmnet(X, Y, alpha = 1, lambda = cv$lambda.min)
  ridge.delta = m$beta[1]
    
  
  return(list(unmatched.delta=unmatched.delta, reg.delta=reg.delta, reg.delta.nom=reg.delta.nom, ridge.delta=ridge.delta,
              matched.delta=delta, matched.delta.nom=delta.nom, wgtd.delta=wgtd.delta, wgtd.delta.nom=wgtd.delta.nom, 
              wgtd.delta.cap=wgtd.delta.cap, wgtd.delta.cap.nom=wgtd.delta.cap.nom,  
              true.wgtd.delta=true.wgtd.delta, true.wgtd.delta.nom=true.wgtd.delta.nom,
              true.wgtd.delta.cap=true.wgtd.delta.cap, true.wgtd.delta.cap.nom=true.wgtd.delta.cap.nom))
  
}
