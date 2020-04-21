# testing
data = gen.data(n=2000, delta=5)
ps = psm(data$x, data$y, data$t, data$true.pscore)
ps

data = gen.data(n=2000, delta=5)
ps = psm(data$x, data$y, data$t, which.betas=c(1,2,3,4,5,6,7,8,9,10), data$true.pscore)
ps




# test
which.betas=c(1,2,3,4,5,6,7,8,9,10)
dta = gen.data(n=2000, delta=5)
data = as.data.frame(cbind(dta$t, dta$y, dta$x, data$true.pscore))
betas=which.betas+2

xnam = paste0("V", betas)
(fmla.ps = as.formula(paste("V1 ~ ", paste(xnam, collapse= "+"))))

out = matchit(fmla.ps, data=data, ratio=1)





# assemble dataframe 
data = as.data.frame(cbind(data$t, data$y, data$x))

# formula
betas = which.betas + 2 # X matrix begins in column 3
xnam = paste0("V", betas)
(fmla.ps = as.formula(paste("V1 ~ ", paste(xnam, collapse= "+"))))
(fmla.reg = as.formula(paste("V2 ~ V1 +", paste(xnam, collapse= "+"))))

# match using logistic regression
if (method == "logit") {
  
  ###################
  ### PS MATCHING ###
  ###################
  
  # unmatched treatment effect
  unmatched.delta = mean(Y[t==1]) - mean(Y[t==0]) 
  
  # test whether unmatched treatment effect significant
  unmatched.delta.p = t.test(Y[t==1], Y[t==0], conf.level=0.95)$p.value
  
  # match 
out$weights







out = matchit(data[,1] ~ data[,3] + data[,4] + data[,5] + data[,6] + data[,7] + 
                data[,8] + data[,9] + data[,10] + data[,11] + data[,12], 
              data=data, method="nearest", distance="logit",
              replace=FALSE, ratio=k)







lm(fmla, data=data)
matchit(fmla, data=data, method="nearest", distance="logit",
        replace=FALSE, ratio=1)































# propensity score matching
psm = function(X, Y, t, which.betas=c(1,2,3,4,5,6,7,8,9,10), k=1, method="logit") {
  
  # assemble dataframe 
  data = as.data.frame(cbind(t, Y, X))
  
  # formula
  betas = which.betas + 2 # X matrix begins in column 3
  xnam = paste0("V", betas)
  (fmla.ps = as.formula(paste("V1 ~ ", paste(xnam, collapse= "+"))))
  (fmla.reg = as.formula(paste("V2 ~ V1 +", paste(xnam, collapse= "+"))))
  
  # match using logistic regression
  if (method == "logit") {
    
    ###################
    ### PS MATCHING ###
    ###################
    

    # match 
    out = matchit(fmla.ps, data=data, ratio=k)
  }
  return(out)
}


