# call user written functions
source("C:/Users/lamps/Documents/Research/Code/1.0_Generate_data.R")
source("C:/Users/lamps/Documents/Research/Code/2.0_PS_weight_cap.R")
source("C:/Users/lamps/Documents/Research/Code/3.0_Functions_robust_SE.R")

dta = gen.data(n=200, delta=1)
data = as.data.frame(cbind(dta$t, dta$y, dta$x))

which.betas=c(1,2,3,4)
betas = which.betas + 2 # X matrix begins in column 3 (X matrix doesn't include column of 1s)
xnam = paste0("V", betas)
(fmla.ps = as.formula(paste("V1 ~ ", paste(xnam, collapse= "+"))))
(fmla.reg = as.formula(paste("V2 ~ V1 +", paste(xnam, collapse= "+"))))


# weights from logistic regression
out.psw = glm(fmla.ps, data=data, family=binomial())
est.pscores = predict(out.psw, type = "response")
psw = cbind(dta$t, est.pscores)
w.est = ifelse(psw[,1]==1, 1/psw[,2], 1/(1-psw[,2])) # estimated pscores


# regression adjustment with capped weights
best.cap = psw.cap.selection(data=data, weights=w.est, fmla.reg)
cap = quantile(w.est, probs=(best.cap))
w.est.cap = ifelse(w.est > cap, cap, w.est)


cbind(w.est, w.est.cap)
max(w.est.cap)
