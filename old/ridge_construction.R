

dta = gen.data(n=100, delta=1)
data = as.data.frame(cbind(dta$t, dta$Y, dta$X))
library(glmnet)

# ridge penalization
X = as.matrix(cbind(dta$t, dta$x))
Y = as.matrix(dta$y)
lambdas = 10^seq(3, -2, by = -.1)
cv = cv.glmnet(X, Y, alpha = 0, lambda = lambdas)  
m = glmnet(X, Y, alpha = 0, lambda = cv$lambda.min)
est = m$beta$1
