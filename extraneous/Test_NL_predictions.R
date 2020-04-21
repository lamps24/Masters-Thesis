source("C:/Users/lamps/Documents/Research/Code/NL_1.0_Generate_data.R")
source("C:/Users/lamps/Documents/Research/Code/1.0_Generate_data.R")
require(xgboost)

# assemble data
data=gen.data.nl(2000, delta=1)
#data=gen.data(2000,delta=1)
mean(data$true.pscore)
t = data$t
X = data$x
label = as.factor(t)
df = data.frame(X, label)

# RANDOM FOREST
#X[,2] = X[,2]^2
rf = randomForest(x=X, y=label, data=df, ntree=1000)
#nn = randomForest(label~X[,1]+X[,2]+X[,3]+X[,4], data=df, ntree=500)
est.pscores.rf = rf$votes[,2]

# LOGISTIC REGRESSION
out.psw = glm(label ~ X1 + X2 + X3 + X4, data=df, family=binomial())
est.pscores.logit = predict(out.psw, type = "response")

# NEURAL NETWORK
softplus = function(x) log(1+exp(x))
nn = neuralnet(label ~ X1 + X2 + X3 + X4, data=df, hidden=c(1), act.fct = "logistic",
               linear.output = FALSE)
est.pscores.nn = compute(nn, df)$net.result[,2]

# XGBOOST
matrix = as.matrix(X)
xgb = xgboost(data=matrix, label = t, max.depth = 2, eta = 1,
                    nthread = 2, nrounds = 10, subsample=0.8, objective = "binary:logistic")
est.pscores.xgb = predict(xgb, X)


# check abs err
err.rf = mean(abs(data$true.pscore - est.pscores.rf))
err.logit = mean(abs(data$true.pscore - est.pscores.logit))
err.nn = mean(abs(data$true.pscore - est.pscores.nn))
err.xgb = mean(abs(data$true.pscore - est.pscores.xgb))

err.avg = mean(abs(data$true.pscore - (est.pscores.rf+est.pscores.logit+est.pscores.nn+est.pscores.xgb)/4))
err.avg.xgb.rf = mean(abs(data$true.pscore - (est.pscores.rf+est.pscores.xgb)/2))

c(err.rf, err.logit, err.nn, err.xgb, err.avg, err.avg.xgb.rf)

comp = cbind(data$true.pscore, est.pscores)