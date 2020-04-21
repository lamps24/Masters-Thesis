# Machine learning to estimate pscores
library(xgboost)

dta = gen.data(100, delta=1)
##################################################
### CROSS VALIDATION FOR HYPERPARAMETER TUNING ###
##################################################

xgb = function(X, Y, t, which.betas=c(1,2,3,4,5,6,7,8,9,10), true.pscores, rounds, depth, eta, folds, rep) {
  
  auc.fold = rep(0,nrep)
  pred = numeric(n)
  for(i in 1:rep) {
    
    train = sample(1:n, n*(folds-1)/folds)
    x.train = xgb.DMatrix(data = data[training,], label=[train])
    x.test = xgb.DMatrix(data = data[-train,], label=response[-train])
    x_hidden = xgb.DMatrix(data = hidden)
    y_train=response[train]
    y_test=response[-train]
    watchlist = list(train=x_train, test=x_test)
    
    m=xgb.train(data=x.train, max.depth=depth, eta=eta, nthread=2, nrounds=rounds, watchlist=watchlist, objective = "binary:logistic", eval_metric="auc")
    pred_test=predict(m,x_test)
    auc_fold[i]=auc(y_test,pred_test)
    predictions=predict(m,x_hidden)
    pred=(pred+predictions)
    
    ##### CONSTRUCTION
    data = cbind(X, t)
    xgb.data = xgb.DMatrix(data, label=data[,12])
    param = list(subsample = 0.8, max_depth = 3, eta = 0.11, verbose = 0, nthread = 2, objective = "binary:logistic", 
                 eval_metric = "auc")
    m.xgb = xgb.train(param, xgb.data, nrounds = 100)
    pscores = predict(m.xgb, X)
    pscores
    
    data(agaricus.train, package='xgboost')
    data(agaricus.test, package='xgboost')
    
    
    ## An 'xgboost' interface example:
    dta = gen.data(200, delta=1)
    data = cbind(dta$x, dta$t)
    bst = xgboost(data = data, label = data[,12], max_depth = 2, eta = 1, nthread = 2, nrounds = 2,
                   objective = "binary:logistic")
    pred = predict(bst, X)    
    
    
    
    
    
    

    
    
  }
  ans=c(mean(auc_fold), rounds, depth, eta, nfolds, pred/nrep)
  return(ans)
}
test_run=auc_cv_xgboost(100,2,.12,10, 200)

#Grid search for parameter tuning
rounds = c(100)
depth = c(2, 3, 4)
eta = c(.08, .11, .14, .17)
nfolds = c(5, 10)
nrep = c(50)

results = c()
for (a in 1:length(rounds)) {
  for (b in 1:length(depth)) {
    for (c in 1:length(eta)) {
      for (d in 1:length(nfolds)) {
        for (e in 1:length(nrep)) {
          results = append(results,auc_cv_xgboost(rounds[a],depth[b],eta[c],nfolds[d],nrep[e]))
        }
      }
    }
  }
}
scores = matrix(results,ncol=length(rounds)*length(depth)*length(eta)*length(nfolds)*length(nrep), byrow=FALSE)
write.csv(scores, file="C:\\Users\\lamps\\Documents\\STAT 8051 - Advanced Regression\\Final Project\\Predictions\\cv_results_after_eda.csv")


#########################
### FINAL PREDICTIONS ###
#########################

auc_cv_xgboost=function(rounds,depth,eta)
{
  auc_fold=rep(0,1)
  pred=rep(0,12002)
  for(i in 1:1)
  {
    x_full = xgb.DMatrix(data=data, label=response)
    x_hidden = xgb.DMatrix(data=hidden)
    watchlist = list(full=x_full)
    
    m=xgb.train(data=x_full, max.depth=depth, eta=eta, nthread=2, nrounds=rounds, watchlist=watchlist, objective = "binary:logistic", eval_metric="auc")
    predictions=predict(m, x_full)
    auc_fold[i]=auc(response, predictions)
    final_predictions=predict(m, x_hidden)
  }
  ans=c(mean(auc_fold), rounds, depth, eta, final_predictions)
  return(ans)
}
