
# gen data
data = gen.data.nl(200, delta=0.1)
X = data$x
t = data$t

# set up formula
which.betas=c(1,2,3,4)
betas = which.betas + 2 # X matrix begins in column 3 (X matrix doesn't include column of 1s)
xnam = paste0("V", betas)
(fmla.ps = as.formula(paste("V1 ~ ", paste(xnam, collapse= "+"))))

# set up data
df = data.frame(X, t)
colnames(df) = c("V3", "V4", "V5", "V6", "V1")

softplus = function(x) log(1 + exp(x))
nn = neuralnet(fmla.ps, data=df, hidden=c(3,2), act.fct = "logistic",
               linear.output = FALSE)
#preds = compute(nn, df)
#preds$net.result

plot(nn)