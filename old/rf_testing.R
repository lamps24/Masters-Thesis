source("C:/Users/lamps/Documents/Research/Code/NL_1.0_Generate_data.R")

data = gen.data.nl(200, delta=5)

mean(data$t)
mean(data$true.pscore)

truth = data$true.pscore

label = as.factor(data$t)
x = data$x
df = data.frame(x, label)

# run RF model, extract predicted probabilities
rf = randomForest(x=x, y=label, data=df)
est.pscores = rf$votes[,2]

mean(est.pscores)

