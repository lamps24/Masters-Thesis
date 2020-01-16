# Research
Explanation of Code used for Simulation Study
Background
I am performing a simulation study to compare methods in their ability to estimate the causal impact of a hypothetical treatment when the assignment to the treatment group is not randomized. 


1.0	Generate_data

This program is a function that randomly generates the data. It is setup to create a design matrix with four variables, a user-specified sample size, treatment assignment based on the user-specified selection model coefficients, and a vector of the response variable based on the user-specified causal model coefficients. I also hand-picked a covariance matrix based on various research papers that I read prior to the study. 

The way that subjects are assigned to groups is based on a “propensity score.” The propensity score is the probability that a unit is assigned to the treatment group. If a subject has a propensity score of 75%, the function will assign it to the treatment group with 75% probability. The propensity scores are generated through an inverse logit model. The function returns the propensity scores as the “true” propensity scores, which are ultimately estimated as part of the propensity score matching process. 


2.0	PS_weight_cap

This code was a late addition as I discovered that units with very small propensity scores (<0.1%) were resulting in massive weights. ATE (average treatment effect) weights are derived as [1 / propensity score]. This code uses cross validation to pick a cap for the weights. It cycles through percentiles and picks the percentile that minimizes the out of sample MSE. 


3.0	Functions_robust_SE

This program is the meat of the code and takes the data and uses various methods to estimate the treatment effect. Variants of OLS, weighted OLS, and ridge penalized regression are used.


4.0	Simulation

The simulation pulls everything together and returns a matrix of results. There is a row for each of the combinations of treatment effect size and sample size. Each row contains (1) an estimate of the treatment effect from each method, (2) the average of the theoretical standard errors from each method, and (3) the simulation based standard error for each method. Due to the cross validation for the capped weights, the simulation can take upwards of 2-3 hours with 1,000 replications. Without cross validation it is much quicker. 

