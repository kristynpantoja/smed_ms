home = "/Users/kristyn/Documents/research/smed_ms"
functions_home = paste(home, "/functions", sep="")


###################
# Loading Designs #
###################

load(paste(home, "/designs_slides.RData", sep = ""))
# object loaded is list "all_designs"
  
# parameters/settings
mean_beta0 = c(0, 0, 0) # null model prior mean on beta
mean_beta1 = c(0, 0.2, -0.2, 0.2, -0.2) # alternative model prior mean on beta
var_e = 0.005 # variance on error
var_beta0 = diag(rep(var_e, 3)) # null model prior variance on beta
var_beta1 = diag(rep(var_e, 5)) # alternative model prior variance on beta
xmin = 0
xmax = 1
f0 = function(x) mean_beta0[1] + mean_beta0[2] * x[1] + mean_beta0[3] * x[2] # null
f1 = function(x) mean_beta1[1] + mean_beta1[2] * x[1] + mean_beta1[3] * x[1]^2 + mean_beta1[4] * x[2] + mean_beta1[3] * x[2]^2
N = 100
type = c(4, 5) # type (4,5) model is my indicator to the function of what kind of models are being compared so that marginal variance is correctly computed
p = 3


# for one-at-a-time algorithm, k = 1
numCandidates = 10e3
k = c(1, 4, 50)
one_at_a_time_k1 = all_designs$one_at_a_time_k1

# for one-at-a-time algorithm, k = 1
one_at_a_time_k4 = all_designs$one_at_a_time_k4

# for one-at-a-time algorithm, k = 50
one_at_a_time_k50 = all_designs$one_at_a_time_k50

# for fast algorithm, S = 5
S = 5
fast_S5 = all_designs$fast_S5

# D-optimal
#rep1 = rep(c(0, 1), each = 25); rep2 = rep(0, each = 50); rep3 = rep(1, each = 50)
#doptimal = rbind(cbind(rep1, rep2), cbind(rep1, rep3))
doptimal = all_designs$doptimal

# Space-filling (grid)
#axis1 = seq(from = xmin, to = xmax, length.out = 10); axis2 = seq(from = xmin, to = xmax, length.out = 10)
#space_filing = cbind(rep(axis1, each = 10), rep(axis2, times = 10))
space_filling = all_designs$space_filling

# Random
#set.seed(1)
#axis1 = runif(100, xmin, xmax); axis2 = runif(100, xmin, xmax)
#random_design = cbind(axis1, axis2)
random_design = all_designs$random_design

####################
# Evaluate Designs #
####################



library(expm)
library(matrixStats)
library(scatterplot3d)
library(knitr)
library(mvtnorm)
library(here)
source(paste(functions_home, "/construct_design_matrix.R", sep = ""))
source(paste(functions_home, "/posterior_variance.R", sep = ""))
source(paste(functions_home, "/variance_marginal_y.R", sep = ""))
source(paste(functions_home, "/wasserstein_distance.R", sep = ""))
source(paste(functions_home, "/charge_function_q.R", sep = ""))
source(paste(functions_home, "/total_potential_energy_criteria.R", sep = ""))
source(paste(functions_home, "/fast_criteria.R", sep = ""))
source(paste(functions_home, "/oneatatime_criteria.R", sep = ""))
source(paste(functions_home, "/simulate_y.R", sep = ""))
source(paste(functions_home, "/postprob_hypotheses.R", sep = ""))
source(paste(functions_home, "/postmean_mse_closedform.R", sep = ""))
source(paste(functions_home, "/postmean_mse_mc.R", sep = ""))

#

calculateEvals = function(D, N, true_beta, mean_beta0, mean_beta1, var_beta0, var_beta1,
                          var_e, f0, f1, true_type, type, var_margy0 = NULL, var_margy1 = NULL, 
                          p, numSims){
  
  # posterior variance of beta for the alternative hypothesis
  v = postvar(D, N, var_e, var_beta1, type[2])
  # Total Potential Energy criterion
  TPE = totalPE_2d(D, N, mean_beta0, mean_beta1, var_beta0, var_beta1, var_e, 
                   f0, f1, type, var_margy0, var_margy1, p)
  # Fast Algorithm criterion
  cFast = crit_fast_2d(D, N, mean_beta0, mean_beta1, var_beta0, var_beta1, var_e, 
                       f0, f1, type, var_margy0, var_margy1, p)
  # One-at-a-Time Algorithm criterion
  c1atT = crit_1atatime_2d(D, N, k = 4, mean_beta0, mean_beta1, var_beta0, var_beta1, 
                           var_e, f0, f1, type, var_margy0, var_margy1, p)
  
  # sumstats
  means = apply(D, 2, mean)
  sds = apply(D, 2, sd)
  
  # expected posterior probabilities of hypotheses and bayes factors
  exppostprobs = calcExpPostProbH_2d(D, N, true_beta, mean_beta0, mean_beta1, 
                                     var_beta0, var_beta1, var_e, numSims, true_type, type)
  exppostmeanMSE = calcExpPostMeanMSE(D, N, true_beta, mean_beta0, mean_beta1, 
                                      var_beta0, var_beta1, var_e,
                                      numSims, true_type, type)
  
  # export
  evals = list("v" = v, "TPE" = TPE, "cFast" = cFast, "c1atT" = c1atT, 
               "mean" = means, "sd" = sds, "exppostprobs" = exppostprobs, "expMSEs" = exppostmeanMSE)
  return(evals)
}

# We want to find an example where MEDs do best i.t.o. Expected Posterior Probabilities of Hypotheses
# i.e. where MEDs give the highest value of E[P(H1|Y,D)|BT]

# Expected Posterior Probabilities of Hypotheses, BT = c(0.0, 0.4, -0.4, 0.4, -0.4)
# here, space-filling design does best, and d-optimal design does worst.
# does it make sense that the MEDs do in-between, though?
# why does that happen?

# with steeper values, e.g. BT = c(0, 1, -1, 1, -1) the true model is way off from H1, 
# and all of them do fine except for the d-optimal design. not much to analyze in that case.

numSims = 100
set.seed(123)

true_beta_v1 = c(0.0, 0.4, -0.4, 0.4, -0.4)
true_type_v1 = 5
calculateEvals_v1 = function(D){
  calculateEvals(D, N, true_beta_v1, mean_beta0, mean_beta1, 
                 var_beta0, var_beta1, var_e, f0, f1, true_type_v1, type, 
                 var_margy0 = NULL, var_margy1 = NULL, p, numSims)
}
k1_evals_v1 = calculateEvals_v1(one_at_a_time_k1)
k4_evals_v1 = calculateEvals_v1(one_at_a_time_k4)
k50_evals_v1 = calculateEvals_v1(one_at_a_time_k50)
fast_evals_v1 = calculateEvals_v1(fast_S5)
dopt_evals_v1 = calculateEvals_v1(doptimal)
space_evals_v1 = calculateEvals_v1(space_filling)
random_evals_v1 = calculateEvals_v1(random_design)

exppostprobs_v1 = data.frame(cbind(k1_evals_v1$exppostprobs, 
                                   k4_evals_v1$exppostprobs, 
                                   k50_evals_v1$exppostprobs, 
                                   fast_evals_v1$exppostprobs, 
                                   dopt_evals_v1$exppostprobs, 
                                   space_evals_v1$exppostprobs, 
                                   random_evals_v1$exppostprobs))
colnames(exppostprobs_v1) = c("1atT,k=1","1atT,k=4","1atT,k=50","Fast","DOptimal","Space","Random")
rownames(exppostprobs_v1) = c("E[P(H0|Y,D)|BT]", "E[P(H1|Y,D)|BT]", "E[BF01|BT]")
round(exppostprobs_v1, 5)

# Expected Posterior Probabilities of Hypotheses, c(0,2, 0, 0)
# to seei if other type works

numSims = 100
set.seed(123)

true_beta_v2 = c(0.2, 0, 0)
true_type_v2 = 4
calculateEvals_v2 = function(D){
  calculateEvals(D, N, true_beta_v2, mean_beta0, mean_beta1, 
                 var_beta0, var_beta1, var_e, f0, f1, true_type_v2, type, 
                 var_margy0 = NULL, var_margy1 = NULL, p, numSims)
}
k1_evals_v2 = calculateEvals_v2(one_at_a_time_k1)
k4_evals_v2 = calculateEvals_v2(one_at_a_time_k4)
k50_evals_v2 = calculateEvals_v2(one_at_a_time_k50)
fast_evals_v2 = calculateEvals_v2(fast_S5)
dopt_evals_v2 = calculateEvals_v2(doptimal)
space_evals_v2 = calculateEvals_v2(space_filling)
random_evals_v2 = calculateEvals_v2(random_design)

expMSEs_v2 = data.frame(cbind(k1_evals_v2$expMSEs, 
                                   k4_evals_v2$expMSEs, 
                                   k50_evals_v2$expMSEs, 
                                   fast_evals_v2$expMSEs, 
                                   dopt_evals_v2$expMSEs, 
                                   space_evals_v2$expMSEs, 
                                   random_evals_v2$expMSEs))
colnames(expMSEs_v2) = c("1atT,k=1","1atT,k=4","1atT,k=50","Fast","DOptimal","Space","Random")
rownames(expMSEs_v2) = c("E[MSE|H0]", "E[MSE|H1]")
round(expMSEs_v2, 5)

# Expected Posterior Probabilities of Hypotheses, BT = c(0, 0.2, 0, 0.2, 0)
# here, we expect d-optimal design to do best, and it does. MEDs do better than the space-filling design
# and random design. This might be worth noting.

numSims = 100
set.seed(123)

true_beta_v3 = c(0.2, 0, 0, 0, 0)
true_type_v3 = 5
calculateEvals_v3 = function(D){
  calculateEvals(D, N, true_beta_v3, mean_beta0, mean_beta1, 
                 var_beta0, var_beta1, var_e, f0, f1, true_type_v3, type, 
                 var_margy0 = NULL, var_margy1 = NULL, p, numSims)
}
k1_evals_v3 = calculateEvals_v3(one_at_a_time_k1)
k4_evals_v3 = calculateEvals_v3(one_at_a_time_k4)
k50_evals_v3 = calculateEvals_v3(one_at_a_time_k50)
fast_evals_v3 = calculateEvals_v3(fast_S5)
dopt_evals_v3 = calculateEvals_v3(doptimal)
space_evals_v3 = calculateEvals_v3(space_filling)
random_evals_v3 = calculateEvals_v3(random_design)

expMSEs_v3 = data.frame(cbind(k1_evals_v3$expMSEs, 
                              k4_evals_v3$expMSEs, 
                              k50_evals_v3$expMSEs, 
                              fast_evals_v3$expMSEs, 
                              dopt_evals_v3$expMSEs, 
                              space_evals_v3$expMSEs, 
                              random_evals_v3$expMSEs))
colnames(expMSEs_v3) = c("1atT,k=1","1atT,k=4","1atT,k=50","Fast","DOptimal","Space","Random")
rownames(expMSEs_v3) = c("E[MSE|H0]", "E[MSE|H1]")
round(expMSEs_v3, 5)





