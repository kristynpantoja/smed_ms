library(mvtnorm)

# --- sources to generate MEDs --- #
home = "/Users/kristyn/Documents/research/smed_ms"
functions_home = paste(home, "/functions", sep="")
source(paste(functions_home, "/wasserstein_distance.R", sep = ""))
source(paste(functions_home, "/charge_function_q.R", sep = ""))
source(paste(functions_home, "/variance_marginal_y.R", sep = ""))
source(paste(functions_home, "/variance_marginal_y.R", sep = ""))
source(paste(functions_home, "/generate_MED_oneatatime.R", sep = ""))
source(paste(functions_home, "/generate_MED_fast.R", sep = ""))

# --- sources to designs : MSE(Bn), E[P(H1|Y,D)] --- #
source(paste(functions_home, "/construct_design_matrix.R", sep = ""))
source(paste(functions_home, "/simulate_y.R", sep = ""))
source(paste(functions_home, "/postprob_hypotheses.R", sep = ""))
source(paste(functions_home, "/postmean_mse_mc.R", sep = ""))
source(paste(functions_home, "/posterior_variance.R", sep = ""))
source(paste(functions_home, "/postmean_mse_closedform.R", sep = ""))
source(paste(functions_home, "/plot_EPH1.R", sep = ""))

load(paste(home, "/designs/designs2d_0905.RData", sep = ""))
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
#space_filling = cbind(rep(axis1, each = 10), rep(axis2, times = 10))
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

# calculateEvals = function(D, N, true_beta, mean_beta0, mean_beta1, var_beta0, var_beta1,
#                           var_e, f0, f1, true_type, type, var_margy0 = NULL, var_margy1 = NULL, 
#                           p, numSims){
#   
#   # posterior variance of beta for the alternative hypothesis
#   v = postvar(D, N, var_e, var_beta1, type[2])
#   # Total Potential Energy criterion
#   TPE = totalPE_2d(D, N, mean_beta0, mean_beta1, var_beta0, var_beta1, var_e, 
#                    f0, f1, type, var_margy0, var_margy1, p)
#   # Fast Algorithm criterion
#   cFast = crit_fast_2d(D, N, mean_beta0, mean_beta1, var_beta0, var_beta1, var_e, 
#                        f0, f1, type, var_margy0, var_margy1, p)
#   # One-at-a-Time Algorithm criterion
#   c1atT = crit_1atatime_2d(D, N, k = 4, mean_beta0, mean_beta1, var_beta0, var_beta1, 
#                            var_e, f0, f1, type, var_margy0, var_margy1, p)
#   
#   # sumstats
#   means = apply(D, 2, mean)
#   sds = apply(D, 2, sd)
#   
#   # expected posterior probabilities of hypotheses and bayes factors
#   exppostprobs = calcExpPostProbH_2d(D, N, true_beta, mean_beta0, mean_beta1, 
#                                      var_beta0, var_beta1, var_e, numSims, true_type, type)
#   exppostmeanMSE = calcExpPostMeanMSE(D, N, true_beta, mean_beta0, mean_beta1, 
#                                       var_beta0, var_beta1, var_e,
#                                       numSims, true_type, type)
#   
#   # export
#   evals = list("v" = v, "TPE" = TPE, "cFast" = cFast, "c1atT" = c1atT, 
#                "mean" = means, "sd" = sds, "exppostprobs" = exppostprobs, "expMSEs" = exppostmeanMSE)
#   return(evals)
# }

# calculateEvals_v1 = function(D){
#   calculateEvals(D, N, true_beta_v1, mean_beta0, mean_beta1, 
#                  var_beta0, var_beta1, var_e, f0, f1, true_type_v1, type, 
#                  var_margy0 = NULL, var_margy1 = NULL, p, numSims)
# }
# k1_evals_v1 = calculateEvals_v1(one_at_a_time_k1)
# k4_evals_v1 = calculateEvals_v1(one_at_a_time_k4)
# k50_evals_v1 = calculateEvals_v1(one_at_a_time_k50)
# fast_evals_v1 = calculateEvals_v1(fast_S5)
# dopt_evals_v1 = calculateEvals_v1(doptimal)
# space_evals_v1 = calculateEvals_v1(space_filling)
# random_evals_v1 = calculateEvals_v1(random_design)
# 
# expMSEs_v1 = data.frame(cbind(k1_evals_v1$expMSEs, 
#                                    k4_evals_v1$expMSEs, 
#                                    k50_evals_v1$expMSEs, 
#                                    fast_evals_v1$expMSEs, 
#                                    dopt_evals_v1$expMSEs, 
#                                    space_evals_v1$expMSEs, 
#                                    random_evals_v1$expMSEs))
# colnames(expMSEs_v1) = c("1atT,k=1","1atT,k=4","1atT,k=50","Fast","DOptimal","Space","Random")
# rownames(expMSEs_v1) = c("MSEH0", "expEmpiricalMSEH1")
# round(expMSEs_v1, 5)



#########################################s
#########################################
#########################################
# TESTING MSE OF POSTERIOR MEAN OF BETA #
#########################################
#########################################
#########################################


# our functions:

#################################
### EMPIRICAL ESTIMATE OF MSE ###
#################################
getPostMean = function(y, D, N, beta_prior_mean, beta_prior_var, var_e, 
                       hypothesis_model_type, diagPrior = TRUE){
  X = constructDesignX(D, N, hypothesis_model_type)
  D_postvar = postvar(D, N, var_e, beta_prior_var, hypothesis_model_type, diagPrior)
  D_postmean = (1 / var_e) * D_postvar %*% (t(X) %*% y + var_e * solve(beta_prior_var, beta_prior_mean))
  return(D_postmean)
}


getDevianceSq = function(postmean_beta, true_beta){
  return((postmean_beta - true_beta)^2)
}

calcExpPostMeanMSE = function(D, N, true_beta, beta_prior_mean, beta_prior_var, var_e, type, numSims, diagPrior = TRUE, seed = 123){
  set.seed(seed)
  Ysims = simulateY(D, N, true_beta, var_e, numSims, type, seed)
  # calculating posterior mean
  post_means = apply(Ysims, 2, FUN = function(y) getPostMean(y, D, N, beta_prior_mean, beta_prior_var, 
                                                             var_e, type, diagPrior))
  # calculate squared deviances for each parameter Bni in Bn = (Bn1, ..., Bnp) from the true Bi in B = (B1, ..., Bp)
  empMSEs = apply(post_means, 2, FUN = function(x) getDevianceSq(x, true_beta))
  # get the mean squared deviances
  expEmpiricalMSE = apply(empMSEs, 1, mean)
  
  # also look at the variance of posterior mean...
  # var_postmeans = apply(post_means, 1, var)
  # and their means, to see if they're centered... i.e. if i calculated Bn correctly
  # avg_postmeans = apply(post_means, 1, mean)
  # return(list("expEmpiricalMSE" = expEmpiricalMSE, "var_postmeans" = var_postmeans, "avg_postmeans" = avg_postmeans))
  return("expEmpiricalMSE" = expEmpiricalMSE)
}

# calcExpPostMeanMSE_old = function(D, N, true_beta, beta_prior_mean0, beta_prior_mean1, 
#                               beta_prior_var0, beta_prior_var1, var_e,
#                               numSims = 100, true_model_type = NULL, H01_model_types = NULL, 
#                               seed = 123, diagPrior = TRUE){
#   set.seed(seed)
#   Ysims = simulateY(D, N, true_beta, var_e, numSims, true_model_type, seed)
#   # calculating posterior means from each hypothesis' prior on beta
#   if(true_model_type == 5 & H01_model_types[1] == 4 & H01_model_types[2] == 5){
#     # posterior mean, given H0 prior on beta
#     # we assume linear model to calculate MSE, since this model doesn't give estimates of quadratic terms
#     post_meansH0 = apply(Ysims, 2, FUN = function(y) getPostMean(y, D, N, beta_prior_mean0, beta_prior_var0, 
#                                                                  var_e, type[1], diagPrior))
#     empMSEsH0 = apply(post_meansH0, 2, FUN = function(x) getDevianceSq(x, true_beta[c(1, 2, 4)])) # Only MSE if assume linear model for estimation
#     # expEmpiricalMSEH0 = mean(empMSEsH0)
#     expEmpiricalMSEH0 = apply(empMSEsH0, 1, mean)
#     # posterior mean, given H1 prior on beta
#     post_meansH1 = apply(Ysims, 2, FUN = function(y) getPostMean(y, D, N, beta_prior_mean1, beta_prior_var1, 
#                                                                  var_e, type[2], diagPrior))
#     empMSEsH1 = apply(post_meansH1, 2, FUN = function(x) getDevianceSq(x, true_beta))
#     # expEmpiricalMSEH1 = mean(empMSEsH1)
#     expEmpiricalMSEH1 = apply(empMSEsH1, 1, mean)
#   } else if(true_model_type == 4 & H01_model_types[1] == 4 & H01_model_types[2] == 5){
#     # posterior mean, given H0 prior on beta
#     post_meansH0 = apply(Ysims, 2, FUN = function(y) getPostMean(y, D, N, beta_prior_mean0, beta_prior_var0, 
#                                                                  var_e, type[1], diagPrior))
#     # get squared deviances for each parameter Bi
#     empMSEsH0 = apply(post_meansH0, 2, FUN = function(x) getDevianceSq(x, true_beta)) # Only MSE if assume linear model for estimation
#     expEmpiricalMSEH0 = mean(empMSEsH0) # average squared deviances of each Bi
#     # posterior mean, given H1 prior on beta
#     # post_meansH1 = apply(Ysims, 2, FUN = function(y) getPostMean(y, D, N, beta_prior_mean1[c(1, 2, 4)], 
#     #                                                              diag(diag(beta_prior_var1)[c(1, 2, 4)]), 
#     #                                                              var_e, type[1], diagPrior))
#     # we ignore estimates for quadratic terms
#     post_meansH1 = apply(Ysims, 2, FUN = function(y) getPostMean(y, D, N, beta_prior_mean1, beta_prior_var1, 
#                                                                  var_e, type[2], diagPrior))
#     post_meansH1 = post_meansH1[c(1, 2, 4), ]
#     empMSEsH1 = apply(post_meansH1, 2, FUN = function(x) getDevianceSq(x, true_beta))
#     expEmpiricalMSEH1 = mean(empMSEsH1)
#   } else{
#     expEmpiricalMSEH0 = rep(NA, length(beta_prior_mean0))
#     expEmpiricalMSEH1 = rep(NA, length(beta_prior_mean1))
#   }
#   return(list("expEmpiricalMSEH0" = expEmpiricalMSEH0, "expEmpiricalMSEH1" = expEmpiricalMSEH1))
# }

####################################
### ANALYTICAL / CLOSED FORM MSE ###
####################################
getClosedMSE = function(D, N, true_beta, beta_prior_mean, beta_prior_var, var_e, type, diagPrior = TRUE){
  X = constructDesignX(D, N, type) # design matrix, depends on the model (type)
  Sigma_B = postvar(D, N, var_e, beta_prior_var, type, diagPrior) # posterior variance
  # 1. calculate the variance of the posterior mean
  XtX = crossprod(X)
  var_postmean_term1 = (1/var_e) * Sigma_B %*% XtX %*% Sigma_B
  # var_postmean_term2 = (1/var_e)^2 * Sigma_B %*% XtX %*% beta_prior_var %*% XtX %*% Sigma_B
  var_postmean_term2 = 0
  # grab the variance terms in the variance-covariance matrix
  var_postmean = diag(var_postmean_term1 + var_postmean_term2)
  # 2. calculate expectation of posterior mean
  # expect_postmean = (1/var_e) * Sigma_B %*% (XtX + var_e * solve(beta_prior_var)) %*% matrix(beta_prior_mean)
  expect_postmean = (1/var_e) * Sigma_B %*% XtX %*% matrix(true_beta) + Sigma_B %*% solve(beta_prior_var) %*% matrix(beta_prior_mean)
  # 3. MSE is a function of variance and expectation of posterior mean, as well as the true mean
  biassq_postmean = expect_postmean^2 - 2 * true_beta * expect_postmean + true_beta^2
  MSE_postmean = var_postmean + biassq_postmean
  return(list("var_term" = var_postmean, "biassq_term" = biassq_postmean, "MSE_postmean" = as.vector(MSE_postmean)))
}

getClosedMSEold = function(D, N, true_beta, beta_prior_mean, beta_prior_var, var_e, type, diagPrior = TRUE){
  X = constructDesignX(D, N, type) # design matrix, depends on the model (type)
  Sigma_B = postvar(D, N, var_e, beta_prior_var, type, diagPrior) # posterior variance
  # 1. calculate the variance of the posterior mean
  XtX = crossprod(X)
  var_postmean_term1 = (1/var_e) * Sigma_B %*% XtX %*% Sigma_B
  var_postmean_term2 = (1/var_e)^2 * Sigma_B %*% XtX %*% beta_prior_var %*% XtX %*% Sigma_B
  # grab the variance terms in the variance-covariance matrix
  var_postmean = diag(var_postmean_term1 + var_postmean_term2)
  # 2. calculate expectation of posterior mean
  expect_postmean = (1/var_e) * Sigma_B %*% (XtX + var_e * solve(beta_prior_var)) %*% matrix(beta_prior_mean)
  # 3. MSE is a function of variance and expectation of posterior mean, as well as the true mean
  biassq_postmean = expect_postmean^2 - 2 * true_beta * expect_postmean + true_beta^2
  MSE_postmean = var_postmean + biassq_postmean
  return(list("var_term" = var_postmean, "biassq_term" = biassq_postmean, "MSE_postmean" = as.vector(MSE_postmean)))
}


# First, we look at the case where the bias^2 should be 0, i.e. prior on beta is the true beta for each of H0 and H1

calcEmpMSE_v1 = function(D){
  expEmpMSEH0 = calcExpPostMeanMSE(D, N, c(0,0,0), mean_beta0, var_beta0, var_e, type[1], numSims)
  expEmpMSEH1 = calcExpPostMeanMSE(D, N, c(0, 0.2, -0.2, 0.2, -0.2), mean_beta1, var_beta1, var_e, type[2], numSims)
  return(list("expEmpMSEH0" = expEmpMSEH0, "expEmpMSEH1" = expEmpMSEH1))
}
calcClosedMSE_v1 = function(D){
  # closedMSEH0 = getClosedMSE(D, N, true_beta_v1[c(1,2,4)], mean_beta0, var_beta0, var_e, type[1])
  closedMSEH0 = getClosedMSE(D, N, c(0, 0, 0), mean_beta0, var_beta0, var_e, type[1])
  closedMSEH1 = getClosedMSE(D, N, c(0, 0.2, -0.2, 0.2, -0.2), mean_beta1, var_beta1, var_e, type[2])
  return(list("closedMSEH0" = closedMSEH0, "closedMSEH1" = closedMSEH1))
}

set.seed(123)
numSims = 1000
# Here, I've calculated MSEs for H0 and H1 with the true beta matching up with each of their priors
space_EmpMSE_v1 = calcEmpMSE_v1(space_filling)
space_closedMSE_v1 = calcClosedMSE_v1(space_filling)
# Some observations:
# the bias is appropriately 0 for both examples (H0 and H1)
# the variances don't match, which is why they aren't equal
# in empMSE, the variance = postmean, but also in closedMSE, the variance = postmean
#   so the question is, which one is correct and which one is not?
# Since empirical MSE should approach the true values, I am tempted to say that the closed MSE is wrong 
#   in its variance calculation.
#   But where? if bias calculation is correct, where could the variance calculation have been incorrect?
set.seed(123)
numSims = 5000
# Here, I've calculated MSEs for H0 and H1 with the true beta matching up with each of their priors
space_EmpMSE_v1_ns5k = calcEmpMSE_v1(space_filling)

# evaluate : if they look similar
# H0
space_closedMSE_v1$closedMSEH0$MSE_postmean
space_EmpMSE_v1_ns5k$expEmpMSEH0
# H1
space_closedMSE_v1$closedMSEH1$MSE_postmean
space_EmpMSE_v1_ns5k$expEmpMSEH1

# evaluate : if squared differences go to 0
max_index = 50
sqdiffH0_v1 = matrix(NA, 3, max_index)
sqdiffH1_v1 = matrix(NA, 5, max_index)
for(i in 1:max_index){
  numSims = 100*i
  space_EmpMSE_v1 = calcEmpMSE_v1(space_filling)
  space_closedMSE_v1 = calcClosedMSE_v1(space_filling)
  sqdiffH0_v1[,i] = (space_EmpMSE_v1$expEmpMSEH0 - space_closedMSE_v1$closedMSEH0$MSE_postmean)^2
  sqdiffH1_v1[,i] = (space_EmpMSE_v1$expEmpMSEH1 - space_closedMSE_v1$closedMSEH1$MSE_postmean)^2
}

beta_index = 2 #1:5
plot(x = 100*(1:max_index), y = sqdiffH0_v1[beta_index,], col = 2)
plot(x = 100*(1:max_index), y = sqdiffH1_v1[beta_index,], col = 4)

sumsqdiffH0_v1 = apply(sqdiffH0_v1, 2, sum)
sumsqdiffH1_v1 = apply(sqdiffH1_v1, 2, sum)
plot(x = 100*(1:max_index), y = sumsqdiffH0_v1, col = 2)
plot(x = 100*(1:max_index), sumsqdiffH1_v1, col = 4)























# What about when bias isn't 0?

true_beta_v2 = c(0, 0.7, -0.7, 0.7, -0.7)

numSims = 100
set.seed(123)


calcEmpMSE_v2 = function(D){
  expEmpMSEH0 = calcExpPostMeanMSE(D, N, true_beta_v2[c(1,2,4)], mean_beta0, var_beta0, var_e, type[1], numSims)
  expEmpMSEH1 = calcExpPostMeanMSE(D, N, true_beta_v2, mean_beta1, var_beta1, var_e, type[2], numSims)
  return(list("expEmpMSEH0" = expEmpMSEH0, "expEmpMSEH1" = expEmpMSEH1))
}
calcClosedMSE_v2 = function(D){
  closedMSEH0 = getClosedMSE(D, N, true_beta_v2[c(1,2,4)], mean_beta0, var_beta0, var_e, type[1])
  closedMSEH1 = getClosedMSE(D, N, true_beta_v2, mean_beta1, var_beta1, var_e, type[2])
  return(list("closedMSEH0" = closedMSEH0, "closedMSEH1" = closedMSEH1))
}

set.seed(123)
numSims = 1000
# Here, I've calculated MSEs for H0 and H1 with the true beta matching up with each of their priors
space_EmpMSE_v2 = calcEmpMSE_v2(space_filling)
space_closedMSE_v2 = calcClosedMSE_v2(space_filling)
# Some observations:
# the bias is appropriately 0 for both examples (H0 and H1)
# the variances don't match, which is why they aren't equal
# in empMSE, the variance = postmean, but also in closedMSE, the variance = postmean
#   so the question is, which one is correct and which one is not?
# Since empirical MSE should approach the true values, I am tempted to say that the closed MSE is wrong 
#   in its variance calculation.
#   But where? if bias calculation is correct, where could the variance calculation have been incorrect?
set.seed(123)
numSims = 5000
# Here, I've calculated MSEs for H0 and H1 with the true beta matching up with each of their priors
space_EmpMSE_v2_ns5k = calcEmpMSE_v2(space_filling)

# evaluate : if they look similar
# H0
space_closedMSE_v2$closedMSEH0$MSE_postmean
space_EmpMSE_v2_ns5k$expEmpMSEH0
# H1
space_closedMSE_v2$closedMSEH1$MSE_postmean
space_EmpMSE_v2_ns5k$expEmpMSEH1

# evaluate : if squared differences go to 0
max_index = 50
sqdiffH0_v2 = matrix(NA, 3, max_index)
sqdiffH1_v2 = matrix(NA, 5, max_index)
for(i in 1:max_index){
  numSims = 100*i
  space_EmpMSE_v2 = calcEmpMSE_v2(space_filling)
  space_closedMSE_v2 = calcClosedMSE_v2(space_filling)
  sqdiffH0_v2[,i] = (space_EmpMSE_v2$expEmpMSEH0 - space_closedMSE_v2$closedMSEH0$MSE_postmean)^2
  sqdiffH1_v2[,i] = (space_EmpMSE_v2$expEmpMSEH1 - space_closedMSE_v2$closedMSEH1$MSE_postmean)^2
}

beta_index = 5 #1:5
plot(x = 100*(1:max_index), y = sqdiffH0_v2[beta_index,], col = 2)
plot(x = 100*(1:max_index), y = sqdiffH1_v2[beta_index,], col = 4)

sumsqdiffH0_v2 = apply(sqdiffH0_v2, 2, sum)
sumsqdiffH1_v2 = apply(sqdiffH1_v2, 2, sum)
plot(x = 100*(1:max_index), y = sumsqdiffH0_v2, col = 2)
plot(x = 100*(1:max_index), sumsqdiffH1_v2, col = 4)







