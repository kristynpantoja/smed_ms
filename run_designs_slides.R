###################
# Loading Designs #
###################

load("designs_slides.RData")
# object loaded is list "all_designs"

# parameters/settings
mean_beta0 = c(0, 0, 0) # null model prior mean on beta
mean_beta1 = c(0, 0.2, -0.2, 0.2, -0.2) # alternative model prior mean on beta
var_mean0 = diag(c(0.005, 0.005, 0.005)) # null model prior variance on beta
var_mean1 = diag(c(0.005, 0.005, 0.005, 0.005, 0.005)) # alternative model prior variance on beta
var_e = 0.005 # variance on error
xmin = 0
xmax = 1
f0 = function(x) mean_beta0[1] + mean_beta0[2] * x[1] + mean_beta0[3] * x[2] # null
f1 = function(x) mean_beta1[1] + mean_beta1[2] * x[1] + mean_beta1[3] * x[1]^2 + mean_beta1[4] * x[2] + mean_beta1[3] * x[2]^2
N = 100
type = c(4, 5) # type (4,5) model is my indicator to the function of what kind of models are being compared so that marginal variance is correctly computed
p = 3


# for one-at-a-time algorithm, k = 1
numCandidates = 10e3
k = 1
one_at_a_time_k1 = all_designs$one_at_a_time_k1

# for one-at-a-time algorithm, k = 1
k = 4
one_at_a_time_k4 = all_designs$one_at_a_time_k4

# for one-at-a-time algorithm, k = 50
k = 50
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
source("med_ms_functions.R")
source("med_ms_fns_2d.R")

library(mvtnorm)

calculateEvals = function(D, N, mean_beta0, mean_beta1, var_mean0, var_mean1,
                          var_e, f0, f1, type, var_margy0 = NULL, var_margy1 = NULL, p){
  
  # posterior variance of beta for the alternative hypothesis
  v = postvar(D, N, var_e, var_mean1, type[2])
  # Total Potential Energy criterion
  TPE = totalPE_2d(D, N, mean_beta0, mean_beta1, var_mean0, var_mean1, var_e, 
                   f0, f1, type, var_margy0, var_margy1, p)
  # Fast Algorithm criterion
  cFast = crit_fast_2d(D, N, mean_beta0, mean_beta1, var_mean0, var_mean1, var_e, 
                       f0, f1, type, var_margy0, var_margy1, p)
  # One-at-a-Time Algorithm criterion
  c1atT = crit_1atatime_2d(D, N, k = 4, mean_beta0, mean_beta1, var_mean0, var_mean1, 
                           var_e, f0, f1, type, var_margy0, var_margy1, p)
  
  # sumstats
  means = apply(D, 2, mean)
  sds = apply(D, 2, sd)
  
  # expected posterior probabilities of hypotheses and bayes factors
  exppostprobs = suppressWarnings(calcExpPostProbH_2d(D, N, mean_beta0, mean_beta1, var_mean0, var_mean1, var_e,
                                                      numSims = 500, type = type, log_space = FALSE))
  exppostmeanMSE = calcExpPostMeanMSE(D, N, true_beta = mean_beta1, numSims = 1000, mean_beta1, var_e, var_mean1, type[2], diagPrior = TRUE)
  
  # export
  evals = list("v" = v, "TPE" = TPE, "cFast" = cFast, "c1atT" = c1atT, 
               "mean" = means, "sd" = sds, "exppostprobs" = exppostprobs, "expMSEs" = exppostmeanMSE)
  return(evals)
}

one_at_a_time_k1_evaluations = suppressWarnings(calculateEvals(one_at_a_time_k1, N, mean_beta0, mean_beta1, var_mean0, var_mean1, var_e, f0, f1, type, var_margy0 = NULL, var_margy1 = NULL, p))
# the warnings were just to remind myself that the hypotheses do not have models of the same dimensionality.
one_at_a_time_k1_evaluations$exppostprobs # these are the expected posterior probabilities of the models from simulated responses (simulated from H0 and H1)

# the rest of the design evaluations:
one_at_a_time_k4_evaluations = suppressWarnings(calculateEvals(one_at_a_time_k4, N, mean_beta0, mean_beta1, var_mean0, var_mean1, var_e, f0, f1, type, var_margy0 = NULL, var_margy1 = NULL, p))
one_at_a_time_k50_evaluations = suppressWarnings(calculateEvals(one_at_a_time_k50, N, mean_beta0, mean_beta1, var_mean0, var_mean1, var_e, f0, f1, type, var_margy0 = NULL, var_margy1 = NULL, p))
fast_S5_evaluations = suppressWarnings(calculateEvals(fast_S5, N, mean_beta0, mean_beta1, var_mean0, var_mean1, var_e, f0, f1, type, var_margy0 = NULL, var_margy1 = NULL, p))
doptimal_evaluations = suppressWarnings(calculateEvals(doptimal, N, mean_beta0, mean_beta1, var_mean0, var_mean1, var_e, f0, f1, type, var_margy0 = NULL, var_margy1 = NULL, p))
space_filling_evaluations = suppressWarnings(calculateEvals(space_filling, N, mean_beta0, mean_beta1, var_mean0, var_mean1, var_e, f0, f1, type, var_margy0 = NULL, var_margy1 = NULL, p))
random_design_evaluations = suppressWarnings(calculateEvals(random_design, N, mean_beta0, mean_beta1, var_mean0, var_mean1, var_e, f0, f1, type, var_margy0 = NULL, var_margy1 = NULL, p))


expectedPosteriorProbs = data.frame(cbind(one_at_a_time_k1_evaluations$exppostprobs, 
                                          one_at_a_time_k4_evaluations$exppostprobs, 
                                          one_at_a_time_k50_evaluations$exppostprobs, 
                                          fast_S5_evaluations$exppostprobs, 
                                          doptimal_evaluations$exppostprobs, 
                                          space_filling_evaluations$exppostprobs, 
                                          random_design_evaluations$exppostprobs))
colnames(expectedPosteriorProbs) = c("1atT,k=1","1atT,k=4","1atT,k=50","Fast","DOptimal","Space","Random")
rownames(expectedPosteriorProbs) = c("E[P(H0|Y,D)|H0,D]", "E[P(H1|Y,D)|H0,D]", "E[BF01 | H0,D]", 
                                     "E[P(H0|Y,D)|H1,D]", "E[P(H1|Y,D)|H1,D]", "E[BF01|H1,D]")
expectedPosteriorProbs


