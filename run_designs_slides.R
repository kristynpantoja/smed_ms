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
# one_at_a_time_k1_evaluations$exppostprobs # these are the expected posterior probabilities of the models from simulated responses (simulated from H0 and H1)

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
round(expectedPosteriorProbs,5)


# Trying to figure out why the MEDs do not have higher expected posterior probabilities than all other designs
# => why E[P(H0|Y,D)|H0,D], E[P(H1|Y,D)|H1,D] are not highest for MEDs.

# Dr. Jones: Check why the Bayesian probability is not higher for our methods than for the space filling design 
# e.g. check the mean density of the points in the space filling design and our designs (under H1). 
# If the points in our design have higher density then I think there must be a mistake in the computation of the posterior probability. 
# If the mean density is not higher in our design then there is something wrong with the criterion we are optimizing.  

# note: (I think he meant to say "If the points in our design DON'T have higher density")

# checking mean density: #


numSims = 100

# for data simulated from H1:
# fast MED
fast_simY_H1 = simulateY(fast_S5, N, mean_beta1, var_mean1, var_e, numSims, type = type[2])
fast_simEvidenceH1_YH1 = rep(NA, numSims)
for(j in 1:numSims){
  Y = fast_simY_H1[, j]
  fast_simEvidenceH1_YH1[j] = model_evidence(Y, fast_S5, N, mean_beta1, var_mean1, var_e, type = type[2])
}
# space-filling
space_simY_H1 = simulateY(space_filling, N, mean_beta1, var_mean1, var_e, numSims, type = type[2])
space_simEvidenceH1_YH1 = rep(NA, numSims)
for(j in 1:numSims){
  Y = space_simY_H1[, j]
  space_simEvidenceH1_YH1[j] = model_evidence(Y, space_filling, N, mean_beta1, var_mean1, var_e, type = type[2])
}
mean(fast_simEvidenceH1_YH1) > mean(space_simEvidenceH1_YH1) # fast MED gives higher mean density than space MED for true hypothesis (here, H1), as expected.



# is this true for H0-simulated data also? it should be.
# fast MED
fast_simY_H0 = simulateY(fast_S5, N, mean_beta0, var_mean0, var_e, numSims, type = type[1])
fast_simEvidenceH0_YH0 = rep(NA, numSims) # added for bf01 calc
fast_simEvidenceH1_YH0 = rep(NA, numSims)
for(j in 1:numSims){
  Y = fast_simY_H0[, j]
  fast_simEvidenceH0_YH0[j] = model_evidence(Y, fast_S5, N, mean_beta0, var_mean0, var_e, type = type[1])
  fast_simEvidenceH1_YH0[j] = model_evidence(Y, fast_S5, N, mean_beta1, var_mean1, var_e, type = type[2])
}
# space-filling
space_simY_H0 = simulateY(space_filling, N, mean_beta0, var_mean0, var_e, numSims, type = type[1])
space_simEvidenceH0_YH0 = rep(NA, numSims) # added for bf01 calc
space_simEvidenceH1_YH0 = rep(NA, numSims)
for(j in 1:numSims){
  Y = space_simY_H0[, j]
  space_simEvidenceH0_YH0[j] = model_evidence(Y, space_filling, N, mean_beta0, var_mean0, var_e, type = type[1])
  space_simEvidenceH1_YH0[j] = model_evidence(Y, space_filling, N, mean_beta1, var_mean1, var_e, type = type[2])
}

mean(fast_simEvidenceH1_YH0) > mean(space_simEvidenceH1_YH0) # again, fast MED gives higher mean density than space MED for true hypothesis (here, H0), as expected.

# and...
mean(fast_simEvidenceH1_YH0)/mean(fast_simEvidenceH1_YH1) > mean(space_simEvidenceH1_YH0)/mean(space_simEvidenceH1_YH1)

# but what about bf? evidence_H0/evidence_H1
mean(fast_simEvidenceH0_YH0)/mean(fast_simEvidenceH1_YH0) > mean(space_simEvidenceH0_YH0)/mean(space_simEvidenceH1_YH0)
# why is it false? if they had the same denominator, this would not be true. but evidence for H1 is higher for MED. why?
mean(fast_simEvidenceH1_YH0) > mean(space_simEvidenceH1_YH0) # well why is the evidence for H1 in MED greater than the evidence for H1 in space-filling design?
# this is the problem! why does it exist?
sum(round(fast_simEvidenceH1_YH0, 5) > round(space_simEvidenceH1_YH0, 5))/numSims # it's not from chance!
# and it's significant! look:
mean(fast_simEvidenceH1_YH0) / 1e54
mean(space_simEvidenceH1_YH0) / 1e54

# is evidence even higher for H0 than for H1 for a given design?
mean(fast_simEvidenceH0_YH0) > mean(fast_simEvidenceH1_YH0)
mean(space_simEvidenceH0_YH0) > mean(space_simEvidenceH1_YH0)
# okay well at least there's that.

# gotta look at model evidence calculation:

j = 10
model_evidence(fast_simY_H0[, j], fast_S5, N, mean_beta1, var_mean1, var_e, type = type[2])
model_evidence(space_simY_H0[, j], space_filling, N, mean_beta1, var_mean1, var_e, type = type[2])

model_evidence = function(Y, D, N, mean_beta, var_mean, var_e, type){
  # Y is a vector
  # X is a matrix
  # var_mean is a matrix
  # var_e is a scalar
  if(N != length(Y)) stop("N is not the same length as Y")
  X = constructDesignX(D, N, type)
  marginaly_mean = X %*% mean_beta
  if(dim(X)[1] > 1){
    marginaly_var = diag(rep(var_e, N)) + (X %*% var_mean %*% t(X))
  } else{
    marginaly_var = var_e + (X %*% var_mean %*% t(X))
  }
  return(dmvnorm(Y, mean = marginaly_mean, sigma = marginaly_var, log = FALSE))
}






















































































































































































































# Then what's the issue? Look into the rest of the expected posterior probability calculation...

# for H0:
# fast
fast_simY_H0 = simulateY(fast_S5, N, mean_beta0, var_mean0, var_e, numSims, type = type[1])
fast_simEvidenceH0_YH0 = rep(NA, numSims)
fast_simEvidenceH1_YH0 = rep(NA, numSims)
fast_simPostH0 = rep(NA, numSims)
fast_simPostH1 = rep(NA, numSims)
fast_simBF01 = rep(NA, numSims)
for(j in 1:numSims){
  Y = fast_simY_H0[, j]
  # get model evidences
  fast_simEvidenceH0_YH0[j] = model_evidence(Y, fast_S5, N, mean_beta0, var_mean0, var_e, type = type[1])
  fast_simEvidenceH1_YH0[j] = model_evidence(Y, fast_S5, N, mean_beta1, var_mean1, var_e, type = type[2])
  # calculate posterior probabilities of models
  fast_simPostH0[j] = fast_simEvidenceH0_YH0[j] / (fast_simEvidenceH0_YH0[j] + fast_simEvidenceH1_YH0[j])
  fast_simPostH1[j] = fast_simEvidenceH1_YH0[j] / (fast_simEvidenceH0_YH0[j] + fast_simEvidenceH1_YH0[j])
  # calculate bayes factor
  fast_simBF01[j] = fast_simPostH0[j] / fast_simPostH1[j]
}
fast_expected_postH0_YH0 = mean(fast_simPostH0)
fast_expected_postH1_YH0 = mean(fast_simPostH1)
fast_expected_BF01_YH0 = mean(fast_simBF01)
# space_filling
space_simY_H0 = simulateY(space_filling, N, mean_beta0, var_mean0, var_e, numSims, type = type[1])
space_simEvidenceH0_YH0 = rep(NA, numSims)
space_simEvidenceH1_YH0 = rep(NA, numSims)
space_simPostH0 = rep(NA, numSims)
space_simPostH1 = rep(NA, numSims)
space_simBF01 = rep(NA, numSims)
for(j in 1:numSims){
  Y = space_simY_H0[, j]
  # get model evidences
  space_simEvidenceH0_YH0[j] = model_evidence(Y, space_filling, N, mean_beta0, var_mean0, var_e, type = type[1])
  space_simEvidenceH1_YH0[j] = model_evidence(Y, space_filling, N, mean_beta1, var_mean1, var_e, type = type[2])
  # calculate posterior probabilities of models
  space_simPostH0[j] = space_simEvidenceH0_YH0[j] / (space_simEvidenceH0_YH0[j] + space_simEvidenceH1_YH0[j])
  space_simPostH1[j] = space_simEvidenceH1_YH0[j] / (space_simEvidenceH0_YH0[j] + space_simEvidenceH1_YH0[j])
  # calculate bayes factor
  space_simBF01[j] = space_simPostH0[j] / space_simPostH1[j]
}
space_expected_postH0_YH0 = mean(space_simPostH0)
space_expected_postH1_YH0 = mean(space_simPostH1)
space_expected_BF01_YH0 = mean(space_simBF01)


fast_expected_postH0_YH0 > space_expected_postH0_YH0 # shouldn't be false :(

mean(fast_simEvidenceH1_YH0)/mean(fast_simEvidenceH1_YH0) > mean(space_simEvidenceH1_YH0)/mean(space_simEvidenceH1_YH0)

