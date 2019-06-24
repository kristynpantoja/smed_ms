library(expm)

# --- Kristyn --- #
#home = "/Users/kristyn/Documents/research/smed_ms/"
#output_home = paste(home,"output/",sep="")

# --- Cluster --- #
home = "/scratch/user/kristynp/med_ms/"
output_home = paste(home,"output/",sep="")

med_fns = paste(home,"med_ms_functions.R",sep="")
source(med_fns)



###########################################
###########################################
# Simple Linear Regression: Unknown Slope #
###########################################
###########################################


########################
# One-at-a-Time, k = 1 #
########################

mean_beta0 = 1 # slope of null model
mean_beta1 = 1 / 2 # slope of alternative model
var_mean0 = 0.005; var_mean1 = var_mean0; # variance on beta
var_e = 0.025 # variance on error
xmin = 0
xmax = 1
f0 = function(x) mean_beta0 * x # null regression model
f1 = function(x) mean_beta1 * x # alternative regression model
N = 67
type = c(1, 1)
p = 1
# for one-at-a-time algorithm:
numCandidates = 10^4 # suggested 10^5
k = 1
# One-at-a-Time Algorithm
D_oaat_k1_slope = MED_ms(mean_beta0 = mean_beta0, mean_beta1 = mean_beta1, 
                         var_mean0 = var_mean0, var_mean1 = var_mean1, var_e = var_e, 
                         f0 = f0, f1 = f1, type = type, var_margy0 = NULL, var_margy1 = NULL, 
                         N = N, numCandidates = numCandidates, k = k, p = p, xmin = xmin, xmax = xmax, 
                         genCandidates = 1, initialpt = 1)
saveRDS(D_oaat_k1_slope, paste(output_home,"D_oaat_k1_slope.rds",sep=""))

########################
# One-at-a-Time, k = 4 #
########################

mean_beta0 = 1 # slope of null model
mean_beta1 = 1 / 2 # slope of alternative model
var_mean0 = 0.005; var_mean1 = var_mean0; # variance on beta
var_e = 0.025 # variance on error
xmin = 0
xmax = 1
f0 = function(x) mean_beta0 * x # null regression model
f1 = function(x) mean_beta1 * x # alternative regression model
N = 67
type = c(1, 1)
p = 1
# for one-at-a-time algorithm:
numCandidates = 10^4 # suggested 10^5
k = 4
# One-at-a-Time Algorithm
D_oaat_k4_slope = MED_ms(mean_beta0 = mean_beta0, mean_beta1 = mean_beta1, 
                         var_mean0 = var_mean0, var_mean1 = var_mean1, var_e = var_e, 
                         f0 = f0, f1 = f1, type = type, var_margy0 = NULL, var_margy1 = NULL, 
                         N = N, numCandidates = numCandidates, k = k, p = p, xmin = xmin, xmax = xmax, 
                         genCandidates = 1, initialpt = 1)
saveRDS(D_oaat_k4_slope, paste(output_home,"D_oaat_k4_slope.rds",sep=""))

################
# Fast, S = 20 #
################

mean_beta0 = 1 # slope of null model
mean_beta1 = 1 / 2 # slope of alternative model
var_mean0 = 0.005; var_mean1 = var_mean0; # variance on beta
var_e = 0.025 # variance on error
xmin = 0
xmax = 1
f0 = function(x) mean_beta0 * x # null regression model
f1 = function(x) mean_beta1 * x # alternative regression model
N = 67
type = c(1, 1)
p = 1
# for fast algorithm:
S = 20 # ceiling(4* sqrt(p))
numParameters = 1 # number of parameters (just slope!)
p = numParameters
## Fast Algorithm
D_fast_S20_slope = MED_ms_fast(mean_beta0 = mean_beta0, mean_beta1 = mean_beta1, 
                               var_mean0 = var_mean0, var_mean1 = var_mean1, var_e = var_e, 
                               f0 = f0, f1 = f1, type = type, var_margy0 = NULL, var_margy1 = NULL, 
                               N = N, K = S, p = p, xmin = xmin, xmax = xmax, 
                               genCandidates = 1, initialpt = 1)
saveRDS(D_fast_S20_slope, paste(output_home,"D_fast_S20_slope.rds",sep=""))



#######################################################
#######################################################
# Simple Linear Regression: Unknown Intercept & Slope #
#######################################################
#######################################################


########################
# One-at-a-Time, k = 1 #
########################

mean_beta0 = c(0, 1) # slope of null model
mean_beta1 = c(0, 1/2) # slope of alternative model
var_mean0 = diag(c(0.005, 0.005)); var_mean1 = var_mean0 # variance on beta
var_e = 0.025 # variance on error
xmin = 0
xmax = 1
f0 = function(x) mean_beta0[1] + mean_beta0[2] * x # null regression model
f1 = function(x) mean_beta1[1] + mean_beta1[2] * x # alternative regression model
N = 67
type = c(2, 2)
p = 2
# for one-at-a-time algorithm:
numCandidates = 10^4 # suggested 10^5
k = 1
# One-at-a-Time Algorithm
D_oaat_k1_lin = MED_ms(mean_beta0 = mean_beta0, mean_beta1 = mean_beta1, 
                       var_mean0 = var_mean0, var_mean1 = var_mean1, var_e = var_e, 
                       f0 = f0, f1 = f1, type = type, var_margy0 = NULL, var_margy1 = NULL, 
                       N = N, numCandidates = numCandidates, p = p, k = k, xmin = xmin, xmax = xmax, 
                       genCandidates = 1, initialpt = 1)
saveRDS(D_oaat_k1_lin, paste(output_home,"D_oaat_k1_lin.rds",sep=""))

########################
# One-at-a-Time, k = 4 #
########################

mean_beta0 = c(0, 1) # slope of null model
mean_beta1 = c(0, 1/2) # slope of alternative model
var_mean0 = diag(c(0.005, 0.005)); var_mean1 = var_mean0 # variance on beta
var_e = 0.025 # variance on error
xmin = 0
xmax = 1
f0 = function(x) mean_beta0[1] + mean_beta0[2] * x # null regression model
f1 = function(x) mean_beta1[1] + mean_beta1[2] * x # alternative regression model
N = 67
type = c(2, 2)
p = 2
# for one-at-a-time algorithm:
numCandidates = 10^4 # suggested 10^5
k = 4
# One-at-a-Time Algorithm
D_oaat_k4_lin = MED_ms(mean_beta0 = mean_beta0, mean_beta1 = mean_beta1, 
                       var_mean0 = var_mean0, var_mean1 = var_mean1, var_e = var_e, 
                       f0 = f0, f1 = f1, type = type, var_margy0 = NULL, var_margy1 = NULL, 
                       N = N, numCandidates = numCandidates, k = k, p = p, xmin = xmin, xmax = xmax, 
                       genCandidates = 1, initialpt = 1)
saveRDS(D_oaat_k4_lin, paste(output_home,"D_oaat_k4_lin.rds",sep=""))

################
# Fast, S = 20 #
################

mean_beta0 = c(0, 1) # slope of null model
mean_beta1 = c(0, 1/2) # slope of alternative model
var_mean0 = diag(c(0.005, 0.005)); var_mean1 = var_mean0 # variance on beta
var_e = 0.025 # variance on error
xmin = 0
xmax = 1 
f0 = function(x) mean_beta0[1] + mean_beta0[2] * x # null regression model
f1 = function(x) mean_beta1[1] + mean_beta1[2] * x # alternative regression model
N = 67
type = c(2, 2)
p = 2
# for fast algorithm:
S = 20 # ceiling(4* sqrt(p))
numParameters = 1 # number of parameters (just slope!)
p = numParameters
## Fast Algorithm
D_fast_S20_lin = MED_ms_fast(mean_beta0 = mean_beta0, mean_beta1 = mean_beta1, 
                             var_mean0 = var_mean0, var_mean1 = var_mean1, var_e = var_e, 
                             f0 = f0, f1 = f1, type = type, var_margy0 = NULL, var_margy1 = NULL, 
                             N = N, K = S, p = p, xmin = xmin, xmax = xmax, 
                             genCandidates = 1, initialpt = 1)
saveRDS(D_fast_S20_lin, paste(output_home,"D_fast_S20_lin.rds",sep=""))


#######################
#######################
# Linear vs Quadratic #
#######################
#######################


########################
# One-at-a-Time, k = 1 #
########################

mean_beta0 = c(1, 0) # slope of null model
mean_beta1 = c(1, 4, -4) # slope of alternative model
var_mean0 = diag(c(0.005, 0.005)) # variance on beta0
var_mean1 = diag(c(0.005, 0.005, 0.005)) # variance on beta1
var_e = 0.025 # variance on error
xmin = 0
xmax = 1
f0 = function(x) mean_beta0[1] + mean_beta0[2] * x # null regression model
f1 = function(x) mean_beta1[1] + mean_beta1[2] * x + mean_beta1[3] * x^2 # alternative regression model
N = 67
type = c(2, 3)
p = 3
# for one-at-a-time algorithm:
numCandidates = 10^4 # suggested 10^5
k = 1
# One-at-a-Time Algorithm
D_oaat_k1_lq = MED_ms(mean_beta0 = mean_beta0, mean_beta1 = mean_beta1, 
                      var_mean0 = var_mean0, var_mean1 = var_mean1, var_e = var_e, 
                      f0 = f0, f1 = f1, type = type, var_margy0 = NULL, var_margy1 = NULL, 
                      N = N, numCandidates = numCandidates, p = p, k = k, xmin = xmin, xmax = xmax, 
                      genCandidates = 1, initialpt = 1)
saveRDS(D_oaat_k1_lq, paste(output_home,"D_oaat_k1_lq.rds",sep=""))

########################
# One-at-a-Time, k = 4 #
########################

mean_beta0 = c(1, 0) # slope of null model
mean_beta1 = c(1, 4, -4) # slope of alternative model
var_mean0 = diag(c(0.005, 0.005)) # variance on beta0
var_mean1 = diag(c(0.005, 0.005, 0.005)) # variance on beta1
var_e = 0.025 # variance on error
xmin = 0
xmax = 1
f0 = function(x) mean_beta0[1] + mean_beta0[2] * x # null regression model
f1 = function(x) mean_beta1[1] + mean_beta1[2] * x + mean_beta1[3] * x^2 # alternative regression model
N = 67
type = c(2, 3)
p = 3
# for one-at-a-time algorithm:
numCandidates = 10^4 # suggested 10^5
k = 4
# One-at-a-Time Algorithm
D_oaat_k4_lq = MED_ms(mean_beta0 = mean_beta0, mean_beta1 = mean_beta1, 
                      var_mean0 = var_mean0, var_mean1 = var_mean1, var_e = var_e, 
                      f0 = f0, f1 = f1, type = type, var_margy0 = NULL, var_margy1 = NULL, 
                      N = N, numCandidates = numCandidates, k = k, p = p, xmin = xmin, xmax = xmax, 
                      genCandidates = 1, initialpt = 1)
saveRDS(D_oaat_k4_lq, paste(output_home,"D_oaat_k4_lq.rds",sep=""))
################
# Fast, S = 20 #
################

mean_beta0 = c(1, 0) # slope of null model
mean_beta1 = c(1, 4, -4) # slope of alternative model
var_mean0 = diag(c(0.005, 0.005)) # variance on beta0
var_mean1 = diag(c(0.005, 0.005, 0.005)) # variance on beta1
var_e = 0.025 # variance on error
xmin = 0
xmax = 1
f0 = function(x) mean_beta0[1] + mean_beta0[2] * x # null regression model
f1 = function(x) mean_beta1[1] + mean_beta1[2] * x + mean_beta1[3] * x^2 # alternative regression model
N = 67
type = c(2, 3)
p = 3
# for fast algorithm:
S = 20 # ceiling(4* sqrt(p))
numParameters = 1 # number of parameters (just slope!)
p = numParameters
## Fast Algorithm
D_fast_S20_lq = MED_ms_fast(mean_beta0 = mean_beta0, mean_beta1 = mean_beta1, 
                            var_mean0 = var_mean0, var_mean1 = var_mean1, var_e = var_e, 
                            f0 = f0, f1 = f1, type = type, var_margy0 = NULL, var_margy1 = NULL, 
                            N = N, K = S, p = p, xmin = xmin, xmax = xmax, 
                            genCandidates = 1, initialpt = 1)
saveRDS(D_fast_S20_lq, paste(output_home,"D_fast_S20_lq.rds",sep=""))


