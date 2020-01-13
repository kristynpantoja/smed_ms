# Originally ran in run_designs/limiting_distribution/limiting_distribution4.R (bottom of script, with k = 2 though)

# --- Directory --- #

# Kristyn
# home = "/Users/kristyn/Documents/smed_ms/"

# Cluster
home = "/scratch/user/kristynp/smed_ms/"
output_home = paste(home,"run_designs/",sep="")

# --- Sources/Libraries --- #
functions_home = paste(home, "/functions", sep="")
source(paste(functions_home, "/wasserstein_distance.R", sep = ""))
source(paste(functions_home, "/charge_function_q.R", sep = ""))
source(paste(functions_home, "/variance_marginal_y.R", sep = ""))
source(paste(functions_home, "/generate_MED_oneatatime.R", sep = ""))
source(paste(functions_home, "/generate_MED_fast.R", sep = ""))
library(expm)
library(matrixStats)

mean_beta0 = c(0, 1 / 2) # slope of null model
mean_beta1 = c(0, 1 / 4) # slope of alternative model
var_e = 0.1
var_mean_elt = 0.25
var_mean0 = diag(rep(var_mean_elt, 2)); var_mean1 = var_mean0 # variance on beta
xmin = 0
xmax = 1
f0 = function(x) mean_beta0[1] + mean_beta0[2] * x # null regression model
f1 = function(x) mean_beta1[1] + mean_beta1[2] * x # alternative regression model
N = 500
type = c(2, 2)
# for one-at-a-time algorithm:
numCandidates = 10^4
k = 4
p = 1
greedy_N500 = MED_ms_oneatatime(mean_beta0, mean_beta1, var_mean0, var_mean1, var_e,
                                f0, f1, type, N = N, numCandidates = numCandidates, k = k,
                                xmin = xmin, xmax = xmax, p = p, alpha = 1, buffer = 0,
                                genCandidates = 1, initialpt = 1)
saveRDS(greedy_N500, paste(output_home,"original_example.rds", sep = ""))
###
# f_dens = function(x, mean_beta0, mean_beta1, var_mean0, var_mean1, var_e, f0, f1, type){
#   if(length(type) != 2) stop("type should be vector with length == 2")
#   mu1 = f0(x) # mean of marginal dist of y | H0
#   mu2 = f1(x) # mean of marginal dist of y | H1
#   var1 = var_marginaly(x, var_mean0, var_e, type = type[1], var_margy0) # variance of marginal dist of y | H0
#   var2 = var_marginaly(x, var_mean1, var_e, type = type[2], var_margy1) # variance of marginal dist of y | H1
#   Wass_dist = Wasserstein_distance(mu1, mu2, var1, var2)
#   return(Wass_dist)
# }
# unnormalized_f_dens_curve = function(x) f_dens(x,mean_beta0, mean_beta1, var_mean0, var_mean1, var_e, f0, f1, type)
# integratefn_int_const = integrate(f = unnormalized_f_dens_curve, lower = xmin, upper = xmax)$value
# f_dens_curve2 = function(x) (1/integratefn_int_const) * unnormalized_f_dens_curve(x)
# ###
# hist(greedy_N500, probability = T, ylim = c(0, 2.5), breaks = 10)
# curve(f_dens_curve2, add = T)