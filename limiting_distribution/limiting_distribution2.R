###########################
# Making W a Distribution #
###########################

# --- Directory --- #
# Kristyn
home = "/Users/kristyn/Documents/research/smed_ms/"

# --- Sources/Libraries --- #
functions_home = paste(home, "/functions", sep="")
source(paste(functions_home, "/wasserstein_distance.R", sep = ""))
source(paste(functions_home, "/charge_function_q.R", sep = ""))
source(paste(functions_home, "/variance_marginal_y.R", sep = ""))
source(paste(functions_home, "/generate_MED_oneatatime.R", sep = ""))
source(paste(functions_home, "/generate_MED_fast.R", sep = ""))
library(expm)
library(matrixStats)

#######################################
# xmin = -1
#######################################

# Model 1 : intercept = 0
mean_beta0 = c(0, 1 / 2) # slope of null model
mean_beta1 = c(0, 1 / 4) # slope of alternative model
var_mean0 = diag(c(0.25, 0.25)); var_mean1 = var_mean0 # variance on beta
var_e = 0.1 # variance on error
xmin = 0
xmax = 1
f0 = function(x) mean_beta0[1] + mean_beta0[2] * x # null regression model
f1 = function(x) mean_beta1[1] + mean_beta1[2] * x # alternative regression model
N = 100
type = c(2, 2)
p = 2
# for fast algorithm:
S = 5
# mmed = MED_ms_oneatatime(mean_beta0, mean_beta1, var_mean0, var_mean1, var_e, 
#                          f0, f1, type, N = N, numCandidates = 10^3, k = 4, 
#                          xmin = -1, xmax = 1, p = 2, alpha = 1, buffer = 0, 
#                          genCandidates = 1, initialpt = 1)
# hist(mmed)
# mmed2 = MED_ms_oneatatime(mean_beta0, mean_beta1, var_mean0, var_mean1, var_e, 
#                          f0, f1, type, N = 300, numCandidates = 10^3, k = 4, 
#                          xmin = -1, xmax = 1, p = 2, alpha = 1, buffer = 0, 
#                          genCandidates = 1, initialpt = 1)
mmed = MED_ms_oneatatime(mean_beta0, mean_beta1, var_mean0, var_mean1, var_e, 
                          f0, f1, type, N = 100, numCandidates = 10e3, k = 4, 
                          xmin = xmin, xmax = xmax, p = 2, alpha = 1, buffer = 0, 
                          genCandidates = 1, initialpt = 1)

###
f_dens = function(x, mean_beta0, mean_beta1, var_mean0, var_mean1, var_e, f0, f1, type){
  if(length(type) != 2) stop("type should be vector with length == 2")
  mu1 = f0(x) # mean of marginal dist of y | H0
  mu2 = f1(x) # mean of marginal dist of y | H1
  var1 = var_marginaly(x, var_mean0, var_e, type = type[1], var_margy0) # variance of marginal dist of y | H0
  var2 = var_marginaly(x, var_mean1, var_e, type = type[2], var_margy1) # variance of marginal dist of y | H1
  Wass_dist = Wasserstein_distance(mu1, mu2, var1, var2)
  return(Wass_dist)
}

unnormalized_f_dens_curve = function(x) f_dens(x,mean_beta0, mean_beta1, var_mean0, var_mean1, var_e, f0, f1, type)
integratefn_int_const = integrate(f = unnormalized_f_dens_curve, lower = xmin, upper = xmax)$value


f_dens_curve2 = function(x) (1/integratefn_int_const) * unnormalized_f_dens_curve(x)



###
hist(mmed, probability = T, ylim = c(0, 2.5), breaks = 10)
curve(f_dens_curve2, add = T)

length(mmed)
saveRDS(mmed, "mmed_N100.rds")



#######################################
# xmin = -1
#######################################

# Model 1 : intercept = 0
mean_beta0 = c(0, 1 / 2) # slope of null model
mean_beta1 = c(0, 1 / 4) # slope of alternative model
var_mean0 = diag(c(0.005, 0.005)); var_mean1 = var_mean0 # variance on beta
var_e = 0.025 # variance on error
xmin = -1
xmax = 1
f0 = function(x) mean_beta0[1] + mean_beta0[2] * x # null regression model
f1 = function(x) mean_beta1[1] + mean_beta1[2] * x # alternative regression model
N = 100
type = c(2, 2)
p = 2
# for fast algorithm:
S = 5
# mmed = MED_ms_oneatatime(mean_beta0, mean_beta1, var_mean0, var_mean1, var_e, 
#                          f0, f1, type, N = N, numCandidates = 10^3, k = 4, 
#                          xmin = -1, xmax = 1, p = 2, alpha = 1, buffer = 0, 
#                          genCandidates = 1, initialpt = 1)
# hist(mmed)
# mmed2 = MED_ms_oneatatime(mean_beta0, mean_beta1, var_mean0, var_mean1, var_e, 
#                          f0, f1, type, N = 300, numCandidates = 10^3, k = 4, 
#                          xmin = -1, xmax = 1, p = 2, alpha = 1, buffer = 0, 
#                          genCandidates = 1, initialpt = 1)
mmed3 = MED_ms_oneatatime(mean_beta0, mean_beta1, var_mean0, var_mean1, var_e, 
                          f0, f1, type, N = 300, numCandidates = 5000, k = 4, 
                          xmin = -1, xmax = 1, p = 2, alpha = 1, buffer = 0, 
                          genCandidates = 1, initialpt = 1)

###
int_const = mean_beta0[1]^2 - mean_beta0[1] * mean_beta0[2] + (1/3) * mean_beta0[2]^2 + 
  mean_beta1[1]^2 - mean_beta1[1] * mean_beta1[2] + (1/3) * mean_beta1[2]^2 - 
  2 * mean_beta0[1] * mean_beta1[1] - mean_beta0[1] * mean_beta1[2] - mean_beta0[2] * mean_beta1[1] - 
  (2/3) * mean_beta0[2] * mean_beta1[2]
int_const = sqrt(int_const)

# check, where mean_beta0[1] and mean_beta1[1] are both 0
#int_const2 = (1/3) * (mean_beta0[2] - mean_beta1[2])^2
#int_const2 = sqrt(int_const2)
# it worked

###
f_dens = function(x, mean_beta0, mean_beta1, var_mean0, var_mean1, var_e, f0, f1, type){
  if(length(type) != 2) stop("type should be vector with length == 2")
  mu1 = f0(x) # mean of marginal dist of y | H0
  mu2 = f1(x) # mean of marginal dist of y | H1
  var1 = var_marginaly(x, var_mean0, var_e, type = type[1], var_margy0) # variance of marginal dist of y | H0
  var2 = var_marginaly(x, var_mean1, var_e, type = type[2], var_margy1) # variance of marginal dist of y | H1
  Wass_dist = Wasserstein_distance(mu1, mu2, var1, var2)
  return(Wass_dist)
}

unnormalized_f_dens_curve = function(x) f_dens(x,mean_beta0, mean_beta1, var_mean0, var_mean1, var_e, f0, f1, type)
integratefn_int_const = integrate(f = unnormalized_f_dens_curve, lower = -1, upper = 1)$value


f_dens_curve = function(x) (1/int_const) * unnormalized_f_dens_curve(x)
f_dens_curve2 = function(x) (1/integratefn_int_const) * unnormalized_f_dens_curve(x)



###
hist(mmed, probability = T, ylim = c(0, 2.5), breaks = 20)
curve(f_dens_curve, add = T)
curve(f_dens_curve2, add = T)
###
hist(mmed2, probability = T, ylim = c(0, 2.5), breaks = 20)
curve(f_dens_curve, add = T)
curve(f_dens_curve2, add = T)
###
hist(mmed3, probability = T, ylim = c(0, 2.5), breaks = 20)
curve(f_dens_curve, add = T)
curve(f_dens_curve2, add = T)

# saveRDS(mmed3, "mmed_N300.rds")
# saveRDS(mmed, "mmed_N100.rds")










