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




### RAN. ###
N = 500
numCandidates = 10^4
k = 2
p = 1
# --- example 3 --- #

sigmasq01 = 0.005
mu0 = c(0, 0)
V0 = diag(rep(sigmasq01,length(mu0)))
mu1 = c(1, 0, -1)
V1 = diag(rep(sigmasq01,length(mu1)))
sigmasq = 0.025
# function settings (including and based on prior settings above)
f0 = function(x) mu0[1] + mu0[2] * x
f1 = function(x) mu1[1] + mu1[2] * x + mu1[3] * x^2
type = c(2, 3)
xmin = -1
xmax = 1
alpha = 1
greedy_N500 = MED_ms_oneatatime(mu0, mu1, V0, V1, sigmasq, f0, f1, type, N,
                                          numCandidates, k, xmin, xmax, p = p, alpha = alpha)

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
unnormalized_f_dens_curve = function(x) f_dens(x,mu0, mu1, V0, V1, sigmasq, f0, f1, type)
integratefn_int_const = integrate(f = unnormalized_f_dens_curve, lower = xmin, upper = xmax)$value
f_dens_curve2 = function(x) (1/integratefn_int_const) * unnormalized_f_dens_curve(x)
###
hist(greedy_N500, probability = T, ylim = c(0, 2.5), breaks = 10)
curve(f_dens_curve2, add = T)
saveRDS(greedy_N500, "mmed_greedy_N500_k2p1_ex3.rds")




hist(greedy_N500[1:50], probability = T, ylim = c(0, 2.5), breaks = 10)
curve(f_dens_curve2, add = T)
hist(greedy_N500[1:100], probability = T, ylim = c(0, 2.5), breaks = 10)
curve(f_dens_curve2, add = T)
hist(greedy_N500[1:150], probability = T, ylim = c(0, 2.5), breaks = 10)
curve(f_dens_curve2, add = T)
hist(greedy_N500[1:200], probability = T, ylim = c(0, 2.5), breaks = 10)
curve(f_dens_curve2, add = T)
hist(greedy_N500[1:250], probability = T, ylim = c(0, 2.5), breaks = 10)
curve(f_dens_curve2, add = T)
hist(greedy_N500[1:300], probability = T, ylim = c(0, 2.5), breaks = 10)
curve(f_dens_curve2, add = T)
hist(greedy_N500[1:350], probability = T, ylim = c(0, 2.5), breaks = 10)
curve(f_dens_curve2, add = T)
hist(greedy_N500[1:400], probability = T, ylim = c(0, 2.5), breaks = 10)
curve(f_dens_curve2, add = T)
hist(greedy_N500[1:450], probability = T, ylim = c(0, 2.5), breaks = 10)
curve(f_dens_curve2, add = T)
hist(greedy_N500[1:500], probability = T, ylim = c(0, 2.5), breaks = 10)
curve(f_dens_curve2, add = T)

#########
#########
#########
#########
#########
#########
#########
#########
#########
#########
#########
#########
#########
#########
#########



#########
#########
#########
#########
#########
#########
#########
#########
#########
#########
#########
#########
#########
#########
#########



###############
## TRY AGAIN ##
###############


### NOT RAN. ###
mean_beta0 = c(0, 1 / 2) # slope of null model
mean_beta1 = c(0, 1 / 4) # slope of alternative model
var_mean0 = diag(c(0.005, 0.005)); var_mean1 = var_mean0 # variance on beta
var_e = 0.025 # variance on error
xmin = 0
xmax = 1
f0 = function(x) mean_beta0[1] + mean_beta0[2] * x # null regression model
f1 = function(x) mean_beta1[1] + mean_beta1[2] * x # alternative regression model
N = 500
type = c(2, 2)
# for one-at-a-time algorithm:
numCandidates = 10^4
k = 2
p = 1
greedy_N500v2.4 = MED_ms_oneatatime(mean_beta0, mean_beta1, var_mean0, var_mean1, var_e,
                                    f0, f1, type, N = N, numCandidates = numCandidates, k = k,
                                    xmin = xmin, xmax = xmax, p = p, alpha = 1, buffer = 0,
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
hist(greedy_N500v2.4, probability = T, ylim = c(0, 2.5), breaks = 10)
curve(f_dens_curve2, add = T)

saveRDS(greedy_N500v2.4, "mmed_greedy_N500_k4_p1.rds")

















