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


mean_beta0 = c(0, 1 / 2) # slope of null model
mean_beta1 = c(0, 1 / 4) # slope of alternative model
var_mean0 = diag(c(0.005, 0.005)); var_mean1 = var_mean0 # variance on beta
var_e = 0.025 # variance on error
xmin = 0
xmax = 1
f0 = function(x) mean_beta0[1] + mean_beta0[2] * x # null regression model
f1 = function(x) mean_beta1[1] + mean_beta1[2] * x # alternative regression model
N = 300
type = c(2, 2)
p = 2
# for fast algorithm:
S = 5
fast_N300 = MED_ms_fast(mean_beta0, mean_beta1, var_mean0, var_mean1, var_e,
                        f0 = f0, f1 = f1, type = type,
                        N = N, K = S, xmin = xmin, xmax = xmax, p = p, alpha = 1,
                        genCandidates = 1, initialpt = 1)
# for one-at-a-time algorithm:
numCandidates = 10^4
k = 4
greedy_N300 = MED_ms_oneatatime(mean_beta0, mean_beta1, var_mean0, var_mean1, var_e,
                         f0, f1, type, N = N, numCandidates = numCandidates, k = k,
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
hist(fast_N300$D[,S], probability = T, ylim = c(0, 2.5), breaks = 10)
curve(f_dens_curve2, add = T)
hist(greedy_N300, probability = T, ylim = c(0, 2.5), breaks = 10)
curve(f_dens_curve2, add = T)




saveRDS(fast_N300, "mmed_fast_N300_S5.rds")
saveRDS(greedy_N300, "mmed_greedy_N300_k4.rds")
# 
# 
# ###############
# ## TRY AGAIN ##
# ###############
# 
# 
# 
# mean_beta0 = c(0, 1 / 2) # slope of null model
# mean_beta1 = c(0, 1 / 4) # slope of alternative model
# var_mean0 = diag(c(0.005, 0.005)); var_mean1 = var_mean0 # variance on beta
# var_e = 0.025 # variance on error
# xmin = 0
# xmax = 1
# f0 = function(x) mean_beta0[1] + mean_beta0[2] * x # null regression model
# f1 = function(x) mean_beta1[1] + mean_beta1[2] * x # alternative regression model
# N = 300
# type = c(2, 2)
# p = 2
# # for fast algorithm:
# S = 5
# numCandidats = 1000
# fast_N300v2 = MED_ms_fast(mean_beta0, mean_beta1, var_mean0, var_mean1, var_e,
#                         f0 = f0, f1 = f1, type = type, 
#                         N = N, K = S, xmin = xmin, xmax = xmax, p = p, alpha = 1, 
#                         numCandidates = numCandidates, genCandidates = 1, initialpt = 1)
# # for one-at-a-time algorithm:
# numCandidates = 10^4
# k = 8
# greedy_N300v2 = MED_ms_oneatatime(mean_beta0, mean_beta1, var_mean0, var_mean1, var_e, 
#                                 f0, f1, type, N = N, numCandidates = numCandidates, k = k, 
#                                 xmin = xmin, xmax = xmax, p = 2, alpha = 1, buffer = 0, 
#                                 genCandidates = 1, initialpt = 1)
# ###
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
# hist(fast_N300v2$D[,S], probability = T, ylim = c(0, 2.5), breaks = 10)
# curve(f_dens_curve2, add = T)
# hist(greedy_N300v2, probability = T, ylim = c(0, 2.5), breaks = 10)
# curve(f_dens_curve2, add = T)
# 
# 
# 
# 
# saveRDS(fast_N300v2, "mmed_fast_N300_S5_nC1k.rds")
# saveRDS(greedy_N300v2, "mmed_greedy_N300_k8.rds")
# 
# 
# 
# 
# 
# 
# 
# ###############
# ## TRY AGAIN ##
# ###############
# 
# 
# 
# mean_beta0 = c(0, 1 / 2) # slope of null model
# mean_beta1 = c(0, 1 / 4) # slope of alternative model
# var_mean0 = diag(c(0.005, 0.005)); var_mean1 = var_mean0 # variance on beta
# var_e = 0.025 # variance on error
# xmin = 0
# xmax = 1
# f0 = function(x) mean_beta0[1] + mean_beta0[2] * x # null regression model
# f1 = function(x) mean_beta1[1] + mean_beta1[2] * x # alternative regression model
# N = 300
# type = c(2, 2)
# p = 2
# numCandidates = 10^4
# k = 25
# greedy_N300v2.1 = MED_ms_oneatatime(mean_beta0, mean_beta1, var_mean0, var_mean1, var_e, 
#                                   f0, f1, type, N = N, numCandidates = numCandidates, k = k, 
#                                   xmin = xmin, xmax = xmax, p = 2, alpha = 1, buffer = 0, 
#                                   genCandidates = 1, initialpt = 1)
# ###
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
# 
# #  check integrates to 1
# integraldenscurve2 = integrate(f = f_dens_curve2, lower = xmin, upper = xmax)$value
# 
# 
# ###
# 
# hist(greedy_N300v2.1, probability = T, ylim = c(0, 2.5), breaks = 10)
# curve(f_dens_curve2, add = T)
# 
# 
# saveRDS(greedy_N300v2.1, "mmed_greedy_N300_k25.rds")
# 
# 
# 
# 
# 
# ###############
# ## TRY AGAIN ##
# ###############
# 
# 
# 
# mean_beta0 = c(0, 1 / 2) # slope of null model
# mean_beta1 = c(0, 1 / 4) # slope of alternative model
# var_mean0 = diag(c(0.005, 0.005)); var_mean1 = var_mean0 # variance on beta
# var_e = 0.025 # variance on error
# xmin = 0
# xmax = 1
# f0 = function(x) mean_beta0[1] + mean_beta0[2] * x # null regression model
# f1 = function(x) mean_beta1[1] + mean_beta1[2] * x # alternative regression model
# N = 500
# type = c(2, 2)
# p = 2
# numCandidates = 10^4
# k = 25
# greedy_N500v2.2 = MED_ms_oneatatime(mean_beta0, mean_beta1, var_mean0, var_mean1, var_e, 
#                                     f0, f1, type, N = N, numCandidates = numCandidates, k = k, 
#                                     xmin = xmin, xmax = xmax, p = 2, alpha = 1, buffer = 0, 
#                                     genCandidates = 1, initialpt = 1)
# ###
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
# 
# #  check integrates to 1
# integraldenscurve2 = integrate(f = f_dens_curve2, lower = xmin, upper = xmax)$value
# 
# 
# ###
# 
# hist(greedy_N500v2.2, probability = T, ylim = c(0, 2.5), breaks = 10)
# curve(f_dens_curve2, add = T)
# 
# 
# saveRDS(greedy_N500v2.2, "mmed_greedy_N500_k25.rds")
# 
# 
# 
# 
# 
# 
# ###############
# ## TRY AGAIN ##
# ###############
# 
# 
# 
# mean_beta0 = c(0, 1 / 2) # slope of null model
# mean_beta1 = c(0, 1 / 4) # slope of alternative model
# var_mean0 = diag(c(0.005, 0.005)); var_mean1 = var_mean0 # variance on beta
# var_e = 0.025 # variance on error
# xmin = 0
# xmax = 1
# f0 = function(x) mean_beta0[1] + mean_beta0[2] * x # null regression model
# f1 = function(x) mean_beta1[1] + mean_beta1[2] * x # alternative regression model
# N = 500
# type = c(2, 2)
# p = 2
# # for one-at-a-time algorithm:
# numCandidates = 10^4
# k = 8
# greedy_N500v2.3 = MED_ms_oneatatime(mean_beta0, mean_beta1, var_mean0, var_mean1, var_e, 
#                                   f0, f1, type, N = N, numCandidates = numCandidates, k = k, 
#                                   xmin = xmin, xmax = xmax, p = 2, alpha = 1, buffer = 0, 
#                                   genCandidates = 1, initialpt = 1)
# ###
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
# hist(greedy_N500v2.3, probability = T, ylim = c(0, 2.5), breaks = 10)
# curve(f_dens_curve2, add = T)
# 
# saveRDS(greedy_N500v2.3, "mmed_greedy_N500_k8.rds")
# 
# 
# 
# 
# 
# 
# 
###############
## TRY AGAIN ##
###############



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
p = 2
# for one-at-a-time algorithm:
numCandidates = 10^4
k = 4
greedy_N500v2.4 = MED_ms_oneatatime(mean_beta0, mean_beta1, var_mean0, var_mean1, var_e,
                                    f0, f1, type, N = N, numCandidates = numCandidates, k = k,
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
hist(greedy_N500v2.4, probability = T, ylim = c(0, 2.5), breaks = 10)
curve(f_dens_curve2, add = T)

saveRDS(greedy_N500v2.4, "mmed_greedy_N500_k4.rds")










###############
## TRY AGAIN ##
###############



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
p = 2
# for one-at-a-time algorithm:
numCandidates = 10^4
k = 1
greedy_N500v2.4.1 = MED_ms_oneatatime(mean_beta0, mean_beta1, var_mean0, var_mean1, var_e,
                                    f0, f1, type, N = N, numCandidates = numCandidates, k = k,
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
hist(greedy_N500v2.4.1, probability = T, ylim = c(0, 2.5), breaks = 10)
curve(f_dens_curve2, add = T)

saveRDS(greedy_N500v2.4.1, "mmed_greedy_N500_k1.rds")
# 
# 
# 
# 
# 
# 
# 
# 
# 
# ###############
# ## TRY AGAIN ##
# ###############
# # Dr. Jones' uniformly distributed jitter
# 
# 
# mean_beta0 = c(0, 1 / 2) # slope of null model
# mean_beta1 = c(0, 1 / 4) # slope of alternative model
# var_mean0 = diag(c(0.005, 0.005)); var_mean1 = var_mean0 # variance on beta
# var_e = 0.025 # variance on error
# xmin = 0
# xmax = 1
# f0 = function(x) mean_beta0[1] + mean_beta0[2] * x # null regression model
# f1 = function(x) mean_beta1[1] + mean_beta1[2] * x # alternative regression model
# N = 500
# type = c(2, 2)
# p = 2
# # for one-at-a-time algorithm:
# numCandidates = 10^4
# k = 4
# greedy_N500v2.5 = MED_ms_oneatatime(mean_beta0, mean_beta1, var_mean0, var_mean1, var_e, 
#                                     f0, f1, type, N = N, numCandidates = numCandidates, k = k, 
#                                     xmin = xmin, xmax = xmax, p = 2, alpha = 1, buffer = 0, 
#                                     genCandidates = 1, initialpt = 1, jitter = TRUE)
# ###
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
# hist(greedy_N500v2.5, probability = T, ylim = c(0, 2.5), breaks = 10)
# curve(f_dens_curve2, add = T)
# 
# saveRDS(greedy_N500v2.5, "mmed_greedy_N500_k4_jitter.rds")
# 
# 
# 
# ###############
# ## TRY AGAIN ##
# ###############
# # Dr. Jones' uniformly distributed jitter
# 
# 
# mean_beta0 = c(0, 1 / 2) # slope of null model
# mean_beta1 = c(0, 1 / 4) # slope of alternative model
# var_mean0 = diag(c(0.005, 0.005)); var_mean1 = var_mean0 # variance on beta
# var_e = 0.025 # variance on error
# xmin = 0
# xmax = 1
# f0 = function(x) mean_beta0[1] + mean_beta0[2] * x # null regression model
# f1 = function(x) mean_beta1[1] + mean_beta1[2] * x # alternative regression model
# N = 500
# type = c(2, 2)
# p = 2
# # for one-at-a-time algorithm:
# numCandidates = 10^4
# k = 8
# greedy_N500v2.6= MED_ms_oneatatime(mean_beta0, mean_beta1, var_mean0, var_mean1, var_e, 
#                                     f0, f1, type, N = N, numCandidates = numCandidates, k = k, 
#                                     xmin = xmin, xmax = xmax, p = 2, alpha = 1, buffer = 0, 
#                                     genCandidates = 1, initialpt = 1, jitter = TRUE)
# ###
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
# hist(greedy_N500v2.6, probability = T, ylim = c(0, 2.5), breaks = 10)
# curve(f_dens_curve2, add = T)
# 
# saveRDS(greedy_N500v2.6, "mmed_greedy_N500_k8_jitter.rds")
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# ###############
# ## TRY AGAIN ##
# ###############
# # Dr. Jones' uniformly distributed jitter
# 
# 
# mean_beta0 = c(0, 1 / 2) # slope of null model
# mean_beta1 = c(0, 1 / 4) # slope of alternative model
# var_mean0 = diag(c(0.005, 0.005)); var_mean1 = var_mean0 # variance on beta
# var_e = 0.025 # variance on error
# xmin = 0
# xmax = 1
# f0 = function(x) mean_beta0[1] + mean_beta0[2] * x # null regression model
# f1 = function(x) mean_beta1[1] + mean_beta1[2] * x # alternative regression model
# N = 500
# type = c(2, 2)
# p = 2
# # for one-at-a-time algorithm:
# numCandidates = 10^4
# k = 4
# greedy_N500v2.7 = MED_ms_oneatatime(mean_beta0, mean_beta1, var_mean0, var_mean1, var_e, 
#                                     f0, f1, type, N = N, numCandidates = numCandidates, k = k, 
#                                     xmin = xmin, xmax = xmax, p = 2, alpha = 1, buffer = 0, 
#                                     genCandidates = 1, initialpt = 1, jitter = TRUE)
# ###
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
# hist(greedy_N500v2.7, probability = T, ylim = c(0, 2.5), breaks = 10)
# curve(f_dens_curve2, add = T)
# 
# saveRDS(greedy_N500v2.7, "mmed_greedy_N500_k4_jitter2.rds")
# 
# 
# 
# ###############
# ## TRY AGAIN ##
# ###############
# # Dr. Jones' uniformly distributed jitter
# 
# 
# mean_beta0 = c(0, 1 / 2) # slope of null model
# mean_beta1 = c(0, 1 / 4) # slope of alternative model
# var_mean0 = diag(c(0.005, 0.005)); var_mean1 = var_mean0 # variance on beta
# var_e = 0.025 # variance on error
# xmin = 0
# xmax = 1
# f0 = function(x) mean_beta0[1] + mean_beta0[2] * x # null regression model
# f1 = function(x) mean_beta1[1] + mean_beta1[2] * x # alternative regression model
# N = 500
# type = c(2, 2)
# p = 2
# # for one-at-a-time algorithm:
# numCandidates = 10^4
# k = 8
# greedy_N500v2.8= MED_ms_oneatatime(mean_beta0, mean_beta1, var_mean0, var_mean1, var_e, 
#                                    f0, f1, type, N = N, numCandidates = numCandidates, k = k, 
#                                    xmin = xmin, xmax = xmax, p = 2, alpha = 1, buffer = 0, 
#                                    genCandidates = 1, initialpt = 1, jitter = TRUE)
# ###
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
# hist(greedy_N500v2.8, probability = T, ylim = c(0, 2.5), breaks = 10)
# curve(f_dens_curve2, add = T)
# 
# saveRDS(greedy_N500v2.8, "mmed_greedy_N500_k8_jitter2.rds")










###########################################################################
# does mmed with N300 match first 300 pts in mmed with N500?
###########################################################################


mmed_lim1 = readRDS("/Users/kristyn/Documents/research/smed_ms/limiting_distribution/mmed_greedy_N300_k4.rds")
hist(mmed_lim, breaks = 15, ylim = c(0, 2), probability = T, 
     main = "MED, N = 300, q = 1/W^(1/2p)", xlab = "design points")
curve(f_dens_curve2, add = T)

mmed_lim2 = readRDS("/Users/kristyn/Documents/research/smed_ms/limiting_distribution/mmed_greedy_N500_k4.rds")
hist(mmed_lim, breaks = 15, ylim = c(0, 2), probability = T, 
     main = "MED, N = 500, q = 1/W^(1/2p)", xlab = "design points")
curve(f_dens_curve2, add = T)

length(mmed_lim1)
length(mmed_lim2)
all(mmed_lim2[1:300] == mmed_lim1)

hist(mmed_lim2[1:50], breaks = 15, ylim = c(0, 2), probability = T)
hist(mmed_lim2[1:100], breaks = 15, ylim = c(0, 2), probability = T)
hist(mmed_lim2[1:150], breaks = 15, ylim = c(0, 2), probability = T)
hist(mmed_lim2[1:200], breaks = 15, ylim = c(0, 2), probability = T)
hist(mmed_lim2[1:250], breaks = 15, ylim = c(0, 2), probability = T)
hist(mmed_lim2[1:300], breaks = 15, ylim = c(0, 2), probability = T)
hist(mmed_lim2[1:350], breaks = 15, ylim = c(0, 2), probability = T)
hist(mmed_lim2[1:400], breaks = 15, ylim = c(0, 2), probability = T)
hist(mmed_lim2[1:450], breaks = 15, ylim = c(0, 2), probability = T)
hist(mmed_lim2[1:500], breaks = 15, ylim = c(0, 2), probability = T)

mmed_lim3 = readRDS("/Users/kristyn/Documents/research/smed_ms/limiting_distribution/mmed_greedy_N500_k1.rds")
hist(mmed_lim, breaks = 15, ylim = c(0, 2), probability = T, 
     main = "MED, N = 500, q = 1/W^(1/2p)", xlab = "design points")
curve(f_dens_curve2, add = T)

hist(mmed_lim3[1:50], breaks = 15, ylim = c(0, 2), probability = T)
hist(mmed_lim3[1:100], breaks = 15, ylim = c(0, 2), probability = T)
hist(mmed_lim3[1:150], breaks = 15, ylim = c(0, 2), probability = T)
hist(mmed_lim3[1:200], breaks = 15, ylim = c(0, 2), probability = T)
hist(mmed_lim3[1:250], breaks = 15, ylim = c(0, 2), probability = T)
hist(mmed_lim3[1:300], breaks = 15, ylim = c(0, 2), probability = T)
hist(mmed_lim3[1:350], breaks = 15, ylim = c(0, 2), probability = T)
hist(mmed_lim3[1:400], breaks = 15, ylim = c(0, 2), probability = T)
hist(mmed_lim3[1:450], breaks = 15, ylim = c(0, 2), probability = T)
hist(mmed_lim3[1:500], breaks = 15, ylim = c(0, 2), probability = T)
curve(f_dens_curve2, add = T)


###########################################################################
# the optimize q thing
###########################################################################

home = "/Users/kristyn/Documents/research/smed_ms"
functions_home = paste(home, "/functions", sep="")
source(paste(functions_home, "/variance_marginal_y.R", sep = ""))
source(paste(functions_home, "/wasserstein_distance.R", sep = ""))
source(paste(functions_home, "/wasserstein_distance.R", sep = ""))
source(paste(functions_home, "/charge_function_q.R", sep = ""))
source(paste(functions_home, "/generate_MED_oneatatime.R", sep = ""))

# mmed_lim3[1:488] # 0.31313131 picked ind 488, with value 0.31313131, as limit bc
# next point is 0.14911491, even though we think it should probably be higher for some reason
mean_beta0 = c(0, 1 / 2) # slope of null model
mean_beta1 = c(0, 1 / 4) # slope of alternative model
var_beta0 = diag(c(0.005, 0.005)); var_beta1 = var_beta0 # variance on beta
var_e = 0.025 # variance on error
xmin = 0
xmax = 1
f0 = function(x) mean_beta0[1] + mean_beta0[2] * x # null regression model
f1 = function(x) mean_beta1[1] + mean_beta1[2] * x # alternative regression model
N = 500
type = c(2, 2)
p = 2
# for one-at-a-time algorithm:
numCandidates = 10^4
k = 1
D = mmed_lim3[1:488]
candidates = seq(from = xmin, to = xmax, length.out = numCandidates)
f_min_candidates = sapply(candidates, function(x) f_min(x, D, k, mean_beta0, mean_beta1, 
                                                        var_beta0, var_beta1, var_e, f0, f1, 
                                                        type, var_margy0 = NULL, var_margy1 = NULL, 
                                                        p, alpha = 1, buffer = 0, log_space = FALSE))
mmed_lim3[489] == candidates[which.min(f_min_candidates)]
plot(x = candidates, y = f_min_candidates, type = "l")
# it looks like it becomes a matter of balancing the empty areas, in the places where it should be more sparse,
# i.e. closer to 0, with the areas that are more interesting, i.e. closer to 1.


# for one-at-a-time algorithm:
hist(mmed_lim2)
numCandidates = 10^4
candidates = seq(from = xmin, to = xmax, length.out = numCandidates)
k = 4
D = mmed_lim2[1:402]
f_min_candidates2 = sapply(candidates, function(x) f_min(x, D, k, mean_beta0, mean_beta1, 
                                                        var_beta0, var_beta1, var_e, f0, f1, 
                                                        type, var_margy0 = NULL, var_margy1 = NULL, 
                                                        p, alpha = 1, buffer = 0, log_space = FALSE))
mmed_lim2[403] == candidates[which.min(f_min_candidates2)]
plot(x = candidates, y = f_min_candidates2, type = "l")
