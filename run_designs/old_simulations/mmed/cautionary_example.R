# Originally ran in run_designs/limiting_distribution/limiting_distribution3.R (bottom of script, with p = 2 though)

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

mu0 = c(0, 0)
mu1 = c(0, 0, 0)
typeT = 3
betaT = c(-0.2, -0.4, 0.4)
sigmasq01 = 0.25
sigmasq = 0.1

# MED design #
typeT = 3
xmin = -1
xmax = 1
# MED design #
V0 = diag(rep(sigmasq01,length(mu0)))
V1 = diag(rep(sigmasq01,length(mu1)))
# function settings (including and based on prior settings above)
f0 = function(x) mu0[1] + mu0[2] * x
f1 = function(x) mu1[1] + mu1[2] * x + mu1[3] * x^2
type01 = c(2, 3)
numCandidates = 10^4
k = 4
xmin = -1
xmax = 1
p = 1
N = 100

# alpha = 1, i.e. original, more space-filling-design
mmed = MED_ms_oneatatime(mu0, mu1, V0, V1, sigmasq, f0, f1, type01, N, 
                           numCandidates, k, xmin, xmax, p = p, alpha = 1)

saveRDS(mmed, paste(output_home,"cautionary_example.rds", sep = ""))
###
# f_dens = function(x, mean_beta0, mean_beta1, var_mean0, var_mean1, var_e, f0, f1, type){
#   if(length(type) != 2) stop("type should be vector with length == 2")
#   mu1 = f0(x) # mean of marginal dist of y | H0
#   mu2 = f1(x) # mean of marginal dist of y | H1
#   var1 = var_marginaly(x, var_mean0, var_e, type = type[1]) # variance of marginal dist of y | H0
#   var2 = var_marginaly(x, var_mean1, var_e, type = type[2]) # variance of marginal dist of y | H1
#   Wass_dist = Wasserstein_distance(mu1, mu2, var1, var2)
#   return(Wass_dist)
# }
# 
# unnormalized_f_dens_curve = function(x) f_dens(x, mu0, mu1, V0, V1, sigmasq, f0, f1, type01)
# integratefn_int_const = integrate(f = unnormalized_f_dens_curve, lower = -1, upper = 1)$value
# 
# f_dens_curve2 = function(x) (1/integratefn_int_const) * unnormalized_f_dens_curve(x)
# 
# ###
# hist(mmedv2, probability = T, ylim = c(0, 2.5), breaks = 20)
# curve(f_dens_curve2, add = T)