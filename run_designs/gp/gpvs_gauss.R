# gaussian process variable selection between 1D and 2D with no noise and with noise
# testing higher dimensions

# --- Working Directory --- #

# Computer
# home = "/home/kristyn/Documents/smed_ms"
# output_home = paste(home, "/", sep = "")

# Cluster
home = "/scratch/user/kristynp/smed_ms"
output_home = paste(home,"/run_designs_gpvs2/",sep="")

# --- Sources/Libraries --- #
functions_home = paste(home, "/functions", sep="")
source(paste(functions_home, "/generate_MMEDgp_oneatatime.R", sep = ""))
source(paste(functions_home, "/wasserstein_distance.R", sep = ""))
source(paste(functions_home, "/charge_function_q.R", sep = ""))
source(paste(functions_home, "/covariance_functions.R", sep = ""))
source(paste(functions_home, "/gp_predictive.R", sep = ""))
source(paste(functions_home, "/simulate_SMMED_gpvs.R", sep = ""))

library(expm)
library(matrixStats)
library(MASS)
library(mvtnorm)
library(fields)
library(knitr)

# --- simulations  --- #
numSims = 25

# settings
xmin = 0; xmax = 1
l01= c(0.1, 0.1)
type01 = c(1, 1)

N = 6
p = 2
k = 4 * p
alpha = 1
subdim = 1
N2 = 7

# read in grid info
grid12 = readRDS(paste(home, "/grid_gpvs12.rds", sep = ""))
# the grid of inputs
numx = grid12[[1]]
x_seq = grid12[[2]]
x_grid = as.matrix(grid12[[3]])
null_cov1d = grid12[[4]]
null_mean1d = grid12[[5]]
null_cov2d = grid12[[6]]
null_mean2d = grid12[[7]]

# generate matern functions (1d f)
set.seed(1)
y_seq_mat1d_raw = t(rmvnorm(n = numSims, mean = null_mean1d, sigma = null_cov1d)) # the function values
# expand to 2 dims
y_seq_mat1d = matrix(NA, nrow = numx^2, ncol = numSims)
for(i in 1:numSims){
  expanded = expand.grid(y_seq_mat1d_raw[ , i], y_seq_mat1d_raw[ , i])
  y_seq_mat1d[ , i] = expanded[ , 1] # only need 1st column, since restricting to V1
}

# generate matern functions (2d f)
set.seed(1)
y_seq_mat2d = t(rmvnorm(n = numSims, mean = null_mean2d, sigma = null_cov2d)) # the function values

#############
# f 1d sims #
#############
y_seq_mat = y_seq_mat1d

f1dsims_noise = list()
x_input_ind_mat = matrix(NA, N, numSims)

seed = 1
for(i in 1:numSims){
  set.seed(seed + i)
  y_seq = y_seq_mat[ , i]
  
  # get input points
  x_input_ind = sample(1:dim(x_grid)[1], N)
  x_input_ind_mat[ , i] = x_input_ind
  x_input = x_grid[x_input_ind, ]
  y_input = y_seq[x_input_ind]

  nugget2 = 1e-1
  f1dsims_noise[[i]] = add_MMEDgpvs_oneatatime(x_input, y_input, type01, l01, subdim = subdim, var_e = 1, N2 = N2, k = k, p = p, xmin = xmin, xmax = xmax, nugget = nugget2, alpha = alpha, buffer = 0, candidates = x_grid)
}



#############
# f 2d sims #
#############
y_seq_mat = y_seq_mat2d

f2dsims = list()
f2dsims_noise = list()
seed = 1
for(i in 1:numSims){
  set.seed(seed + i)
  y_seq = y_seq_mat[ , i]
  
  # get input points
  x_input_ind = x_input_ind_mat[ , i]
  x_input = x_grid[x_input_ind, ]
  y_input = y_seq[x_input_ind]
  
  nugget2 = 1e-1
  f2dsims_noise[[i]] = add_MMEDgpvs_oneatatime(x_input, y_input, type01, l01, subdim = 1, var_e = 1, N2 = 7, k = k, p = p, xmin = xmin, xmax = xmax, nugget = nugget2, alpha = alpha, buffer = 0, candidates = x_grid)
}

sims = list("x_input_ind" = x_input_ind_mat, 
            "null_mean1d" = null_mean1d, 
            "null_cov1d" = null_cov1d, 
            "y_seq_mat1d" = y_seq_mat1d, 
            "null_mean2d" = null_mean2d, 
            "null_cov2d" = null_cov2d, 
            "y_seq_mat2d" = y_seq_mat2d,
            "f1dsims_noise" = f1dsims_noise, 
            "f2dsims_noise" = f2dsims_noise)

saveRDS(sims, paste(output_home, "gpvs_gaussian_sims.rds", sep = ""))


# 
# # for 1d --- #
# sim_ind = 1
# mmedgp_sim = f1dsims_noise[[sim_ind]]
# true_y = y_seq_mat1d[ , sim_ind]
# x_input_ind = x_input_ind_mat[ , sim_ind]
# y_input = true_y[x_input_ind]
# 
# mmed_gp_vs = mmedgp_sim
# initD = x_grid[x_input_ind, ]
# subinitD = initD[ , subdim, drop = FALSE]
# Kinv0 = solve(getCov(subinitD, subinitD, type01[1], l01[1]) + diag(rep(nugget2, dim(initD)[1])))
# Kinv1 = solve(getCov(initD, initD, type01[2], l01[2]) + diag(rep(nugget2, dim(initD)[1])))
# w_seq = apply(x_grid, 1, FUN = function(x) Wasserstein_distance_postpred_gpvs(x, Kinv0, Kinv1, subdim, subinitD, initD, y_input, var_e = 1, type01, l01))
# quilt.plot(x_grid, w_seq, main = "W(x)", xlab = "", ylab = "")
# points(x_grid[ mmed_gp_vs$indices, 1], x_grid[ mmed_gp_vs$indices, 2], pch = "o", cex = 1.5, col = "magenta")
# points(x_grid[ x_input_ind, 1], x_grid[ x_input_ind, 2], pch = "o", cex = 1.5, col = "darkgreen")
# 
# # for 2d --- #
# sim_ind = 4
# mmedgp_sim = f2dsims_noise[[sim_ind]]
# true_y = y_seq_mat2d[ , sim_ind]
# x_input_ind = x_input_ind_mat[ , sim_ind]
# y_input = true_y[x_input_ind]
# 
# mmed_gp_vs = mmedgp_sim
# initD = x_grid[x_input_ind, ]
# subinitD = initD[ , subdim, drop = FALSE]
# Kinv0 = solve(getCov(subinitD, subinitD, type01[1], l01[1]) + diag(rep(nugget2, dim(initD)[1])))
# Kinv1 = solve(getCov(initD, initD, type01[2], l01[2]) + diag(rep(nugget2, dim(initD)[1])))
# w_seq = apply(x_grid, 1, FUN = function(x) Wasserstein_distance_postpred_gpvs(x, Kinv0, Kinv1, subdim, subinitD, initD, y_input, var_e = 1, type01, l01))
# quilt.plot(x_grid, w_seq, main = "W(x)", xlab = "", ylab = "")
# points(x_grid[ mmed_gp_vs$indices, 1], x_grid[ mmed_gp_vs$indices, 2], pch = "o", cex = 1.5, col = "magenta")
# points(x_grid[ x_input_ind, 1], x_grid[ x_input_ind, 2], pch = "o", cex = 1.5, col = "darkgreen")
