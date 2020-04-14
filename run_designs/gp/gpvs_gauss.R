# gaussian process variable selection between 1D and 2D with no noise and with noise
# testing higher dimensions

# --- Working Directory --- #

# Computer
home = "/home/kristyn/Documents/smed_ms"
output_home = paste(home, "/", sep = "")

# Cluster
# home = "/scratch/user/kristynp/smed_ms"
# output_home = paste(home,"/run_designs_gpvs3/",sep="")

# --- Sources/Libraries --- #
functions_home = paste(home, "/functions", sep="")
source(paste(functions_home, "/generate_MMEDgp_oneatatime.R", sep = ""))
source(paste(functions_home, "/wasserstein_distance.R", sep = ""))
source(paste(functions_home, "/charge_function_q.R", sep = ""))
source(paste(functions_home, "/covariance_functions.R", sep = ""))
source(paste(functions_home, "/gp_predictive.R", sep = ""))

library(expm)
library(matrixStats)
library(MASS)
library(mvtnorm)

# --- simulations  --- #
numSims = 1

# settings
xmin = 0; xmax = 1
l01= c(0.1, 0.1)
type01 = c(1, 1)

N = 6
p = 2
k = 4 * p
alpha = 1
subdim = 1

# the grid of inputs
numx = 31
x_seq = seq(from = xmin, to = xmax, length.out = numx)
x_grid = expand.grid(x_seq, x_seq)

# generate matern functions (1d f)
set.seed(1)
null_cov1d = getCov(x_seq, x_seq, type01[2], l01[2])
null_mean1d = rep(0, numx)
y_seq_mat1d_raw = t(rmvnorm(n = numSims, mean = null_mean1d, sigma = null_cov1d)) # the function values
# expand to 2 dims
y_seq_mat1d = matrix(NA, nrow = numx^2, ncol = numSims)
for(i in 1:numSims){
  expanded = expand.grid(y_seq_mat1d_raw[ , i], y_seq_mat1d_raw[ , i])
  y_seq_mat1d[ , i] = expanded[ , 1] # only need 1st column, since restricting to V1
}

# generate matern functions (2d f)
set.seed(1)
null_cov2d = getCov(x_grid, x_grid, type01[2], l01[2])
null_mean2d = rep(0, numx^2)
y_seq_mat2d = t(rmvnorm(n = numSims, mean = null_mean2d, sigma = null_cov2d)) # the function values

#############
# f 1d sims #
#############
y_seq_mat = y_seq_mat1d

f1dsims = list()
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

  # generate mmeds (no noise and noise cases)
  # nugget = NULL
  # f1dsims[[i]] = add_MED_ms_oneatatime_data_gp2(x_input, y_input, type01, l01, subdim = 1, var_e = 1, N2 = 7, k = k, p = p, xmin = xmin, xmax = xmax, nugget = nugget, alpha = alpha, buffer = 0, candidates = x_grid)
  nugget2 = 1e-1
  f1dsims_noise[[i]] = add_MED_ms_oneatatime_data_gp2(x_input, y_input, type01, l01, subdim = 1, var_e = 1, N2 = 7, k = k, p = p, xmin = xmin, xmax = xmax, nugget = nugget2, alpha = alpha, buffer = 0, candidates = x_grid)
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
  
  # generate mmeds (no noise and noise cases)
  # nugget = NULL
  # f2dsims[[i]] = add_MED_ms_oneatatime_data_gp2(x_input, y_input, type01, l01, subdim = 1, var_e = 1, N2 = 7, k = k, p = p, xmin = xmin, xmax = xmax, nugget = nugget, alpha = alpha, buffer = 0, candidates = x_grid)
  nugget2 = 1e-1
  f2dsims_noise[[i]] = add_MED_ms_oneatatime_data_gp2(x_input, y_input, type01, l01, subdim = 1, var_e = 1, N2 = 7, k = k, p = p, xmin = xmin, xmax = xmax, nugget = nugget2, alpha = alpha, buffer = 0, candidates = x_grid)
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






