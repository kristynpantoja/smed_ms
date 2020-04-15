# things to do:
# 1. gp sequential
#     issue: expand.grid() vs. unname(as.matrix(expand.grid())): 
#           give different results in wasserstein fn
# 2. rerun run_designs_lmvs4 but with N_seq = 27 (currently running)


# --- Working Directory --- #
home = "/home/kristyn/Documents/smed_ms"

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

# read in grid info
grid12 = readRDS(paste(home, "/grid_gpvs12.rds", sep = ""))
# the grid of inputs
numx = grid12[[1]]
x_seq = grid12[[2]]
x_grid = as.matrix(grid12[[3]])
null_cov = grid12[[4]]
null_mean = grid12[[5]]

############################################################
# x_grid2 = expand.grid(x_seq, x_seq)
# null_cov2 = getCov(x_seq, x_seq, type01[2], l01[2])
# null_mean2 = rep(0, numx)
############################################################

# --- simulations  --- #
numSims = 1

# x_seq, grid over which to generate subsequent functions
xmin = 0; xmax = 1
l01= c(0.1, 0.1)
type01 = c(1, 1)

# input set designs
N = 6
p = 1
k = 4 * p
alpha = 1
subdim = 1

## smmed settings
N2 = 7
subdim = 1
nugget2 = 1e-1

# input points (randomly selected)
# x_input
set.seed(1997)
x_input_ind = sample(1:dim(x_grid)[1], N)
x_input = x_grid[x_input_ind, ]
# y_seq and y_input
set.seed(125243363)
y_seq_mat_1d = t(rmvnorm(n = numSims, mean = null_mean, sigma = null_cov)) # the function values
y_seq_mat = matrix(NA, nrow = numx^2, ncol = numSims) # expand to 2 dims, to have a y for each x in x_grid
for(i in 1:numSims){
  expanded = expand.grid(y_seq_mat_1d[ , i], y_seq_mat_1d[ , i])
  y_seq_mat[ , i] = expanded[ , 1] # only need 1st column, since restricting to V1
}
sim_index = 1
y_seq = y_seq_mat[ , sim_index]
y_input = y_seq[x_input_ind]

# generate mmed
# mmedgp = add_MMEDgpvs_oneatatime(x_input, y_input, type01, l01, subdim = subdim, var_e = 1, N2 = N2, k = k, p = p, xmin = xmin, xmax = xmax, nugget = nugget2, alpha = alpha, buffer = 0, candidates = x_grid)
mmedgp2 = add_MMEDgpvs_oneatatime2(x_input, y_input, type01, l01, subdim = subdim, var_e = 1, N2 = N2, k = k, p = p, xmin = xmin, xmax = xmax, nugget = nugget2, alpha = alpha, buffer = 0, candidates = x_grid)
# initD = x_input
# y = y_input
# type = type01
# l = l01
# subdim = subdim
# var_e = 1
# N2 = 7
# numCandidates = NULL
# # k = 4
# # p = 1
# # xmin = 0
# # xmax = 1
# nugget = nugget2
# # alpha = alpha
# buffer = 0
# genCandidates = 1
# candidates = x_grid

# true f (in 1d)
order_truef = order(x_grid[ , 1])
plot(x_grid[order_truef,1], y_seq[order_truef], type = "l")
points(x_grid[ x_input_ind, 1], y_seq[ x_input_ind])
points(x_grid[ mmedgp$indices, 1], y_seq[mmedgp$indices], pch = 3, col = 2)
points(x_grid[ mmedgp2$indices, 1], y_seq[mmedgp2$indices], pch = 3, col = 4)

# f in 2d
quilt.plot(x_grid, y_seq, main = "MMED")
points(x_grid[ mmedgp$indices, 1], x_grid[ mmedgp$indices, 2], pch = "o", cex = 1.5, col = "magenta")
points(x_grid[ x_input_ind, 1], x_grid[ x_input_ind, 2], pch = "o", cex = 1.5, col = "darkgreen")

# wasserstein(x)
initD = x_input
subinitD = initD[ , subdim, drop = FALSE]
Kinv0 = solve(getCov(subinitD, subinitD, type01[1], l01[1]) + diag(rep(nugget2, dim(initD)[1])))
Kinv1 = solve(getCov(initD, initD, type01[2], l01[2]) + diag(rep(nugget2, dim(initD)[1])))
w_seq2 = apply(x_grid, 1, FUN = function(x) Wasserstein_distance_postpred_gpvs2(x, Kinv0, Kinv1, subdim, subinitD, initD, y_input, var_e = 1, type01, l01))
# w_seq = apply(x_grid, 1, FUN = function(x) Wasserstein_distance_postpred_gpvs(x, Kinv0, Kinv1, subdim, subinitD, initD, y_input, var_e = 1, type01, l01))

par(mfrow=c(1,2))
# quilt.plot(x_grid, w_seq, main = "W(x)", xlab = "", ylab = "")
# points(x_grid[ mmedgp$indices, 1], x_grid[ mmedgp$indices, 2], pch = "o", cex = 1.5, col = "magenta")
# points(x_grid[ x_input_ind, 1], x_grid[ x_input_ind, 2], pch = "o", cex = 1.5, col = "darkgreen")
#
quilt.plot(x_grid, w_seq2, main = "W(x)", xlab = "", ylab = "")
points(x_grid[ mmedgp2$indices, 1], x_grid[ mmedgp2$indices, 2], pch = "o", cex = 1.5, col = "magenta")
points(x_grid[ x_input_ind, 1], x_grid[ x_input_ind, 2], pch = "o", cex = 1.5, col = "darkgreen")
#
w_max = which.max(w_seq)
points(x_grid[ w_max, 1], x_grid[ w_max, 2], pch = "x", cex = 1.5, col = 2)
text(x_grid[ mmedgp$indices, 1], x_grid[ mmedgp$indices, 2], 1:7)
w_order = order(w_seq, decreasing = TRUE)
w_order[1:7]
mmedgp$indices
mmedgp2$indices
points(x_grid[ w_order[1:20], 1], x_grid[ w_order[1:20], 2], pch = "^")
### testing Wasserstein function


