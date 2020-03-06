# gaussian process variable selection between 1D and 2D with no noise and with noise

# testing higher dimensions

# --- Working Directory --- #
home = "/home/kristyn/Documents/smed_ms"

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
library(fields)

# --- simulations  --- #

numSims = 3

# input set designs
xmin = 0; xmax = 1
l01= c(0.1, 0.1)
type01 = c(1, 1)
N = 6
k = 4
p = 1
alpha = 1

# create grid
numx2 = 21
x_seq2 = seq(from = xmin, to = xmax, length.out = numx2)
x_grid = expand.grid(x_seq2, x_seq2)

# generate fns! (2D)
seed = 1
set.seed(seed)
null_cov = getCov(x_grid, x_grid, type01[2], l01[2])
null_mean = rep(0, numx2^2)
y_seq_mat2d = t(rmvnorm(n = numSims, mean = null_mean, sigma = null_cov))

# generate fns! (1D)
seed = 1
set.seed(seed)
null_cov = getCov(x_seq2, x_seq2, type01[2], l01[2])
null_mean = rep(0, numx2)
y_seq_mat1d = t(rmvnorm(n = numSims, mean = null_mean, sigma = null_cov)) # the function values
# expand to 2 dims
y_seq_mat1d_expanded = matrix(NA, nrow = numx2^2, ncol = numSims)
for(i in 1:numSims){
  expanded = expand.grid(y_seq_mat1d[ , i], y_seq_mat1d[ , i])
  y_seq_mat1d_expanded[ , i] = expanded[ , 1] # only need 1st col, since only subsetting for V1
}

# apply mmed! (no noise)
nugget = NULL

# generate MMED, where true f is 1D
mmed_gp2d = list()
x_input_ind_mat = matrix(NA, nrow = N, ncol = numSims)
y_input_mat = matrix(NA, nrow = N, ncol = numSims)

# generate mmed
set.seed(1)
for(i in 1:numSims){
  # get x_input
  # can't select points that will have the same X1 value, if no noise!
  x_input_ind = sample(1:numx2, N) # so we do it this way
  x_input = x_grid[x_input_ind + numx2*(1:N), ]
  # get y_input
  y_seq = y_seq_mat1d_expanded[ , i]
  y_input = y_seq[x_input_ind]
  # save them
  x_input_ind_mat[ , i] = x_input_ind
  y_input_mat[ , i] = y_input
  # generate mmed
  mmed_gp2d[[i]] = add_MED_ms_oneatatime_data_gp2(x_input, y_input, type01, l01, subdim = 1, var_e = 1, 
                                                  N2 = 7, k = k, p = p, xmin = xmin, xmax = xmax, 
                                                  nugget = nugget, alpha = alpha, buffer = 0, 
                                                  candidates = x_grid)
}

gpvs12 = list("x_grid" = x_grid,
              "y_seq_mat" = y_seq_mat1d_expanded,
              "x_input_ind_mat" = x_input_ind_mat,
              "y_input_mat" = y_input_mat,
              "mmed_list" = mmed_gp2d)
saveRDS(gpvs12, file = "gpvs12_nonoise.rds")

# apply mmed! (yes noise)
nugget2 = 1e-1

# generate MMED, where true f is 1D
mmed_gp2d_noise = list()
x_input_ind_mat2 = matrix(NA, nrow = N, ncol = numSims)
y_input_mat2 = matrix(NA, nrow = N, ncol = numSims)

# generate mmed
set.seed(1)
for(i in 1:numSims){
  # get x_input
  # can't select points that will have the same X1 value, if no noise!
  x_input_ind2 = sample(1:dim(x_grid)[1], N)
  x_input2 = x_grid[x_input_ind2, ]
  # get y_input
  y_seq = y_seq_mat1d_expanded[ , i]
  y_input2 = y_seq[x_input_ind2]
  # save them
  x_input_ind_mat2[ , i] = x_input_ind2
  y_input_mat2[ , i] = y_input2
  # generate mmed
  mmed_gp2d_noise[[i]] = add_MED_ms_oneatatime_data_gp2(x_input2, y_input2, type01, l01, subdim = 1, var_e = 1, 
                                                  N2 = 7, k = k, p = p, xmin = xmin, xmax = xmax, 
                                                  nugget = nugget2, alpha = alpha, buffer = 0, 
                                                  candidates = x_grid)
}

gpvs12_noise = list("x_grid" = x_grid,
              "y_seq_mat" = y_seq_mat1d_expanded,
              "x_input_ind_mat" = x_input_ind_mat2,
              "y_input_mat" = y_input_mat2,
              "mmed_list" = mmed_gp2d_noise)
saveRDS(gpvs12_noise, file = "gpvs12_noise.rds")





# --- look at a particular simulation --- #
sim_index = 3
subdim = 1
mmed_gp_sim = gpvs12_noise
y_seq = mmed_gp_sim$y_seq_mat[ , sim_index]
x_train_ind = mmed_gp_sim$x_input_ind_mat[ , sim_index]
x_train = mmed_gp_sim$x_grid[x_train_ind, ]
y_train = mmed_gp_sim$y_input_mat[ , sim_index] # all.equal(y_train, y_seq[x_train_ind]) is TRUE
mmed_gp = gpvs12_noise$mmed_list[[sim_index]]
nugget = nugget2
# these graphs should be able to adjust according to specs above
# pre-processing for graphing
order_truef = order(x_grid[ , 1])
# plot wasserstein distance
initD = x_train
subinitD = initD[ , subdim, drop = FALSE]
if(is.null(nugget)){
  Kinv0 = solve(getCov(subinitD, subinitD, type01[1], l01[1]))
  Kinv1 = solve(getCov(initD, initD, type01[2], l01[2]))
} else{
  Kinv0 = solve(getCov(subinitD, subinitD, type01[1], l01[1]) + diag(rep(nugget, dim(initD)[1])))
  Kinv1 = solve(getCov(initD, initD, type01[2], l01[2]) + diag(rep(nugget, dim(initD)[1])))
}
w_seq = apply(x_grid, 1, FUN = function(x) Wasserstein_distance_postpred_gp2(x, Kinv0, Kinv1, subinitD, initD, y_train, var_e = 1, type01, l01))
quilt.plot(x_grid, w_seq, main = "Wasserstein(x)", xlab = "", ylab = "")
# plot input points in 2d
y_seq_hackcolor = y_seq
y_seq_hackcolor[x_train_ind] = 3
quilt.plot(x_grid, y_seq_hackcolor)
plot(x_grid[order_truef,1], y_seq[order_truef], type = "l")
points(x_grid[ x_train_ind2, 1], y_seq[ x_train_ind2])
# plot design points in 2d
y_seq_updated = y_seq
y_seq_updated[mmed_gp$indices] = 3
quilt.plot(x_grid, y_seq_updated)
# plot design points in 1d
plot(x_grid[order_truef,1], y_seq[order_truef], type = "l")
points(x_grid[ x_train_ind, 1], y_seq[ x_train_ind])
for(i in 1:7){
  points(mmed_gp$addD[i, 1], y_seq[mmed_gp$indices[i]], col = 2, pch = 3)
}
legend("topleft", legend = c("inputs", "MMED"), col = c(1, 2), pch = c(1, 3))
