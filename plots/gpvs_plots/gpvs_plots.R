################################################################################
### --- Gaussian Process Variable Selection -------------------------------- ###
################################################################################

# --- Working Directory --- #
home = "/home/kristyn/Documents/smed_ms"

# --- Sources/Libraries --- #
functions_home = paste(home, "/functions", sep="")
source(paste(functions_home, "/add_MMEDgp.R", sep = ""))
source(paste(functions_home, "/wasserstein_distance.R", sep = ""))
source(paste(functions_home, "/charge_function_q.R", sep = ""))
source(paste(functions_home, "/covariance_functions.R", sep = ""))
source(paste(functions_home, "/gp_predictive.R", sep = ""))
source(paste(functions_home, "/SMMEDgp.R", sep = ""))

library(expm)
library(matrixStats)
library(MASS)
library(mvtnorm)
library(fields)
library(knitr)

# for plots
library(ggplot2)
library(ggpubr)
library(reshape2)
library(data.table)
gg_color_hue = function(n) {
  hues = seq(15, 275, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
image_path = paste0(home, "/plots/gpvs_plots/gg")

# --- Helper Functions for Evaluation Metrics --- #
add_errorbands = function(xs, ys, MoE, color){
  y_lower = ys - MoE
  y_upper = ys + MoE
  polygon(c(xs,rev(xs)),c(y_lower,rev(y_upper)),col=color, border = NA)
}

getAICBIC = function(model = NULL, p, y, X, beta, sigmasq){
  # if(sum(class(model) == "glmnet") > 0){ # if glm
  #   tLL = model$nulldev - deviance(model)
  #   k = model$df
  #   n = model$nobs
  #   aic = -tLL + 2* k + 2* k * (k + 1) / (n - k - 1)
  #   bic = -tLL + log(n) * k
  # }
  # if(class(model) == "lm"){
  if(length(y) != dim(X)[1]) stop("length of y is not equal to number of rows in X")
  n = length(y)
  eval_criteria = function(crit_const){
    - 2 * dmvnorm(y, mean = X %*% beta, sigma = diag(rep(sigmasq, length(y)))) + crit_const * p
  }
  aic = eval_criteria(2)
  bic = eval_criteria(log(n))
  # }
  return(c("aic" = aic, "bic" = bic))
}

logjointlik = function(x, y, type_arg, l_arg, nugget = NULL){
  joint_var = getCov(x, x, type_arg, l_arg)
  if(!is.null(nugget)) joint_var = joint_var + diag(rep(nugget, length(y)))
  return(dmvnorm(y, mean = rep(0, length(y)), sigma = joint_var, log = TRUE))
}

logjointlik_star = function(x_star, y_star, x_input, y_input, type_arg, l_arg, nugget = NULL){
  y = c(y_input, y_star)
  x = rbind(x_input, x_star)
  logjointlik(x, y, type_arg, l_arg, nugget)
}

calcPPH = function(x, y, indices0, indices1, type_args, l_args, nugget){
  initD0.temp = x[ , indices0, drop = FALSE]
  initD1.temp = x[ , indices1, drop = FALSE]
  loglikH0 = logjointlik(initD0.temp, y, type_args[1], l_args[1], nugget)
  loglikH1 = logjointlik(initD1.temp, y, type_args[2], l_args[2], nugget)
  PPH0 = exp(loglikH0) / (exp(loglikH0) + exp(loglikH1))
  PPH1 = exp(loglikH1) / (exp(loglikH0) + exp(loglikH1))
  return(c("PPH0" = PPH0, "PPH1" = PPH1))
}

calcPPH_star = function(candidates, true_y, design_indices, inputx, inputy, 
                        indices0, indices1, hyp_type, hyp_l, nugget_term){
  loglikH0 = logjointlik_star(x_star = candidates[ design_indices, indices0, drop = FALSE], 
                              y_star = true_y[ design_indices], 
                              inputx[ , indices0, drop = FALSE], inputy, 
                              hyp_type[1], hyp_l[1], nugget = nugget_term)
  loglikH1 = logjointlik_star(x_star = candidates[ design_indices, indices1, drop = FALSE], 
                              y_star = true_y[ design_indices], 
                              inputx[ , indices1, drop = FALSE], inputy, 
                              hyp_type[2], hyp_l[2], nugget = nugget_term)
  PPH0 = exp(loglikH0) / (exp(loglikH0) + exp(loglikH1))
  PPH1 = exp(loglikH1) / (exp(loglikH0) + exp(loglikH1))
  return(c("PPH0" = PPH0, "PPH1" = PPH1))
}



### --- shared settings for both 1d f & 2d f cases --- ###

# --- simulations --- #
numSims = 25

# settings
xmin = 0; xmax = 1
l01= c(0.1, 0.1)
type01 = c(1, 1)

N = 3
p = 2
k = 4 * p
alpha = 1
indices0 = c(1)
indices1 = c(1, 2)
N2 = 9
nugget = 1e-1

# read in grid info
grid12 = readRDS(paste(home, "/plots/gpvs_plots/grid_gpvs12.rds", sep = ""))
# the grid of inputs
numx = grid12[[1]]
x_seq = grid12[[2]]
x_grid = as.matrix(grid12[[3]])
null_cov1d = grid12[[4]]
null_mean1d = grid12[[5]]
null_cov2d = grid12[[6]]
null_mean2d = grid12[[7]]

# input points (randomly selected)
set.seed(2)
# set.seed(6)
x_input_ind = sample(1:dim(x_grid)[1], N)
x_input = x_grid[x_input_ind, ]

order_x = order(x_grid[ , 1])

# read in simulations
gpvs_sims = readRDS(paste(home, "/run_designs/gp/gpsims_vs/run_designs_gpvs6/gpvs_sims.rds", sep = ""))
mmedgpf1d_sims = gpvs_sims$mmedgp_f1dsims
smmedgpf1d_sims = gpvs_sims$smmedgp_f1dsims
mmedgpf2d_sims = gpvs_sims$mmedgp_f2dsims
smmedgpf2d_sims = gpvs_sims$smmedgp_f2dsims

x_rand_ind_mat = matrix(NA, nrow = N2, ncol = numSims)
set.seed(1)
for(j in 1:numSims) x_rand_ind_mat[ , j] = sample(1:dim(x_grid)[1], N2)



# --- other designs to compare to mmed for gp --- #

# random points
set.seed(1990)
# set.seed(1997)
x_rand_ind = sample(1:dim(x_grid)[1], N2)
x_rand = x_grid[x_rand_ind, ]

# points where x2 = 1
gridlen = length(x_seq)
x_at1_ind = (1:N2) * floor(gridlen / N2) + gridlen * (gridlen - 1)
x_at1 = x_grid[x_at1_ind, ]

# points on diagonal
x_diag_ind = (1:N2) * floor(gridlen / N2) + gridlen * ((1:N2) * floor(gridlen / N2) - 1)
x_diag = x_grid[x_diag_ind, ]

# grid, only works when have 9 pts
beg_mid_end = c(1, ceiling(gridlen/(sqrt(9) - 1)) * 1:(sqrt(9) - 2), gridlen)
x_sf_ind = c(beg_mid_end, beg_mid_end + floor((gridlen * (gridlen - 1))) / 2, beg_mid_end + gridlen * (gridlen - 1))
x_sf = x_grid[x_sf_ind, ]



# --- 1d f simulations --- #

smmed_PPH0_mat = matrix(NA, nrow = N2 + 1, ncol = numSims)
mmed_PPH0_mat = matrix(NA, nrow = N2 + 1, ncol = numSims)
rand_PPH0_mat = matrix(NA, nrow = N2 + 1, ncol = numSims)
x_at1_PPH0_mat = matrix(NA, nrow = N2 + 1, ncol = numSims)
x_diag_PPH0_mat = matrix(NA, nrow = N2 + 1, ncol = numSims)
x_sf_PPH0_mat = matrix(NA, nrow = N2 + 1, ncol = numSims)

smmed_PPH1_mat = matrix(NA, nrow = N2 + 1, ncol = numSims)
mmed_PPH1_mat = matrix(NA, nrow = N2 + 1, ncol = numSims)
rand_PPH1_mat = matrix(NA, nrow = N2 + 1, ncol = numSims)
x_at1_PPH1_mat = matrix(NA, nrow = N2 + 1, ncol = numSims)
x_diag_PPH1_mat = matrix(NA, nrow = N2 + 1, ncol = numSims)
x_sf_PPH1_mat = matrix(NA, nrow = N2 + 1, ncol = numSims)

# for each new point in each design
set.seed(1)
for(j in 1:numSims){
  mmed_gp_vs = mmedgpf1d_sims[[j]]
  smmed_gp_vs = smmedgpf1d_sims[[j]]
  x_rand_ind_mat[ , j] = sample(1:dim(x_grid)[1], N2)
  x_rand_ind = x_rand_ind_mat[ , j]
  
  x_input_ind = gpvs_sims$x_input_ind[ , j]
  x_input = x_grid[x_input_ind, ]
  y_seq = gpvs_sims$y_seq_mat1d[ , j]
  y_input = y_seq[x_input_ind]
  
  # from initial data
  initPPH = calcPPH(x_input, y_input, indices0, indices1, type01, l01, nugget)
  
  smmed_PPH0_mat[1, j] = initPPH[1]
  mmed_PPH0_mat[1, j] = initPPH[1]
  rand_PPH0_mat[1, j] = initPPH[1]
  x_at1_PPH0_mat[1, j] = initPPH[1]
  x_diag_PPH0_mat[1, j] = initPPH[1]
  x_sf_PPH0_mat[1, j] = initPPH[1]
  
  smmed_PPH1_mat[1, j] = initPPH[2]
  mmed_PPH1_mat[1, j] = initPPH[2]
  rand_PPH1_mat[1, j] = initPPH[2]
  x_at1_PPH1_mat[1, j] = initPPH[2]
  x_diag_PPH1_mat[1, j] = initPPH[2]
  x_sf_PPH1_mat[1, j] = initPPH[2]
  
  for(i in 1:N2){
    which_ind = 1:i
    
    # smmed postprobs
    x_ind = smmed_gp_vs$indices[which_ind]
    smmedgpPPH.temp = calcPPH_star(x_grid, y_seq, x_ind, x_input, y_input,
                                   indices0, indices1, type01, l01, nugget)
    smmed_PPH0_mat[i + 1, j] = smmedgpPPH.temp[1]
    smmed_PPH1_mat[i + 1, j] = smmedgpPPH.temp[2]
    # mmed postprobs
    x_ind = mmed_gp_vs$indices[which_ind]
    mmedgpPPH.temp = calcPPH_star(x_grid, y_seq, x_ind, x_input, y_input,
                                  indices0, indices1, type01, l01, nugget)
    mmed_PPH0_mat[i + 1, j] = mmedgpPPH.temp[1]
    mmed_PPH1_mat[i + 1, j] = mmedgpPPH.temp[2]
    # random postprobs
    x_ind = x_rand_ind[which_ind]
    randPPH.temp = calcPPH_star(x_grid, y_seq, x_ind, x_input, y_input,
                                indices0, indices1, type01, l01, nugget)
    rand_PPH0_mat[i + 1, j] = randPPH.temp[1]
    rand_PPH1_mat[i + 1, j] = randPPH.temp[2]
    
    # x_at1 postprobs
    x_ind = x_at1_ind[which_ind]
    at1PPH.temp = calcPPH_star(x_grid, y_seq, x_ind, x_input, y_input,
                               indices0, indices1, type01, l01, nugget)
    x_at1_PPH0_mat[i + 1, j] = at1PPH.temp[1]
    x_at1_PPH1_mat[i + 1, j] = at1PPH.temp[2]
    
    # x_diag postprobs
    x_ind = x_diag_ind[which_ind]
    diagPPH.temp = calcPPH_star(x_grid, y_seq, x_ind, x_input, y_input,
                                indices0, indices1, type01, l01, nugget)
    x_diag_PPH0_mat[i + 1, j] = diagPPH.temp[1]
    x_diag_PPH1_mat[i + 1, j] = diagPPH.temp[2]
    
    # x_sf postprobs
    x_ind = x_sf_ind[which_ind]
    sfPPH.temp = calcPPH_star(x_grid, y_seq, x_ind, x_input, y_input,
                              indices0, indices1, type01, l01, nugget)
    x_sf_PPH0_mat[i + 1, j] = sfPPH.temp[1]
    x_sf_PPH1_mat[i + 1, j] = sfPPH.temp[2]
  }
}

par(mfrow = c(1, 2))

# by mean #
smmed_PPH0_avg = apply(smmed_PPH0_mat, 1, mean, na.rm = TRUE)
mmed_PPH0_avg = apply(mmed_PPH0_mat, 1, mean, na.rm = TRUE)
rand_PPH0_avg = apply(rand_PPH0_mat, 1, mean, na.rm = TRUE)
x_at1_PPH0_avg = apply(x_at1_PPH0_mat, 1, mean, na.rm = TRUE)
x_diag_PPH0_avg = apply(x_diag_PPH0_mat, 1, mean, na.rm = TRUE)
x_sf_PPH0_avg = apply(x_sf_PPH0_mat, 1, mean, na.rm = TRUE)

smmed_PPH1_avg = apply(smmed_PPH1_mat, 1, mean, na.rm = TRUE)
mmed_PPH1_avg = apply(mmed_PPH1_mat, 1, mean, na.rm = TRUE)
rand_PPH1_avg = apply(rand_PPH1_mat, 1, mean, na.rm = TRUE)
x_at1_PPH1_avg = apply(x_at1_PPH1_mat, 1, mean, na.rm = TRUE)
x_diag_PPH1_avg = apply(x_diag_PPH1_mat, 1, mean, na.rm = TRUE)
x_sf_PPH1_avg = apply(x_sf_PPH1_mat, 1, mean, na.rm = TRUE)

plot(x = 1:(N2 + 1), y = smmed_PPH0_avg, type = "l", 
     xlab = "# new points", ylab = "", #ylab = "E[P(H0|Y,X)|X]", 
     ylim = c(0, 1))
lines(1:(N2 + 1), rand_PPH0_avg, col = 2)
lines(1:(N2 + 1), x_at1_PPH0_avg, col = 3)
lines(1:(N2 + 1), x_diag_PPH0_avg, col = 4)
lines(1:(N2 + 1), x_sf_PPH0_avg, col = 5)
legend("bottomright", legend = c("smmed", "random", "x2=1", "diagonal", "grid"), 
       col = 1:5, lty = c(1, 1, 1, 1))
plot(x = 1:(N2 + 1), y = smmed_PPH1_avg, type = "l", 
     xlab = "# new points", ylab = "", #ylab = "E[P(H1|Y,X)|X]", 
     ylim = c(0, 1))
lines(1:(N2 + 1), rand_PPH1_avg, col = 2)
lines(1:(N2 + 1), x_at1_PPH1_avg, col = 3)
lines(1:(N2 + 1), x_sf_PPH1_avg, col = 5)

# by median #
smmed_PPH0_median = apply(smmed_PPH0_mat, 1, median, na.rm = TRUE)
mmed_PPH0_median = apply(mmed_PPH0_mat, 1, median, na.rm = TRUE)
rand_PPH0_median = apply(rand_PPH0_mat, 1, median, na.rm = TRUE)
x_at1_PPH0_median = apply(x_at1_PPH0_mat, 1, median, na.rm = TRUE)
x_diag_PPH0_median = apply(x_diag_PPH0_mat, 1, median, na.rm = TRUE)
x_sf_PPH0_median = apply(x_sf_PPH0_mat, 1, median, na.rm = TRUE)

smmed_PPH1_median = apply(smmed_PPH1_mat, 1, median, na.rm = TRUE)
mmed_PPH1_median = apply(mmed_PPH1_mat, 1, median, na.rm = TRUE)
rand_PPH1_median = apply(rand_PPH1_mat, 1, median, na.rm = TRUE)
x_at1_PPH1_median = apply(x_at1_PPH1_mat, 1, median, na.rm = TRUE)
x_diag_PPH1_median = apply(x_diag_PPH1_mat, 1, median, na.rm = TRUE)
x_sf_PPH1_median = apply(x_sf_PPH1_mat, 1, median, na.rm = TRUE)

plot(x = 1:(N2 + 1), y = smmed_PPH0_median, type = "l", 
     xlab = "# new points", ylab = "", #ylab = "E[P(H0|Y,X)|X]", 
     ylim = c(0, 1))
lines(1:(N2 + 1), rand_PPH0_median, col = 2)
lines(1:(N2 + 1), x_at1_PPH0_median, col = 3)
lines(1:(N2 + 1), x_diag_PPH0_median, col = 4)
lines(1:(N2 + 1), x_sf_PPH0_median, col = 5)
legend("bottomright", legend = c("smmed", "random", "x2=1", "diagonal", "grid"), 
       col = 1:5, lty = c(1, 1, 1, 1))
plot(x = 1:(N2 + 1), y = mmed_PPH1_median, type = "l", 
     xlab = "# new points", ylab = "", #ylab = "E[P(H1|Y,X)|X]", 
     ylim = c(0, 1))
lines(1:(N2 + 1), rand_PPH1_median, col = 2)
lines(1:(N2 + 1), x_at1_PPH1_median, col = 3)
lines(1:(N2 + 1), x_sf_PPH1_median, col = 5)


#
#
#

ggdata0 = data.table(
  x = 1:(N2 + 1), 
  Random = rand_PPH0_avg, 
  `X=1` = x_at1_PPH0_avg, 
  Diagonal = x_diag_PPH0_avg, 
  SpaceFilling = x_sf_PPH0_avg,
  SeqMED = smmed_PPH0_avg, 
  Hypothesis = rep("H0", idxlast)
)
ggdata1 = data.table(
  x = 1:(N2 + 1), 
  Random = rand_PPH1_avg, 
  `X=1` = x_at1_PPH1_avg, 
  Diagonal = x_diag_PPH1_avg, 
  SpaceFilling = x_sf_PPH1_avg,
  SeqMED = smmed_PPH1_avg, 
  Hypothesis = rep("H1", idxlast)
)
ggdata = rbind(ggdata0, ggdata1)
ggdata.melted = melt(ggdata, id = c("x", "Hypothesis"), value.name = "epph", 
                     variable.name = "Design")
plt = ggplot(ggdata.melted, aes(x = x, y = epph, color = Design, linetype = Design)) +
  facet_wrap(vars(Hypothesis)) + 
  geom_path() + 
  scale_linetype_manual(values=c(rep("dashed", 4), "solid")) + 
  geom_point(data = ggdata.melted[x == 10], aes(x = x, y = epph)) + 
  theme_bw() + 
  theme(panel.grid.minor = element_blank()) + 
  labs(y = "", x = "Stages")
plt
ggsave("h0_epph.pdf",
       plot = last_plot(),
       device = "pdf",
       path = image_path,
       scale = 1,
       width = 8,
       height = 4,
       units = c("in")
)

ggdata0.melted = melt(ggdata0[, 1:5], id = c("x"), value.name = "epph", 
                      variable.name = "Design")
plt2.0 = ggplot(ggdata0.melted, aes(x = x, y = epph, color = Design, linetype = Design)) +
  geom_path() + 
  scale_linetype_manual(values=c(rep("dashed", 3), "solid")) + 
  geom_point(data = ggdata0.melted[x == 10], aes(x = x, y = epph)) + 
  theme_bw() + 
  theme(panel.grid.minor = element_blank()) + 
  labs(y = "", x = "Stages")
plt2.0
# ggsave("h0_epph0_seq_h3.png",
#        plot = last_plot(),
#        device = "png",
#        path = image_path,
#        scale = 1,
#        width = 4,
#        height = 3,
#        units = c("in")
# )









# --- 2d f simulations --- #

smmed_PPH0_mat = matrix(NA, nrow = N2 + 1, ncol = numSims)
mmed_PPH0_mat = matrix(NA, nrow = N2 + 1, ncol = numSims)
rand_PPH0_mat = matrix(NA, nrow = N2 + 1, ncol = numSims)
x_at1_PPH0_mat = matrix(NA, nrow = N2 + 1, ncol = numSims)
x_diag_PPH0_mat = matrix(NA, nrow = N2 + 1, ncol = numSims)
x_sf_PPH0_mat = matrix(NA, nrow = N2 + 1, ncol = numSims)

smmed_PPH1_mat = matrix(NA, nrow = N2 + 1, ncol = numSims)
mmed_PPH1_mat = matrix(NA, nrow = N2 + 1, ncol = numSims)
rand_PPH1_mat = matrix(NA, nrow = N2 + 1, ncol = numSims)
x_at1_PPH1_mat = matrix(NA, nrow = N2 + 1, ncol = numSims)
x_diag_PPH1_mat = matrix(NA, nrow = N2 + 1, ncol = numSims)
x_sf_PPH1_mat = matrix(NA, nrow = N2 + 1, ncol = numSims)

# for each new point in each design
set.seed(1)
for(j in 1:numSims){
  mmed_gp_vs = mmedgpf2d_sims[[j]]
  smmed_gp_vs = smmedgpf2d_sims[[j]]
  x_rand_ind = x_rand_ind_mat[ , j] # created in the beginning
  
  x_input_ind = gpvs_sims$x_input_ind[ , j]
  x_input = x_grid[x_input_ind, ]
  y_seq = gpvs_sims$y_seq_mat2d[ , j]
  y_input = y_seq[x_input_ind]
  
  # from initial data
  initPPH = calcPPH(x_input, y_input, indices0, indices1, type01, l01, nugget)
  
  smmed_PPH0_mat[1, j] = initPPH[1]
  mmed_PPH0_mat[1, j] = initPPH[1]
  rand_PPH0_mat[1, j] = initPPH[1]
  x_at1_PPH0_mat[1, j] = initPPH[1]
  x_diag_PPH0_mat[1, j] = initPPH[1]
  x_sf_PPH0_mat[1, j] = initPPH[1]
  
  smmed_PPH1_mat[1, j] = initPPH[2]
  mmed_PPH1_mat[1, j] = initPPH[2]
  rand_PPH1_mat[1, j] = initPPH[2]
  x_at1_PPH1_mat[1, j] = initPPH[2]
  x_diag_PPH1_mat[1, j] = initPPH[2]
  x_sf_PPH1_mat[1, j] = initPPH[2]
  
  for(i in 1:N2){
    which_ind = 1:i
    
    # smmed postprobs
    x_ind = smmed_gp_vs$indices[which_ind]
    smmedgpPPH.temp = calcPPH_star(x_grid, y_seq, x_ind, x_input, y_input,
                                   indices0, indices1, type01, l01, nugget)
    smmed_PPH0_mat[i + 1, j] = smmedgpPPH.temp[1]
    smmed_PPH1_mat[i + 1, j] = smmedgpPPH.temp[2]
    
    # mmed postprobs
    x_ind = mmed_gp_vs$indices[which_ind]
    mmedgpPPH.temp = calcPPH_star(x_grid, y_seq, x_ind, x_input, y_input,
                                  indices0, indices1, type01, l01, nugget)
    mmed_PPH0_mat[i + 1, j] = mmedgpPPH.temp[1]
    mmed_PPH1_mat[i + 1, j] = mmedgpPPH.temp[2]
    
    # random postprobs
    x_ind = x_rand_ind[which_ind]
    randPPH.temp = calcPPH_star(x_grid, y_seq, x_ind, x_input, y_input,
                                indices0, indices1, type01, l01, nugget)
    rand_PPH0_mat[i + 1, j] = randPPH.temp[1]
    rand_PPH1_mat[i + 1, j] = randPPH.temp[2]
    
    # x_at1 postprobs
    x_ind = x_at1_ind[which_ind]
    at1PPH.temp = calcPPH_star(x_grid, y_seq, x_ind, x_input, y_input,
                               indices0, indices1, type01, l01, nugget)
    x_at1_PPH0_mat[i + 1, j] = at1PPH.temp[1]
    x_at1_PPH1_mat[i + 1, j] = at1PPH.temp[2]
    
    # x_diag postprobs
    x_ind = x_diag_ind[which_ind]
    diagPPH.temp = calcPPH_star(x_grid, y_seq, x_ind, x_input, y_input,
                                indices0, indices1, type01, l01, nugget)
    x_diag_PPH0_mat[i + 1, j] = diagPPH.temp[1]
    x_diag_PPH1_mat[i + 1, j] = diagPPH.temp[2]
    
    # x_sf postprobs
    x_ind = x_sf_ind[which_ind]
    sfPPH.temp = calcPPH_star(x_grid, y_seq, x_ind, x_input, y_input,
                              indices0, indices1, type01, l01, nugget)
    x_sf_PPH0_mat[i + 1, j] = sfPPH.temp[1]
    x_sf_PPH1_mat[i + 1, j] = sfPPH.temp[2]
  }
}

par(mfrow = c(1, 2))

# by mean #
smmed_PPH0_avg = apply(smmed_PPH0_mat, 1, mean, na.rm = TRUE)
mmed_PPH0_avg = apply(mmed_PPH0_mat, 1, mean, na.rm = TRUE)
rand_PPH0_avg = apply(rand_PPH0_mat, 1, mean, na.rm = TRUE)
x_at1_PPH0_avg = apply(x_at1_PPH0_mat, 1, mean, na.rm = TRUE)
x_diag_PPH0_avg = apply(x_diag_PPH0_mat, 1, mean, na.rm = TRUE)
x_sf_PPH0_avg = apply(x_sf_PPH0_mat, 1, mean, na.rm = TRUE)

smmed_PPH1_avg = apply(smmed_PPH1_mat, 1, mean, na.rm = TRUE)
mmed_PPH1_avg = apply(mmed_PPH1_mat, 1, mean, na.rm = TRUE)
rand_PPH1_avg = apply(rand_PPH1_mat, 1, mean, na.rm = TRUE)
x_at1_PPH1_avg = apply(x_at1_PPH1_mat, 1, mean, na.rm = TRUE)
x_diag_PPH1_avg = apply(x_diag_PPH1_mat, 1, mean, na.rm = TRUE)
x_sf_PPH1_avg = apply(x_sf_PPH1_mat, 1, mean, na.rm = TRUE)

plot(x = 1:(N2 + 1), y = smmed_PPH0_avg, type = "l", 
     xlab = "# new points", ylab = "", #ylab = "E[P(H0|Y,X)|X]", 
     ylim = c(0, 1))
lines(1:(N2 + 1), rand_PPH0_avg, col = 2)
lines(1:(N2 + 1), x_at1_PPH0_avg, col = 3)
lines(1:(N2 + 1), x_diag_PPH0_avg, col = 4)
lines(1:(N2 + 1), x_sf_PPH0_avg, col = 5)
legend("topright", legend = c("smmed", "random", "x2=1", "diagonal", "grid"), 
       col = 1:5, lty = c(1, 1, 1, 1))
plot(x = 1:(N2 + 1), y = smmed_PPH1_avg, type = "l", 
     xlab = "# new points", ylab = "", #ylab = "E[P(H1|Y,X)|X]", 
     ylim = c(0, 1))
lines(1:(N2 + 1), rand_PPH1_avg, col = 2)
lines(1:(N2 + 1), x_at1_PPH1_avg, col = 3)
lines(1:(N2 + 1), x_sf_PPH1_avg, col = 5)

# by median #
smmed_PPH0_median = apply(smmed_PPH0_mat, 1, median, na.rm = TRUE)
mmed_PPH0_median = apply(mmed_PPH0_mat, 1, median, na.rm = TRUE)
rand_PPH0_median = apply(rand_PPH0_mat, 1, median, na.rm = TRUE)
x_at1_PPH0_median = apply(x_at1_PPH0_mat, 1, median, na.rm = TRUE)
x_diag_PPH0_median = apply(x_diag_PPH0_mat, 1, median, na.rm = TRUE)
x_sf_PPH0_median = apply(x_sf_PPH0_mat, 1, median, na.rm = TRUE)

smmed_PPH1_median = apply(smmed_PPH1_mat, 1, median, na.rm = TRUE)
mmed_PPH1_median = apply(mmed_PPH1_mat, 1, median, na.rm = TRUE)
rand_PPH1_median = apply(rand_PPH1_mat, 1, median, na.rm = TRUE)
x_at1_PPH1_median = apply(x_at1_PPH1_mat, 1, median, na.rm = TRUE)
x_diag_PPH1_median = apply(x_diag_PPH1_mat, 1, median, na.rm = TRUE)
x_sf_PPH1_median = apply(x_sf_PPH1_mat, 1, median, na.rm = TRUE)

plot(x = 1:(N2 + 1), y = smmed_PPH0_median, type = "l", 
     xlab = "# new points", ylab = "", #ylab = "E[P(H0|Y,X)|X]", 
     ylim = c(0, 1))
lines(1:(N2 + 1), rand_PPH0_median, col = 2)
lines(1:(N2 + 1), x_at1_PPH0_median, col = 3)
lines(1:(N2 + 1), x_diag_PPH0_median, col = 4)
lines(1:(N2 + 1), x_sf_PPH0_median, col = 5)
legend("topright", legend = c("smmed", "random", "x2=1", "diagonal", "grid"), 
       col = 1:5, lty = c(1, 1, 1, 1))
plot(x = 1:(N2 + 1), y = smmed_PPH1_median, type = "l", 
     xlab = "# new points", ylab = "", #ylab = "E[P(H1|Y,X)|X]", 
     ylim = c(0, 1))
lines(1:(N2 + 1), rand_PPH1_median, col = 2)
lines(1:(N2 + 1), x_at1_PPH1_median, col = 3)
lines(1:(N2 + 1), x_sf_PPH1_median, col = 5)




#
#
#

ggdata0 = data.table(
  x = 1:(N2 + 1), 
  Random = rand_PPH0_avg, 
  `X=1` = x_at1_PPH0_avg, 
  Diagonal = x_diag_PPH0_avg, 
  SpaceFilling = x_sf_PPH0_avg,
  SeqMED = smmed_PPH0_avg, 
  Hypothesis = rep("H0", idxlast)
)
ggdata1 = data.table(
  x = 1:(N2 + 1), 
  Random = rand_PPH1_avg, 
  `X=1` = x_at1_PPH1_avg, 
  Diagonal = x_diag_PPH1_avg, 
  SpaceFilling = x_sf_PPH1_avg,
  SeqMED = smmed_PPH1_avg, 
  Hypothesis = rep("H1", idxlast)
)
ggdata = rbind(ggdata0, ggdata1)
ggdata.melted = melt(ggdata, id = c("x", "Hypothesis"), value.name = "epph", 
                     variable.name = "Design")
plt2 = ggplot(ggdata.melted, aes(x = x, y = epph, color = Design, linetype = Design)) +
  facet_wrap(vars(Hypothesis)) + 
  geom_path() + 
  scale_linetype_manual(values=c(rep("dashed", 4), "solid")) + 
  geom_point(data = ggdata.melted[x == 10], aes(x = x, y = epph)) + 
  theme_bw() + 
  theme(panel.grid.minor = element_blank()) + 
  labs(y = "", x = "Stages")
plt2
ggsave("h1_epph.pdf",
       plot = last_plot(),
       device = "pdf",
       path = image_path,
       scale = 1,
       width = 8,
       height = 4,
       units = c("in")
)

ggdata1.melted = melt(ggdata1[, 1:5], id = c("x"), value.name = "epph", 
                      variable.name = "Design")
plt2.1 = ggplot(ggdata1.melted, aes(x = x, y = epph, color = Design, linetype = Design)) +
  geom_path() + 
  scale_linetype_manual(values=c(rep("dashed", 3), "solid")) + 
  geom_point(data = ggdata1.melted[x == 10], aes(x = x, y = epph)) + 
  theme_bw() + 
  theme(panel.grid.minor = element_blank()) + 
  labs(y = "", x = "Stages")
plt2.1
# ggsave("h1_epph1_seq_h3.png",
#        plot = last_plot(),
#        device = "png",
#        path = image_path,
#        scale = 1,
#        width = 4,
#        height = 3,
#        units = c("in")
# )
