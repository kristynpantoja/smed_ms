################################################################################
### --- Gaussian Process Regression ---------------------------------------- ###
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
image_path = paste0(home, "/plots/gp_plots/gg")

# --- Helper Functions for Evaluation Metrics --- #
add_errorbands = function(xs, ys, MoE, color){
  y_lower = ys - MoE
  y_upper = ys + MoE
  polygon(c(xs,rev(xs)),c(y_lower,rev(y_upper)),col=color, border = NA)
}

sigmoid = function(x) {
  1 / (1 + exp(-x))
}

getMedianLogRSS01 = function(simcase, design = NULL){
  if(is.null(design)) stop("must specify design! 1 = mmed, 2 = spacefilling, 3 = uniform")
  if(design == 1) logRSS01_vec = (simcase$RSS01mmed)
  if(design == 2) logRSS01_vec = (simcase$RSS01sf)
  if(design == 3) logRSS01_vec = (simcase$RSS01unif)
  return(median(logRSS01_vec, na.rm = TRUE))
}

getMedianLogPredEvid01 = function(simcase, design = NULL){
  if(is.null(design)) stop("must specify design! 1 = mmed, 2 = spacefilling, 3 = uniform")
  if(design == 1) logPredEvid01_vec = simcase$logPredEvid01mmed
  if(design == 2) logPredEvid01_vec = simcase$logPredEvid01sf
  if(design == 3) logPredEvid01_vec = simcase$logPredEvid01unif
  return(median(logPredEvid01_vec, na.rm = TRUE))
}

getMedianPostProbH1 = function(simcase, design = NULL){
  if(is.null(design)) stop("must specify design! 1 = mmed, 2 = spacefilling, 3 = uniform")
  if(design == 1) postProbH1_vec = exp(simcase$logJointEvid1mmed_vec) / (exp(simcase$logJointEvid0mmed_vec) + exp(simcase$logJointEvid1mmed_vec))
  if(design == 2) postProbH1_vec = exp(simcase$logJointEvid1sf_vec) / (exp(simcase$logJointEvid0sf_vec) + exp(simcase$logJointEvid1sf_vec))
  if(design == 3) postProbH1_vec = exp(simcase$logJointEvid1unif_vec) / (exp(simcase$logJointEvid0unif_vec) + exp(simcase$logJointEvid1unif_vec))
  return(median(postProbH1_vec, na.rm = TRUE))
}


### --- shared settings for both 1d f & 2d f cases --- ###

# --- simulations --- #

# x_seq, grid over which to generate subsequent functions
xmin = 0; xmax = 1
numx = 1001
x_seq = seq(from = xmin, to = xmax, length.out = numx) # set training points



# train sets

# train set designs
N = 6
N2 = 15
Ntotal = N + N2

# 1. make space-filling design
# space_filling = seq(from = xmin, to = xmax, length.out = Ntotal)
space_filling_ind = c(1, 1 + ((numx - 1)/(Ntotal - 1)) * 1:((numx - 1) / ((numx - 1)/(Ntotal - 1))))
space_filling = x_seq[space_filling_ind]

# train set 1
x_train1_ind = space_filling_ind[1:N]
x_train1 = x_seq[x_train1_ind]
x_spacefill1_ind = space_filling_ind[-c(1:N)]
x_spacefill1 = x_seq[x_spacefill1_ind]
# all.equal(space_filling, c(x_train1, x_spacefill1))

# train set 2
x_train2_ind = space_filling_ind[c(1, 2, 4, 7, 12, 21)]
x_train2 = x_seq[x_train2_ind]
x_spacefill2_ind = space_filling_ind[-c(1, 2, 4, 7, 12, 21)]
x_spacefill2 = x_seq[x_spacefill2_ind]
# all.equal(space_filling, sort(c(x_train2, x_spacefill2)))

# train set 3 (space-filling)
x_train3_ind = c(1, 1 + ((numx - 1)/(N - 1)) * 1:((numx - 1) / ((numx - 1)/(N - 1))))
x_train3 = x_seq[x_train3_ind]
x_spacefill3_ind = space_filling_ind[!(space_filling_ind %in% x_train3_ind)]
x_spacefill3 = x_seq[x_spacefill3_ind]
# all.equal(space_filling, sort(c(x_train3, x_spacefill3)))



# --- Matern vs Gaussian --- #
l01= c(0.1, 0.1)
type01 = c(1, 4)
numCandidates = 1001
k = 4
p = 1
nugget = NULL
alpha = 1

# generate matern function
set.seed(12)
null_cov = getCov(x_seq, x_seq, type01[2], l01[2])
null_mean = rep(0, numx)

# read in sims
gaussianvsmatern_train4sims = readRDS(paste(home, "/run_designs/gp/gpsims/gaussianvsmatern_train4sims.rds", sep = ""))

# --- plots sims --- #

par(mfrow = c(2, 2))

# get sim info
sim_ind = 1
x_train = gaussianvsmatern_train4sims$x_train[ , sim_ind]
x_train_ind = gaussianvsmatern_train4sims$x_train_ind[ , sim_ind]
mmed_gp = gaussianvsmatern_train4sims$mmed_gp_list[[sim_ind]]
y_seq = gaussianvsmatern_train4sims$sim_fns[ , sim_ind]
y_train = y_seq[x_train_ind]

# plot sf
x_spacefill_ind = gaussianvsmatern_train4sims$x_spacefill_ind
x_spacefill = gaussianvsmatern_train4sims$x_spacefill
newpts = x_spacefill
truey = y_seq[x_spacefill_ind]

H0_pred = getGPPredictive(newpts, x_train, y_train, type01[1], l01[1], nugget = NULL)
H1_pred = getGPPredictive(newpts, x_train, y_train, type01[2], l01[2], nugget = NULL)
postpredmu0 = H0_pred$pred_mean
postpredmu1 = H1_pred$pred_mean
plot(x_seq, y_seq, type = "l", ylim = c(-2.5, 2.5), ylab = "y", xlab = "x",
     main = "Space-filling") # plot the function
H0_predfn = getGPPredictive(x_seq, x_train, y_train, type01[1], l01[1], nugget = NULL)
H1_predfn = getGPPredictive(x_seq, x_train, y_train, type01[2], l01[2], nugget = NULL)
lines(x = x_seq, y = H0_predfn$pred_mean, col = 3, lty = 3)
lines(x = x_seq, y = H1_predfn$pred_mean, col = 4, lty = 2)
err0 = 2 * sqrt(diag(H0_predfn$pred_var))
err1 = 2 * sqrt(diag(H1_predfn$pred_var))
add_errorbands(x_seq, H0_predfn$pred_mean, err0, rgb(0, 1, 0, 0.15))
add_errorbands(x_seq, H1_predfn$pred_mean, err1, rgb(0, 0, 1, 0.15))
points(x = x_train, y = y_train, pch = 16)
points(x = x_train, y = rep(-2.5, N), col = rgb(0, 0, 0, 0.5), pch = 16)
points(x = newpts, y = truey, col = 2, pch = 16)
points(x = newpts, y = rep(-2.5, N2), col = rgb(1, 0, 0, 0.5), pch = 16)
points(x = newpts, y = postpredmu0, col = 3, pch = 16)
points(x = newpts, y = postpredmu1, col = 4, pch = 16)
# legend("bottomleft", legend = c("train y", "true y", "H0:gaussian", "H1:matern"),
#        col = c(1:4), pch = rep(16,4))

# plot mmed
newpts = mmed_gp$addD
truey = y_seq[mmed_gp$indices]

H0_pred = getGPPredictive(newpts, x_train, y_train, type01[1], l01[1], 
                          nugget = NULL)
H1_pred = getGPPredictive(newpts, x_train, y_train, type01[2], l01[2], 
                          nugget = NULL)
postpredmu0 = H0_pred$pred_mean
postpredmu1 = H1_pred$pred_mean
plot(x_seq, y_seq, type = "l", ylim = c(-2.5, 2.5), ylab = "y", xlab = "x",
     main = "SeqMED") # plot the function
H0_predfn = getGPPredictive(x_seq, x_train, y_train, type01[1], l01[1],
                            nugget = NULL)
H1_predfn = getGPPredictive(x_seq, x_train, y_train, type01[2], l01[2],
                            nugget = NULL)
lines(x = x_seq, y = H0_predfn$pred_mean, col = 3, lty = 3)
lines(x = x_seq, y = H1_predfn$pred_mean, col = 4, lty = 2)
err0 = 2 * sqrt(diag(H0_predfn$pred_var))
err1 = 2 * sqrt(diag(H1_predfn$pred_var))
add_errorbands(x_seq, H0_predfn$pred_mean, err0, rgb(0, 1, 0, 0.15))
add_errorbands(x_seq, H1_predfn$pred_mean, err1, rgb(0, 0, 1, 0.15))
points(x = x_train, y = y_train, pch = 16)
points(x = x_train, y = rep(-2.5, N), col = rgb(0, 0, 0, 0.5), pch = 16)
points(x = newpts, y = truey, col = 2, pch = 16)
points(x = newpts, y = rep(-2.5, N2), col = rgb(1, 0, 0, 0.5), pch = 16)
points(x = newpts, y = postpredmu0, col = 3, pch = 16)
points(x = newpts, y = postpredmu1, col = 4, pch = 16)
# legend("bottomleft", legend = c("train y", "true y", 
#                                 "H0:gaussian", "H1:matern"),
#        col = c(1:4), pch = rep(16,4))


hist(mmed_gp$addD, main = "SeqMED", xlab = "design points", breaks = 10)

# get w_seq
Kinv0 = solve(getCov(x_train, x_train, type01[1], l01[1]))
Kinv1 = solve(getCov(x_train, x_train, type01[2], l01[2]))
w_seq = sapply(x_seq, FUN = function(x1) 
  Wasserstein_distance_postpred_gp(x1, Kinv0, Kinv1, x_train, y_train, 
                                   var_e = 1, type01, l01))
plot(x_seq, w_seq, type = "l", ylim = c(0, 0.45), main = "Wasserstein(x)", 
     xlab = "", ylab = "")
points(mmed_gp$addD, w_seq[mmed_gp$indices], col = rgb(1, 0, 0, 1), pch = 16)
text(mmed_gp$addD, w_seq[mmed_gp$indices], labels = 1:N2, pos = 3, 
     offset = 0.1, cex = 1, col = 4)


#
#
#
err0 = 2 * sqrt(diag(H0_predfn$pred_var))
err1 = 2 * sqrt(diag(H1_predfn$pred_var))
ggdata = data.table(
  x = x_seq, 
  `True Function` = y_seq, 
  Wasserstein = w_seq, 
  `H0 Predictive Mean` = H0_predfn$pred_mean, 
  `H1 Predictive Mean` = H1_predfn$pred_mean,
  lower0 = H0_predfn$pred_mean - err0, 
  lower1 = H1_predfn$pred_mean - err1, 
  upper0 = H0_predfn$pred_mean + err0,
  upper1 = H1_predfn$pred_mean + err1
)
yrange = range(ggdata$lower0, ggdata$lower1, 
               ggdata$upper0, ggdata$upper1)
yrange[1] = yrange[1] - 1
ggdata$Wasserstein = ggdata$Wasserstein - abs(yrange[1])
ggdata$zero1 = NA
ggdata$zero2 = NA
ggdata.melted = melt(ggdata, id.vars = c("x"), 
                     measure.vars = c("True Function", "Wasserstein", "H0 Predictive Mean", "H1 Predictive Mean"))
ggdata.lower = melt(ggdata, id.vars = c("x"), 
                    measure.vars = c("zero1", "zero2", "lower0", "lower1"))
ggdata.upper = melt(ggdata, id.vars = c("x"), 
                    measure.vars = c("zero1", "zero2", "upper0", "upper1"))
ggdata.melted = cbind(ggdata.melted, 
                      lower = ggdata.lower$value, 
                      upper = ggdata.upper$value)
ggdata_pts = data.table(
  x = c(x_train, newpts), 
  y = c(y_train, truey), 
  color = c(rep(gg_color_hue(2)[2], length(x_train)), 
            rep(gg_color_hue(2)[1], length(newpts))), 
  shape = c(rep(8, length(x_train)), 
            rep(16, length(newpts)))
)
ggplot(data = ggdata.melted, aes(x = x, y =value, color = variable), 
       linetype = 1) + 
  geom_path(size = 2) + 
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = variable), 
              alpha = 0.1, linetype = 0) +
  scale_linetype_manual(values = c(1, 1, 2, 2)) + 
  scale_fill_manual(values = c(NA, NA, "#00BFC4", "#C77CFF")) + 
  scale_color_manual(values = c(1, "gray", "#00BFC4", "#C77CFF")) + 
  geom_point(data = ggdata_pts, mapping = aes(x = x, y = y), 
             inherit.aes = FALSE, color = ggdata_pts$color, 
             shape = ggdata_pts$shape,
             size = 3) +
  geom_point(data = ggdata_pts, mapping = aes(x = x, y = yrange[1]), 
             inherit.aes = FALSE, color = ggdata_pts$color, 
             shape = ggdata_pts$shape, 
             size = 3) +
  scale_y_continuous(limits = yrange) +
  theme_bw(base_size = 20) + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  labs(y = "y", x = "x", fill = "Function", color = "Function")
# ggsave("poster_gvm_h8.png",
#        plot = last_plot(),
#        device = "png",
#        path = image_path,
#        scale = 1,
#        width = 13.5,
#        height = 8,
#        units = c("in")
# )
ggplot(data = ggdata.melted, aes(x = x, y =value, color = variable), 
       linetype = 1) + 
  geom_path() + 
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = variable), 
              alpha = 0.1, linetype = 0) +
  scale_linetype_manual(values = c(1, 1, 2, 2)) + 
  scale_fill_manual(values = c(NA, NA, "#00BFC4", "#C77CFF")) + 
  scale_color_manual(values = c(1, "gray", "#00BFC4", "#C77CFF")) + 
  geom_point(data = ggdata_pts, mapping = aes(x = x, y = y), 
             inherit.aes = FALSE, color = ggdata_pts$color, 
             shape = ggdata_pts$shape, 
             size = 3) +
  geom_point(data = ggdata_pts, mapping = aes(x = x, y = yrange[1]), 
             inherit.aes = FALSE, color = ggdata_pts$color, 
             shape = ggdata_pts$shape, 
             size = 3) +
  scale_y_continuous(limits = yrange) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  labs(y = "y", x = "x", fill = "Function", color = "Function")
# ggsave("gvm.pdf",
#        plot = last_plot(),
#        device = "pdf",
#        path = image_path,
#        scale = 1,
#        width = 8,
#        height = 4,
#        units = c("in")
# )






# --- Matern vs Periodic --- #
# MMED parameters for testing
l01= c(0.1, 0.5); type01 = c(4, 5); numCandidates = 1001; k = 4; p = 1; nugget = NULL; alpha = 1

# generate periodic function
set.seed(13)
null_cov = getCov(x_seq, x_seq, type01[2], l01[2])
null_mean = rep(0, numx)

# read in sims

maternvsperiodic_train4sims = readRDS(paste(home, "/run_designs/gp/gpsims/maternvsperiodic_train4sims.rds", sep = ""))

# --- plot sims --- #

par(mfrow = c(2, 2))

# get sim info
sim_ind = 1
x_train = maternvsperiodic_train4sims$x_train[ , sim_ind]
x_train_ind = maternvsperiodic_train4sims$x_train_ind[ , sim_ind]
mmed_gp = maternvsperiodic_train4sims$mmed_gp_list[[sim_ind]]
y_seq = maternvsperiodic_train4sims$sim_fns[ , sim_ind]
y_train = y_seq[x_train_ind]

# plot sf
x_spacefill_ind = maternvsperiodic_train4sims$x_spacefill_ind
x_spacefill = maternvsperiodic_train4sims$x_spacefill
newpts = x_spacefill
truey = y_seq[x_spacefill_ind]

H0_pred = getGPPredictive(newpts, x_train, y_train, type01[1], l01[1], nugget = NULL)
H1_pred = getGPPredictive(newpts, x_train, y_train, type01[2], l01[2], nugget = NULL)
postpredmu0 = H0_pred$pred_mean
postpredmu1 = H1_pred$pred_mean
plot(x_seq, y_seq, type = "l", ylim = c(-2.5, 2.5), ylab = "y", xlab = "x",
     main = "Space-filling") # plot the function
H0_predfn = getGPPredictive(x_seq, x_train, y_train, type01[1], l01[1], nugget = NULL)
H1_predfn = getGPPredictive(x_seq, x_train, y_train, type01[2], l01[2], nugget = NULL)
lines(x = x_seq, y = H0_predfn$pred_mean, col = 3, lty = 3)
lines(x = x_seq, y = H1_predfn$pred_mean, col = 4, lty = 2)
err0 = 2 * sqrt(diag(H0_predfn$pred_var))
err1 = 2 * sqrt(diag(H1_predfn$pred_var))
add_errorbands(x_seq, H0_predfn$pred_mean, err0, rgb(0, 1, 0, 0.2))
add_errorbands(x_seq, H1_predfn$pred_mean, err1, rgb(0, 0, 1, 0.2))
points(x = x_train, y = y_train, pch = 16)
points(x = x_train, y = rep(-2.5, N), col = rgb(0, 0, 0, 0.5), pch = 16)
points(x = newpts, y = truey, col = 2, pch = 16)
points(x = newpts, y = rep(-2.5, N2), col = rgb(1, 0, 0, 0.5), pch = 16)
points(x = newpts, y = postpredmu0, col = 3, pch = 16)
points(x = newpts, y = postpredmu1, col = 4, pch = 16)
# legend("bottomright", legend = c("train y", "true y", "H0:matern", "H1:periodic"),
#        col = c(1:4), pch = rep(16,4))


# plot mmed
newpts = mmed_gp$addD
truey = y_seq[mmed_gp$indices]

H0_pred = getGPPredictive(newpts, x_train, y_train, type01[1], l01[1], nugget = NULL)
H1_pred = getGPPredictive(newpts, x_train, y_train, type01[2], l01[2], nugget = NULL)
postpredmu0 = H0_pred$pred_mean
postpredmu1 = H1_pred$pred_mean
plot(x_seq, y_seq, type = "l", ylim = c(-2.5, 2.5), ylab = "y", xlab = "x", 
     main = "SeqMED") # plot the function
H0_predfn = getGPPredictive(x_seq, x_train, y_train, type01[1], l01[1], nugget = NULL)
H1_predfn = getGPPredictive(x_seq, x_train, y_train, type01[2], l01[2], nugget = NULL)
lines(x = x_seq, y = H0_predfn$pred_mean, col = 3, lty = 3)
lines(x = x_seq, y = H1_predfn$pred_mean, col = 4, lty = 2)
err0 = 2 * sqrt(diag(H0_predfn$pred_var))
err1 = 2 * sqrt(diag(H1_predfn$pred_var))
add_errorbands(x_seq, H0_predfn$pred_mean, err0, rgb(0, 1, 0, 0.2))
add_errorbands(x_seq, H1_predfn$pred_mean, err1, rgb(0, 0, 1, 0.2))
points(x = x_train, y = y_train, pch = 16)
points(x = x_train, y = rep(-2.5, N), col = rgb(0, 0, 0, 0.5), pch = 16)
points(x = newpts, y = truey, col = 2, pch = 16)
points(x = newpts, y = rep(-2.5, N2), col = rgb(1, 0, 0, 0.5), pch = 16)
points(x = newpts, y = postpredmu0, col = 3, pch = 16)
points(x = newpts, y = postpredmu1, col = 4, pch = 16)
# legend("bottomright", legend = c("train y", "true y", "H0:matern", "H1:periodic"),
#        col = c(1:4), pch = rep(16,4))

hist(mmed_gp$addD, main = "M-MED", xlab = "design points", breaks = 10)

# get w_seq
Kinv0 = solve(getCov(x_train, x_train, type01[1], l01[1]))
Kinv1 = solve(getCov(x_train, x_train, type01[2], l01[2]))
w_seq = sapply(x_seq, FUN = function(x1) Wasserstein_distance_postpred_gp(x1, Kinv0, Kinv1, x_train, y_train, var_e = 1, type01, l01))
plot(x_seq, w_seq, type = "l", ylim = c(0, 3.3), main = "Wasserstein(x)", xlab = "", ylab = "")
points(mmed_gp$addD, w_seq[mmed_gp$indices], col = rgb(1, 0, 0, 1), pch = 16)
text(mmed_gp$addD, w_seq[mmed_gp$indices], labels = 1:N2, pos = 3, offset = 0.1, cex = 1, col = 4)


#
#
#
err0 = 2 * sqrt(diag(H0_predfn$pred_var))
err1 = 2 * sqrt(diag(H1_predfn$pred_var))
ggdata = data.table(
  x = x_seq, 
  `True Function` = y_seq, 
  Wasserstein = w_seq, 
  `H0 Pred Mean` = H0_predfn$pred_mean, 
  `H1 Pred Mean` = H1_predfn$pred_mean,
  lower0 = H0_predfn$pred_mean - err0, 
  lower1 = H1_predfn$pred_mean - err1, 
  upper0 = H0_predfn$pred_mean + err0,
  upper1 = H1_predfn$pred_mean + err1
)
yrange = range(ggdata$lower0, ggdata$lower1, 
               ggdata$upper0, ggdata$upper1)
yrange[1] = yrange[1] - 1
ggdata$Wasserstein = ggdata$Wasserstein * 0.25 - abs(yrange[1])
ggdata$zero1 = NA
ggdata$zero2 = NA
ggdata.melted = melt(ggdata, id.vars = c("x"), 
                     measure.vars = 2:5)
ggdata.lower = melt(ggdata, id.vars = c("x"), 
                    measure.vars = c("zero1", "zero2", "lower0", "lower1"))
ggdata.upper = melt(ggdata, id.vars = c("x"), 
                    measure.vars = c("zero1", "zero2", "upper0", "upper1"))
ggdata.melted = cbind(ggdata.melted, 
                      lower = ggdata.lower$value, 
                      upper = ggdata.upper$value)
ggdata_pts = data.table(
  x = c(x_train, newpts), 
  y = c(y_train, truey), 
  color = c(rep(gg_color_hue(2)[2], length(x_train)), 
            rep(gg_color_hue(2)[1], length(newpts))), 
  shape = c(rep(8, length(x_train)), 
            rep(16, length(newpts)))
)
mvp = ggplot(data = ggdata.melted, aes(x = x, y =value, color = variable), 
       linetype = 1) + 
  geom_path(size = 2) + 
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = variable), 
              alpha = 0.1, linetype = 0) +
  scale_linetype_manual(values = c(1, 1, 2, 2)) + 
  scale_fill_manual(values = c(NA, NA, "#00BFC4", "#C77CFF")) + 
  scale_color_manual(values = c(1, "gray", "#00BFC4", "#C77CFF")) + 
  geom_point(data = ggdata_pts, mapping = aes(x = x, y = y), 
             inherit.aes = FALSE, color = ggdata_pts$color, 
             shape = ggdata_pts$shape, 
             size = 3) +
  geom_point(data = ggdata_pts, mapping = aes(x = x, y = yrange[1]), 
             inherit.aes = FALSE, color = ggdata_pts$color, 
             shape = ggdata_pts$shape, 
             size = 3) +
  scale_y_continuous(limits = yrange) + 
  theme_bw(base_size = 20) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  labs(y = "y", x = "x", fill = "Function", color = "Function")
mvp
# ggsave("mvp.pdf",
#        plot = last_plot(),
#        device = "pdf",
#        path = image_path,
#        scale = 1,
#        width = 8,
#        height = 4,
#        units = c("in")
# )
ggplot(data = ggdata.melted, aes(x = x, y =value, color = variable), 
       linetype = 1) + 
  geom_path() + 
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = variable), 
              alpha = 0.1, linetype = 0) +
  scale_linetype_manual(values = c(1, 1, 2, 2)) + 
  scale_fill_manual(values = c(NA, NA, "#00BFC4", "#C77CFF")) + 
  scale_color_manual(values = c(1, "gray", "#00BFC4", "#C77CFF")) + 
  geom_point(data = ggdata_pts, mapping = aes(x = x, y = y), 
             inherit.aes = FALSE, color = ggdata_pts$color, 
             shape = ggdata_pts$shape, 
             size = 3) +
  geom_point(data = ggdata_pts, mapping = aes(x = x, y = yrange[1]), 
             inherit.aes = FALSE, color = ggdata_pts$color, 
             shape = ggdata_pts$shape, 
             size = 3) +
  scale_y_continuous(limits = yrange) + 
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  labs(y = "y", x = "x", fill = "Function", color = "Function")
# ggsave("mvp.pdf",
#        plot = last_plot(),
#        device = "pdf",
#        path = image_path,
#        scale = 1,
#        width = 8,
#        height = 4,
#        units = c("in")
# )




### --- SMMED --- ###

# different input points
par(mfrow = c(1, 1))
set.seed(1)
inputs = as.data.frame(rbind(x_train1, x_train2, x_train3, runif(N)), 
                       row.names = c("extrapolation", "inc.spread", 
                                     "even.coverage", "random"))
matplot(inputs, pch = 1, col = 1, axes = FALSE)
axis(2)
axis(1, labels = rownames(inputs), at = 1:4)
inputset.names = c("extrap", "incr", "even", "random")

# N
numSteps = 3
stepN = 5


# alternative plot

par(mar=c(0,7,0,0))
# make a list
input_cases <- list("Random" = runif(N), "EvenCoverage" = x_train3, "IncSpread" = x_train2, "Extrapolation" = x_train1)

stripchart(input_cases,
           main="",
           xlab="locations",
           ylab="",
           pch=16, cex = 2, las = 1,
)

#
#
#

ggdata = data.table(
  Extrapolation = x_train1, 
  `Inc Spread` = x_train2, 
  `Even Coverage` = x_train3, 
  `Random` = runif(N)
)
ggdata = melt(ggdata, measure.vars = 1:4)
ggdata$variable = factor(ggdata$variable)
plt0 = ggplot(data = ggdata, aes(x = value, y = variable)) + 
  geom_point(size = 3) +
  theme_bw(base_size = 20) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  labs(y = "", x = "x")
plt0
# ggsave("initial.pdf",
#        plot = last_plot(),
#        device = "pdf",
#        path = image_path,
#        scale = 1,
#        width = 8,
#        height = 4,
#        units = c("in")
# )
ggplot(data = ggdata, aes(x = value, y = variable)) + 
  geom_point(size = 3) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  labs(y = "", x = "x")
# ggsave("initial.pdf",
#        plot = last_plot(),
#        device = "pdf",
#        path = image_path,
#        scale = 1,
#        width = 8,
#        height = 4,
#        units = c("in")
# )

# --- Gaussian vs Matern (gvm) Metrics --- #

# MMED parameters for testing
l01= c(0.01, 0.01)
type01 = c(1, 4)
numCandidates = 1001
i = 4
p = 1
nugget = NULL
alpha = 1

gvm1 = readRDS(paste(home, "/run_designs/gp/gpsims_seq/run_designs_v2/gvm_seq_train1sims.rds", sep = ""))
gvm2 = readRDS(paste(home, "/run_designs/gp/gpsims_seq/run_designs_v2/gvm_seq_train2sims.rds", sep = ""))
gvm3 = readRDS(paste(home, "/run_designs/gp/gpsims_seq/run_designs_v2/gvm_seq_train3sims.rds", sep = ""))
gvm4 = readRDS(paste(home, "/run_designs/gp/gpsims_seq/run_designs_v2/gvm_seq_train4sims.rds", sep = ""))

simcases = list(gvm1, gvm2, gvm3, gvm4)


par(mfrow = c(1, 1))

# RSS Ratio (0/1)

case1_medianLogRSS01_vec = rep(NA, 3)# case 1: extrapolation (points close together)
case2_medianLogRSS01_vec = rep(NA, 3)# case 2: ? (points increasingly spread out more)
case3_medianLogRSS01_vec = rep(NA, 3)# case 3: space-filling
case4_medianLogRSS01_vec = rep(NA, 3)# case 4: uniformly-distributed inputs
for(i in 1:3){
  case1_medianLogRSS01_vec[i] = getMedianLogRSS01(simcases[[1]], i)
  case2_medianLogRSS01_vec[i] = getMedianLogRSS01(simcases[[2]], i)
  case3_medianLogRSS01_vec[i] = getMedianLogRSS01(simcases[[3]], i)
  case4_medianLogRSS01_vec[i] = getMedianLogRSS01(simcases[[4]], i)
}
medianLogRSS01_df = data.frame(rbind(case1_medianLogRSS01_vec, 
                                     case2_medianLogRSS01_vec, 
                                     case3_medianLogRSS01_vec, 
                                     case4_medianLogRSS01_vec), 
                               row.names = row.names(inputs))
colnames(medianLogRSS01_df) = c("smmed", "spacefilling", "random")

matplot(medianLogRSS01_df, type = "b", axes = FALSE, ylab = "Median RSS01", 
        pch = 16:18, lwd = rep(3, 3))
axis(2)
axis(1, labels = rownames(medianLogRSS01_df), at = 1:4)
legend("right", legend = colnames(medianLogRSS01_df), col = 1:3, lty = 1:3,
       pch = 16:18, lwd = rep(3, 3))

#
#
#

ggdata = data.table(
  Extrapolation = case1_medianLogRSS01_vec,
  `Inc Spread` = case2_medianLogRSS01_vec,
  `Even Coverage` = case3_medianLogRSS01_vec,
  `Random` = case4_medianLogRSS01_vec,
  Design = c("SeqMED", "SpaceFilling", "Random")
)
ggdata = melt(ggdata, id.vars = c("Design"))
plt_gvm1 = ggplot(ggdata, aes(x = variable, y = value, group = Design, 
                   color = Design, linetype = Design)) + 
  geom_point(size = 3) + 
  geom_path(size = 2) +
  theme_bw(base_size = 20) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  labs(y = "Median Log(RSS0/RSS1)", x = "Initial Data")
plt_gvm1
# ggsave("gvm_medianlogrss01.pdf",
#        plot = last_plot(),
#        device = "pdf",
#        path = image_path,
#        scale = 1,
#        width = 8,
#        height = 4,
#        units = c("in")
# )
ggplot(ggdata, aes(x = variable, y = value, group = Design, 
                   color = Design, linetype = Design)) + 
  geom_point(size = 3) + 
  geom_path() +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  labs(y = "Median Log(RSS0/RSS1)", x = "Initial Data")
# ggsave("gvm_medianlogrss01.pdf",
#        plot = last_plot(),
#        device = "pdf",
#        path = image_path,
#        scale = 1,
#        width = 8,
#        height = 4,
#        units = c("in")
# )



# Log Predictive Density Ratio (0/1) --- not included

# case1_medianLogPredEvid01_vec = rep(NA, 3)# case 1: extrapolation (points close together)
# case2_medianLogPredEvid01_vec = rep(NA, 3)# case 2: ? (points increasingly spread out more)
# case3_medianLogPredEvid01_vec = rep(NA, 3)# case 3: space-filling
# case4_medianLogPredEvid01_vec = rep(NA, 3)# case 4: uniformly-distributed inputs
# for(i in 1:3){
#   case1_medianLogPredEvid01_vec[i] = getMedianLogPredEvid01(simcases[[1]], i)
#   case2_medianLogPredEvid01_vec[i] = getMedianLogPredEvid01(simcases[[2]], i)
#   case3_medianLogPredEvid01_vec[i] = getMedianLogPredEvid01(simcases[[3]], i)
#   case4_medianLogPredEvid01_vec[i] = getMedianLogPredEvid01(simcases[[4]], i)
# }
# medianLogPredEvid01_df = data.frame(rbind(case1_medianLogPredEvid01_vec, case2_medianLogPredEvid01_vec, case3_medianLogPredEvid01_vec, case4_medianLogPredEvid01_vec), row.names = c("extrapolation", "incr.spread", "even.coverage", "random"))
# colnames(medianLogPredEvid01_df) = c("smmed", "spacefilling", "random")
# 
# matplot(medianLogPredEvid01_df, type = "b", axes = FALSE, 
#         ylab = "Median Log Predictive Density Ratio 01", 
#         pch = 16:18, lwd = rep(3, 3))
# axis(2)
# axis(1, labels = rownames(medianLogPredEvid01_df), at = 1:4)
# legend("right", legend = colnames(medianLogPredEvid01_df), 
#        col = 1:3, lty = 1:3, pch = 16:18, lwd = rep(3, 3))


# posterior probability of H1

case1_medianPPH1_vec = rep(NA, 3)# case 1: extrapolation (points close together)
case2_medianPPH1_vec = rep(NA, 3)# case 2: ? (points increasingly spread out more)
case3_medianPPH1_vec = rep(NA, 3)# case 3: space-filling
case4_medianPPH1_vec = rep(NA, 3)# case 4: uniformly-distributed inputs
for(i in 1:3){
  case1_medianPPH1_vec[i] = getMedianPostProbH1(simcases[[1]], i)
  case2_medianPPH1_vec[i] = getMedianPostProbH1(simcases[[2]], i)
  case3_medianPPH1_vec[i] = getMedianPostProbH1(simcases[[3]], i)
  case4_medianPPH1_vec[i] = getMedianPostProbH1(simcases[[4]], i)
}
medianPPH1_df = data.frame(rbind(case1_medianPPH1_vec, case2_medianPPH1_vec, 
                                 case3_medianPPH1_vec, case4_medianPPH1_vec), 
                           row.names = row.names(inputs))
colnames(medianPPH1_df) = c("smmed", "spacefilling", "random")

matplot(medianPPH1_df, type = "b", ylim = c(0, 1), axes = FALSE, 
        ylab = "Median E[P(H1|X,Y)|X]", pch = 16:18, lwd = rep(3, 3))
axis(2)
axis(1, labels = rownames(medianPPH1_df), at = 1:4)
legend("bottomright", legend = colnames(medianPPH1_df), col = 1:3, lty = 1:3, 
       pch = 16:18, lwd = rep(3, 3))

#
#
#

ggdata = data.table(
  Extrapolation = case1_medianPPH1_vec,
  `Inc Spread` = case2_medianPPH1_vec,
  `Even Coverage` = case3_medianPPH1_vec,
  `Random` = case4_medianPPH1_vec,
  Design = c("SeqMED", "SpaceFilling", "Random")
)
ggdata = melt(ggdata, id.vars = c("Design"))
plt_gvm2 = ggplot(ggdata, aes(x = variable, y = value, group = Design, 
                   color = Design, linetype = Design)) + 
  geom_point(size = 3) + 
  geom_path(size = 2) +
  theme_bw(base_size = 20) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  labs(y = "Median P(H1|X, Y)", x = "Initial Data")
plt_gvm2
# ggsave("gvm_medianpph1.pdf",
#        plot = last_plot(),
#        device = "pdf",
#        path = image_path,
#        scale = 1,
#        width = 8,
#        height = 4,
#        units = c("in")
# )


ggplot(ggdata, aes(x = variable, y = value, group = Design, 
                   color = Design, linetype = Design)) + 
  geom_point(size = 3) + 
  geom_path() +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  labs(y = "Median P(H1|X, Y)", x = "Initial Data")
ggsave("gvm_medianpph1.pdf",
       plot = last_plot(),
       device = "pdf",
       path = image_path,
       scale = 1,
       width = 8,
       height = 4,
       units = c("in")
)


#

#

#



ggdata1 = data.table(
  Extrapolation = case1_medianLogRSS01_vec,
  `Inc Spread` = case2_medianLogRSS01_vec,
  `Even Coverage` = case3_medianLogRSS01_vec,
  `Random` = case4_medianLogRSS01_vec,
  Design = c("SeqMED", "SpaceFilling", "Random")
)
ggdata1 = melt(ggdata1, id.vars = c("Design"))
ggdata1$Metric = "Median Log(RSS0/RSS1)"
ggdata2 = data.table(
  Extrapolation = case1_medianPPH1_vec,
  `Inc Spread` = case2_medianPPH1_vec,
  `Even Coverage` = case3_medianPPH1_vec,
  `Random` = case4_medianPPH1_vec,
  Design = c("SeqMED", "SpaceFilling", "Random")
)
ggdata2 = melt(ggdata2, id.vars = c("Design"))
ggdata2$Metric = "Median P(H1|X,Y)"
ggdata3 = rbind(ggdata1, ggdata2)
ggdata3$Metric = factor(ggdata3$Metric)
plt_gvm3 = ggplot(ggdata3, aes(x = variable, y = value, group = Design, 
                              color = Design, linetype = Design)) + 
  facet_wrap(vars(Metric), scales = "free") + 
  geom_point(size = 3) + 
  geom_path(size = 2) +
  theme_bw(base_size = 20) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.title = element_blank(), 
        axis.text.x = element_text(angle = 45, vjust = 0.5))
plt_gvm3
# ggsave("gvm_medianlogrss01pph1.pdf",
#        plot = last_plot(),
#        device = "pdf",
#        path = image_path,
#        scale = 1,
#        width = 8,
#        height = 4,
#        units = c("in")
# )
# ggsave("poster_gvm_medianlogrss01pph1_h5.png",
#        plot = last_plot(),
#        device = "png",
#        path = image_path,
#        scale = 1,
#        width = 13.5,
#        height = 5,
#        units = c("in")
# )

ggarrange(plt0, plt_gvm3, widths = c(1, 2))
ggsave("poster_gvm_all_h3.png",
       plot = last_plot(),
       device = "png",
       path = image_path,
       scale = 1,
       width = 13.5,
       height = 3,
       units = c("in")
)


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






# --- Matern vs Periodic (gvm) Metrics --- #

# MMED parameters for testing
l01= c(0.01, 0.1); type01 = c(4, 5); numCandidates = 1001; k = 4; p = 1; nugget = NULL; alpha = 1
# gvm sim cases
mvp1 = readRDS(paste(home, "/run_designs/gp/gpsims_seq/run_designs_v2/mvp_seq_train1sims.rds", sep = ""))
mvp2 = readRDS(paste(home, "/run_designs/gp/gpsims_seq/run_designs_v2/mvp_seq_train2sims.rds", sep = ""))
mvp3 = readRDS(paste(home, "/run_designs/gp/gpsims_seq/run_designs_v2/mvp_seq_train3sims.rds", sep = ""))
mvp4 = readRDS(paste(home, "/run_designs/gp/gpsims_seq/run_designs_v2/mvp_seq_train4sims.rds", sep = ""))

simcases = list(mvp1, mvp2, mvp3, mvp4)


# RSS Ratio (0/1)

case1_medianLogRSS01_vec = rep(NA, 3)# case 1: extrapolation (points close together)
case2_medianLogRSS01_vec = rep(NA, 3)# case 2: ? (points increasingly spread out more)
case3_medianLogRSS01_vec = rep(NA, 3)# case 3: space-filling
case4_medianLogRSS01_vec = rep(NA, 3)# case 4: uniformly-distributed inputs
for(i in 1:3){
  case1_medianLogRSS01_vec[i] = getMedianLogRSS01(simcases[[1]], i)
  case2_medianLogRSS01_vec[i] = getMedianLogRSS01(simcases[[2]], i)
  case3_medianLogRSS01_vec[i] = getMedianLogRSS01(simcases[[3]], i)
  case4_medianLogRSS01_vec[i] = getMedianLogRSS01(simcases[[4]], i)
}
medianLogRSS01_df = data.frame(rbind(case1_medianLogRSS01_vec, case2_medianLogRSS01_vec, case3_medianLogRSS01_vec, case4_medianLogRSS01_vec), row.names = c("extrapolation", "incr.spread", "even.coverage", "random"))
colnames(medianLogRSS01_df) = c("smmed", "spacefilling", "random")

matplot(medianLogRSS01_df, type = "b", axes = FALSE, ylab = "Median RSS 01", 
        pch = 16:18, lwd = rep(3, 3))
axis(2)
axis(1, labels = rownames(medianLogRSS01_df), at = 1:4)
legend("topleft", legend = colnames(medianLogRSS01_df), col = 1:3, lty = 1:3, 
       pch = 16:18, lwd = rep(3, 3))

#
#
#

ggdata = data.table(
  Extrapolation = case1_medianLogRSS01_vec,
  `Inc Spread` = case2_medianLogRSS01_vec,
  `Even Coverage` = case3_medianLogRSS01_vec,
  `Random` = case4_medianLogRSS01_vec,
  Design = c("SeqMED", "SpaceFilling", "Random")
)
ggdata = melt(ggdata, id.vars = c("Design"))
plt_mvp1 = ggplot(ggdata, aes(x = variable, y = value, group = Design, 
                              color = Design, linetype = Design)) + 
  geom_point(size = 3) + 
  geom_path(size = 2) +
  theme_bw(base_size = 20) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  labs(y = "Median Log(RSS0/RSS1)", x = "Initial Data")
plt_mvp1
# ggsave("mvp_medianlogrss01.pdf",
#        plot = last_plot(),
#        device = "pdf",
#        path = image_path,
#        scale = 1,
#        width = 8,
#        height = 4,
#        units = c("in")
# )

# Log Predictive Density Ratio (0/1) --- not included

# case1_medianLogPredEvid01_vec = rep(NA, 3)# case 1: extrapolation (points close together)
# case2_medianLogPredEvid01_vec = rep(NA, 3)# case 2: ? (points increasingly spread out more)
# case3_medianLogPredEvid01_vec = rep(NA, 3)# case 3: space-filling
# case4_medianLogPredEvid01_vec = rep(NA, 3)# case 4: uniformly-distributed inputs
# for(i in 1:3){
#   case1_medianLogPredEvid01_vec[i] = getMedianLogPredEvid01(simcases[[1]], i)
#   case2_medianLogPredEvid01_vec[i] = getMedianLogPredEvid01(simcases[[2]], i)
#   case3_medianLogPredEvid01_vec[i] = getMedianLogPredEvid01(simcases[[3]], i)
#   case4_medianLogPredEvid01_vec[i] = getMedianLogPredEvid01(simcases[[4]], i)
# }
# medianLogPredEvid01_df = data.frame(rbind(case1_medianLogPredEvid01_vec, case2_medianLogPredEvid01_vec, case3_medianLogPredEvid01_vec, case4_medianLogPredEvid01_vec), row.names = c("extrapolation", "incr.spread", "even.coverage", "random"))
# colnames(medianLogPredEvid01_df) = c("smmed", "spacefilling", "random")
# 
# matplot(medianLogPredEvid01_df, type = "b", axes = FALSE, ylab = "Median Log Predictive Density Ratio 01", pch = 16:18, lwd = rep(3, 3), cex = rep(2, 3))
# axis(2)
# axis(1, labels = rownames(medianLogPredEvid01_df), at = 1:4)
# legend("topright", legend = colnames(medianLogPredEvid01_df), col = 1:3, lty = 1:3, pch = 16:18, lwd = rep(3, 3), cex = rep(2, 3))

# Posterior Probability of H1

case1_medianPPH1_vec = rep(NA, 3)# case 1: extrapolation (points close together)
case2_medianPPH1_vec = rep(NA, 3)# case 2: ? (points increasingly spread out more)
case3_medianPPH1_vec = rep(NA, 3)# case 3: space-filling
case4_medianPPH1_vec = rep(NA, 3)# case 4: uniformly-distributed inputs
for(i in 1:3){
  case1_medianPPH1_vec[i] = getMedianPostProbH1(simcases[[1]], i)
  case2_medianPPH1_vec[i] = getMedianPostProbH1(simcases[[2]], i)
  case3_medianPPH1_vec[i] = getMedianPostProbH1(simcases[[3]], i)
  case4_medianPPH1_vec[i] = getMedianPostProbH1(simcases[[4]], i)
}
medianPPH1_df = data.frame(rbind(case1_medianPPH1_vec, case2_medianPPH1_vec, case3_medianPPH1_vec, case4_medianPPH1_vec), row.names = c("extrapolation", "incr.spread", "even.coverage", "random"))
colnames(medianPPH1_df) = c("smmed", "spacefilling", "random")

matplot(medianPPH1_df, type = "b", ylim = c(0, 1), axes = FALSE, 
        ylab = "Median E[P(H1|X,Y)|X]", pch = 16:18, lwd = rep(3, 3))
axis(2)
axis(1, labels = rownames(medianPPH1_df), at = 1:4)
legend("bottomright", legend = colnames(medianPPH1_df), col = 1:3, lty = 1:3, 
       pch = 16:18, lwd = rep(3, 3))

#
#
#

ggdata = data.table(
  Extrapolation = case1_medianPPH1_vec,
  `Inc Spread` = case2_medianPPH1_vec,
  `Even Coverage` = case3_medianPPH1_vec,
  `Random` = case4_medianPPH1_vec,
  Design = c("SeqMED", "SpaceFilling", "Random")
)
ggdata = melt(ggdata, id.vars = c("Design"))
plt_mvp2 = ggplot(ggdata, aes(x = variable, y = value, group = Design, 
                              color = Design, linetype = Design)) + 
  geom_point(size = 3) + 
  geom_path(size = 2) +
  theme_bw(base_size = 20) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  labs(y = "Median P(H1|X, Y)", x = "Initial Data")
plt_mvp2
# ggsave("mvp_medianpph1.pdf",
#        plot = last_plot(),
#        device = "pdf",
#        path = image_path,
#        scale = 1,
#        width = 8,
#        height = 4,
#        units = c("in")
# )




#

#

#



ggdata1 = data.table(
  Extrapolation = case1_medianLogRSS01_vec,
  `Inc Spread` = case2_medianLogRSS01_vec,
  `Even Coverage` = case3_medianLogRSS01_vec,
  `Random` = case4_medianLogRSS01_vec,
  Design = c("SeqMED", "SpaceFilling", "Random")
)
ggdata1 = melt(ggdata1, id.vars = c("Design"))
ggdata1$Metric = "Log(RSS0/RSS1)"
ggdata2 = data.table(
  Extrapolation = case1_medianPPH1_vec,
  `Inc Spread` = case2_medianPPH1_vec,
  `Even Coverage` = case3_medianPPH1_vec,
  `Random` = case4_medianPPH1_vec,
  Design = c("SeqMED", "SpaceFilling", "Random")
)
ggdata2 = melt(ggdata2, id.vars = c("Design"))
ggdata2$Metric = "P(H1|X,Y)"
ggdata3 = rbind(ggdata1, ggdata2)
ggdata3$Metric = factor(ggdata3$Metric)
plt_mvp3 = ggplot(ggdata3, aes(x = variable, y = value, group = Design, 
                               color = Design, linetype = Design)) + 
  facet_wrap(vars(Metric), scales = "free") + 
  geom_point(size = 3) + 
  geom_path(size = 2) +
  theme_bw(base_size = 20) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.title = element_blank(), 
        axis.text.x = element_text(angle = 45, vjust = 0.5))
plt_mvp3
ggsave("mvp_medianlogrss01pph1.pdf",
       plot = last_plot(),
       device = "pdf",
       path = image_path,
       scale = 1,
       width = 8,
       height = 4,
       units = c("in")
)
# ggsave("poster_mvp_medianlogrss01pph1_h5.png",
#        plot = last_plot(),
#        device = "png",
#        path = image_path,
#        scale = 1,
#        width = 13.5,
#        height = 5,
#        units = c("in")
# )

ggarrange(mvp, plt_mvp3, widths = c(1.5, 1.5))
# ggsave("poster_mvp_all_h4.png",
#        plot = last_plot(),
#        device = "png",
#        path = image_path,
#        scale = 1,
#        width = 13.5,
#        height = 4,
#        units = c("in")
# )

