################################################################################
# last updated: 03/17/2021
# purpose: to create a list of seqmed simulations for scenario 1:
#   squared exponential vs. matern,
#   where the true function is matern

################################################################################
# Sources/Libraries
################################################################################
output_home = "run_designs/updated_simulations/gp/seqmed"
functions_home = "functions"

# for seqmed design
source(paste(functions_home, "/SeqMEDgp.R", sep = ""))
source(paste(functions_home, "/SeqMEDgp_batch.R", sep = ""))
source(paste(functions_home, "/charge_function_q.R", sep = ""))
source(paste(functions_home, "/covariance_functions.R", sep = ""))
source(paste(functions_home, "/wasserstein_distance.R", sep = ""))
source(paste(functions_home, "/gp_predictive.R", sep = ""))

# for box-hill design
source(paste(functions_home, "/boxhill.R", sep = ""))
source(paste(functions_home, "/boxhill_gp.R", sep = ""))
source(paste(functions_home, "/kl_divergence.R", sep = ""))

library(mvtnorm)

# for plots
library(ggplot2)
library(ggpubr)
library(reshape2)
library(data.table)
gg_color_hue = function(n) {
  hues = seq(15, 275, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

# # set up parallelization
# library(doFuture)
# library(parallel)
# registerDoFuture()
# nworkers = detectCores()
# plan(multisession, workers = nworkers)
# library(doRNG)

gg_color_hue = function(n) {
  hues = seq(15, 275, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

################################################################################
# simulation settings, shared for both scenarios
################################################################################

# simulations settings
numSims = 1 # DEMO SETTING ONLY ###############################################
seed = 1 # DEMO SETTING ONLY ###################################################
Nin = 6
numSeq = 2
seqN = 15
Nnew = (numSeq - 1) * seqN
Nttl = Nin + Nnew
xmin = 0
xmax = 1
numx = 10^3 + 1
x_seq = seq(from = xmin, to = xmax, length.out = numx)

# SeqMED settings
nuggetSM = NULL # DEMO SETTING ONLY ############################################

# boxhill settings
prior_probs = rep(1 / 2, 2)
nuggetBH = 1e-10

################################################################################
# input data
################################################################################

# 1. make space-filling design
# space_filling = seq(from = xmin, to = xmax, length.out = Nttl)
space_filling_idx = c(1, 1 + ((numx - 1)/(Nttl - 1)) * 1:((numx - 1) / ((numx - 1)/(Nttl - 1))))
space_filling = x_seq[space_filling_idx]

# input set 1 (extrapolation)
x_in1_idx = space_filling_idx[1:Nin]
x_in1 = x_seq[x_in1_idx]
x_spacefill1_idx = space_filling_idx[-c(1:Nin)]
x_spacefill1 = x_seq[x_spacefill1_idx]
# all.equal(space_filling, c(x_in1, x_spacefill1))

# input set 2 (increasing spread)
x_in2_idx = space_filling_idx[c(1, 2, 4, 7, 12, 21)]
x_in2 = x_seq[x_in2_idx]
x_spacefill2_idx = space_filling_idx[-c(1, 2, 4, 7, 12, 21)]
x_spacefill2 = x_seq[x_spacefill2_idx]
# all.equal(space_filling, sort(c(x_in2, x_spacefill2)))

# input set 3 (space-filling / even coverage)
x_in3_idx = c(1, 1 + ((numx - 1)/(Nin - 1)) * 1:((numx - 1) / ((numx - 1)/(Nin - 1))))
x_in3 = x_seq[x_in3_idx]
x_spacefill3_idx = space_filling_idx[!(space_filling_idx %in% x_in3_idx)]
x_spacefill3 = x_seq[x_spacefill3_idx]
# all.equal(space_filling, sort(c(x_in3, x_spacefill3)))

# input set 4 (uniform / random)

################################################################################
# Scenario 1: Squared exponential vs. matern, true = matern
################################################################################
type01 = c("squaredexponential", "matern")
l01= c(0.1, 0.1) # DEMO SETTING ONLY ###########################################
# generate matern functions
set.seed(seed)
null_cov = getCov(x_seq, x_seq, type01[2], l01[2])
null_mean = rep(0, numx)
y_seq = t(rmvnorm(n = 1, mean = null_mean, sigma = null_cov)) # the function values

# bh settings
model0 = list(type = type01[1], l = l01[1])
model1 = list(type = type01[2], l = l01[2])

# generate seqmed for demo #####################################################

# input set
x_input_idx = sample(1:numx, Nin)
x_input = x_seq[x_input_idx]
y_input = y_seq[x_input_idx]
seqmed.res = SeqMEDgp(
  y0 = y_input, x0 = x_input, x0.idx = x_input_idx, candidates = x_seq,
  function.values = y_seq, nugget = nuggetSM, type = type01, l = l01,
  numSeq = numSeq, seqN = seqN, prints = TRUE
  , seed = 1234 # DEMO SETTING ONLY ############################################
)

# get sim info
sim_ind = 1
x_in = x_input
x_in_idx = x_input_idx
mmed_gp = seqmed.res
y_in = y_input

newpts = mmed_gp$x.new
truey = y_seq[mmed_gp$x.new.idx]

H0_predfn = getGPPredictive(x_seq, x_in, y_in, type01[1], l01[1],
                            nugget = nuggetSM)
H1_predfn = getGPPredictive(x_seq, x_in, y_in, type01[2], l01[2],
                            nugget = nuggetSM)

# get w_seq
Kinv0 = solve(getCov(x_in, x_in, type01[1], l01[1]))
Kinv1 = solve(getCov(x_in, x_in, type01[2], l01[2]))
w_seq = sapply(x_seq, FUN = function(x1) 
  WNgp(x1, Kinv0, Kinv1, x_in, y_in, var_e = 1, type01, l01))
# plot
err0 = 2 * sqrt(diag(H0_predfn$pred_var))
err1 = 2 * sqrt(diag(H1_predfn$pred_var))
ggdata = data.table(
  x = x_seq, 
  `True Function` = as.vector(y_seq), 
  Wasserstein = w_seq, 
  `H0 Predictive` = H0_predfn$pred_mean, 
  `H1 Predictive` = H1_predfn$pred_mean,
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
ggdata.melted = melt(
  ggdata, id.vars = c("x"), 
  measure.vars = c("True Function", "Wasserstein", "H0 Predictive", "H1 Predictive"))
ggdata.lower = melt(ggdata, id.vars = c("x"), 
                    measure.vars = c("zero1", "zero2", "lower0", "lower1"))
ggdata.upper = melt(ggdata, id.vars = c("x"), 
                    measure.vars = c("zero1", "zero2", "upper0", "upper1"))
ggdata.melted = cbind(ggdata.melted, 
                      lower = ggdata.lower$value, 
                      upper = ggdata.upper$value)
ggdata_pts = data.table(
  x = c(x_in, newpts), 
  y = c(y_in, truey), 
  color = c(rep(gg_color_hue(2)[2], length(x_in)), 
            rep(gg_color_hue(2)[1], length(newpts))), 
  shape = c(rep(8, length(x_in)), 
            rep(16, length(newpts)))
)
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
             size = 2) +
  geom_point(data = ggdata_pts, mapping = aes(x = x, y = yrange[1]), 
             inherit.aes = FALSE, color = ggdata_pts$color, 
             shape = ggdata_pts$shape, 
             size = 2) +
  scale_y_continuous(limits = yrange) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  labs(y = "y", x = "x", fill = "Function", color = "Function")
























## plot!
x_new_idx = seqmed.res$x.new.idx
x_new = seqmed.res$x.new
y_new = seqmed.res$y.new

# plot
x_input.gg = x_input
y_input.gg = y_input
x_new.gg = x_new
y_new.gg = y_new
H0_predfn = getGPPredictive(x_seq, x_input.gg, y_input.gg, type01[1], l01[1],
                            nugget = nuggetSM)
H1_predfn = getGPPredictive(x_seq, x_input.gg, y_input.gg, type01[2], l01[2],
                            nugget = nuggetSM)
err0 = 2 * sqrt(diag(H0_predfn$pred_var))
err1 = 2 * sqrt(diag(H1_predfn$pred_var))
ggdata = data.table(
  x = x_seq,
  `True Function` = y_seq,
  `H0 Predictive Mean` = H0_predfn$pred_mean,
  `H1 Predictive Mean` = H1_predfn$pred_mean,
  lower0 = H0_predfn$pred_mean - err0,
  lower1 = H1_predfn$pred_mean - err1,
  upper0 = H0_predfn$pred_mean + err0,
  upper1 = H1_predfn$pred_mean + err1
)
yrange = range(ggdata$lower0, ggdata$lower1,
               ggdata$upper0, ggdata$upper1,
               na.rm = TRUE)
yrange[1] = yrange[1] - 1
ggdata$zero1 = NA
ggdata.melted = melt(ggdata, id.vars = c("x"),
                     measure.vars = c("True Function", "H0 Predictive Mean", "H1 Predictive Mean"))
ggdata.lower = melt(ggdata, id.vars = c("x"),
                    measure.vars = c("zero1", "lower0", "lower1"))
ggdata.upper = melt(ggdata, id.vars = c("x"),
                    measure.vars = c("zero1", "upper0", "upper1"))
ggdata.melted = cbind(ggdata.melted,
                      lower = ggdata.lower$value,
                      upper = ggdata.upper$value)
ggdata_pts = data.table(
  x = c(x_input.gg, x_new.gg),
  y = c(y_input.gg, y_new.gg),
  color = c(rep(gg_color_hue(2)[2], length(x_input.gg)),
            rep(gg_color_hue(2)[1], length(x_new.gg))),
  shape = c(rep(8, length(x_input.gg)),
            rep(16, length(x_new.gg)))
)
ggplot(data = ggdata.melted, aes(x = x, y =value, color = variable),
       linetype = 1) +
  geom_path() +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = variable),
              alpha = 0.1, linetype = 0) +
  scale_linetype_manual(values = c(1, 1, 2, 2)) +
  scale_fill_manual(values = c(NA, "#00BFC4", "#C77CFF")) +
  scale_color_manual(values = c(1, "#00BFC4", "#C77CFF")) +
  geom_point(data = ggdata_pts, mapping = aes(x = x, y = y),
             inherit.aes = FALSE, color = ggdata_pts$color,
             shape = ggdata_pts$shape,
             size = 3, alpha = 0.25) +
  geom_label(label = c(rep(NA, length(x_input.gg)), c(1:length(x_new.gg))),
             data = ggdata_pts, mapping = aes(x = x, y = y + 0.3),
             inherit.aes = FALSE, alpha = 0.25) +
  geom_point(data = ggdata_pts, mapping = aes(x = x, y = yrange[1]),
             inherit.aes = FALSE, color = ggdata_pts$color,
             shape = ggdata_pts$shape,
             size = 3, alpha = 0.25) +
  scale_y_continuous(limits = yrange) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  labs(y = "y", x = "x", fill = "Function", color = "Function")


