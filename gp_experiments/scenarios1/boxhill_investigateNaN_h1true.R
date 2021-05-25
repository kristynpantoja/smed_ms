################################################################################
# last updated: 05/25/2021
# purpose: to test seqmedgp for scenario 1:
#   squared exponential vs. matern,
#   where the true function is matern
# testing box and hill to see when/why/how it gets NaNs

scenario = 2.1 # scenarios: 1.1, 2.1
input.type = 1 # 1 = extrapolation, 2 = inc spread, 3 = even coverage

################################################################################
# Sources/Libraries
################################################################################
scenario_subtypes = unlist(strsplit(as.character(scenario), split = "\\."))
output_home = paste0( # works for scenarios1 OR scenarios2 simulations
  "gp_experiments/scenarios", 
  scenario_subtypes[2],
  "/scenario", scenario, "/outputs")
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

# set up parallelization
library(foreach)
library(future)
library(doFuture)
library(parallel)
registerDoFuture()
nworkers = detectCores()
plan(multisession, workers = nworkers)

library(rngtools)
library(doRNG)
rng.seed = 123 # 123, 345
registerDoRNG(rng.seed)

library(ggplot2)
library(reshape2)
library(dplyr)
library(data.table)
gg_color_hue = function(n) {
  hues = seq(15, 275, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

################################################################################
# simulation settings, shared for both scenarios
################################################################################

# simulations settings
numSims = 25
Nin = 6
numSeq = 15
seqN = 1
Nnew = numSeq * seqN
Nttl = Nin + Nnew
xmin = 0
xmax = 1
numx = 10^3 + 1
x_seq = seq(from = xmin, to = xmax, length.out = numx)
sigmasq_measuremt = 1e-10

# boxhill settings
sigmasq = 1
if(scenario_subtypes[1] == 1){
  nuggets = c(1e-15, sigmasq_measuremt)
} else if(scenario_subtypes[1] == 2){
  nuggets = c(1e-5, sigmasq_measuremt)
}
prior_probs = rep(1 / 2, 2)

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
# Scenario settings
################################################################################
if(scenario_subtypes[1] == 1){
  type01 = c("squaredexponential", "matern")
} else if(scenario_subtypes[1] == 2){
  type01 = c("matern", "periodic")
}
typeT = type01[2]
l01= c(0.01, 0.01)
lT = l01[2]

################################################################################
# models
model0 = list(type = type01[1], l = l01[1], signal.var = sigmasq, 
              measurement.var = nuggets[1])
model1 = list(type = type01[2], l = l01[2], signal.var = sigmasq, 
              measurement.var = nuggets[2])

################################################################################
# import matern functions
filename_append = ""
if(!is.null(sigmasq_measuremt)){
  filename_append = paste0(
    "_noise", strsplit(as.character(sigmasq_measuremt), "-")[[1]][2])
}
simulated.functions = readRDS(paste0(
  output_home,
  "/scenario", scenario, "_simulated_functions", filename_append,
  "_seed", rng.seed,
  ".rds"))
numSims = simulated.functions$numSims
x_seq = simulated.functions$x
numx = length(x_seq)
null_cov = simulated.functions$null_cov
null_mean = simulated.functions$null_mean
y_seq_mat = simulated.functions$function_values_mat

################################################################################
# generate boxhills
# j = 1:3
j = 1

# j : input setting
input.type = j
# input set
if(input.type == 1){
  x_input = x_in1
  x_input_idx = x_in1_idx
} else if(input.type == 2){
  x_input = x_in2
  x_input_idx = x_in2_idx
} else if(input.type == 3){
  x_input = x_in3
  x_input_idx = x_in3_idx
}

# simulations!
i = 1

y_seq = y_seq_mat[ , i]
y_input = y_seq[x_input_idx]
BHres = BHgp_m2(
  y_input, x_input, x_input_idx, prior_probs, model0, model1, Nnew, 
  x_seq, y_seq, noise = TRUE, measurement.var = sigmasq_measuremt)
BHres$x.new

x_new_idx = BHres$x.new.idx
x_new = BHres$x.new
y_new = y_seq[x_new_idx]

# plot
x_input.gg = x_input
y_input.gg = y_input
x_new.gg = x_new
y_new.gg = y_new
H0_predfn = getGPPredictive(x_seq, x_input.gg, y_input.gg, type01[1], l01[1],
                            signal.var = sigmasq, measurement.var = nuggets[1])
H1_predfn = getGPPredictive(x_seq, x_input.gg, y_input.gg, type01[2], l01[2],
                            signal.var = sigmasq, measurement.var = nuggets[2])
err0 = 2 * sqrt(diag(H0_predfn$pred_var))
err1 = 2 * sqrt(diag(H1_predfn$pred_var))
ggdata = data.table(
  x = x_seq, 
  `True Function` = y_seq, 
  `H0 Predictive Mean` = H0_predfn$pred_mean, 
  `H1 Predictive Mean` = H1_predfn$pred_mean,
  lower0 = H0_predfn$pred_mean - 2 * err0, 
  lower1 = H1_predfn$pred_mean - 2 * err1, 
  upper0 = H0_predfn$pred_mean + 2 * err0,
  upper1 = H1_predfn$pred_mean + 2 * err1
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
# 
# # post probs #
# post.probs.gg = data.frame(
#   x = 1:dim(BHres$post.probs)[1],
#   H0 = BHres$post.probs[, 1], 
#   H1 = BHres$post.probs[, 2]
# )
# post.probs.ggm = reshape2::melt(
#   post.probs.gg, id.vars = c("x"))
# ggplot(post.probs.ggm, aes(x = x, y = value)) + 
#   facet_wrap(vars(variable)) + 
#   geom_path() + 
#   theme_bw() + 
#   ylab("posterior probability of hypothesis")
# 
# # plot criterion for first point #
# BHcrit1 = sapply(x_seq, FUN = function(x) 
#   BHDgp_m2(y_input, x_input, prior_probs, x, model0, model1))
# BHD.gg = data.frame(
#   x = x_seq, 
#   y = BHcrit1
# )
# ggplot(data = ggdata.melted, aes(x = x, y =value, color = variable), 
#        linetype = 1) + 
#   geom_path() + 
#   geom_ribbon(aes(ymin = lower, ymax = upper, fill = variable), 
#               alpha = 0.1, linetype = 0) +
#   scale_linetype_manual(values = c(1, 1, 2, 2)) + 
#   scale_fill_manual(values = c(NA, "#00BFC4", "#C77CFF")) + 
#   scale_color_manual(values = c(1, "#00BFC4", "#C77CFF")) + 
#   geom_point(data = ggdata_pts, mapping = aes(x = x, y = y), 
#              inherit.aes = FALSE, color = ggdata_pts$color, 
#              shape = ggdata_pts$shape, 
#              size = 3, alpha = 0.25) +
#   geom_label(label = c(rep(NA, length(x_input.gg)), c(1:length(x_new.gg))), 
#              data = ggdata_pts, mapping = aes(x = x, y = y + 0.3), 
#              inherit.aes = FALSE, alpha = 0.25) +
#   geom_point(data = ggdata_pts, mapping = aes(x = x, y = yrange[1]), 
#              inherit.aes = FALSE, color = ggdata_pts$color, 
#              shape = ggdata_pts$shape, 
#              size = 3, alpha = 0.25) + 
#   geom_path(data = BHD.gg, aes(x = x, y = (1/1e7) * y - 3), inherit.aes = FALSE, color = 3) + 
#   scale_y_continuous(limits = yrange) +
#   theme_bw() + 
#   theme(panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank()) +
#   labs(y = "y", x = "x", fill = "Function", color = "Function")
# 
# BHDs = matrix(NA, nrow = length(x_seq), ncol = Nnew)
# KLijs = matrix(NA, nrow = length(x_seq), ncol = Nnew)
# KLjis = matrix(NA, nrow = length(x_seq), ncol = Nnew)
# for(i in 1:Nnew){
#   xs = c(x_input, BHres$x.new)
#   ys = c(y_input, BHres$y.new)
#   x.tmp = xs[1:length(x_input) + i - 1]
#   y.tmp = ys[1:length(x_input) + i - 1]
#   probs.tmp = BHres$post.probs[i, ]
#   BHDs[, i] = sapply(x_seq, function(x) BHDgp_m2_testing(
#     y.tmp, x.tmp, probs.tmp, x, model0, model1)$bhd)
#   KLijs[, i] = sapply(x_seq, function(x) BHDgp_m2_testing(
#     y.tmp, x.tmp, probs.tmp, x, model0, model1)$KLij)
#   KLjis[, i] = sapply(x_seq, function(x) BHDgp_m2_testing(
#     y.tmp, x.tmp, probs.tmp, x, model0, model1)$KLji)
# }
# 
# colnames(BHDs) = paste0("X", 1:Nnew)
# colnames(KLijs) = colnames(BHDs)
# colnames(KLjis) = colnames(BHDs)
# BHD.df = data.frame(rbind(BHDs, KLijs, KLjis))
# BHD.df$candidate = rep(x_seq, 3)
# BHD.df$type = c(rep("BH", numx), rep("KLij", numx), rep("KLji", numx))
# BHD.mlt = reshape2::melt(
#   BHD.df, id.vars = c("candidate", "type"), measure.vars = colnames(BHDs))
# 
# # plot BH criterion at particular new points
# ggplot(dplyr::filter(BHD.mlt, type == "BH", variable %in% paste0("X", 1:4)), 
#        aes(x = candidate, y = value, color = type)) + 
#   facet_wrap(vars(variable), scales = "free_y") +
#   geom_vline(xintercept = x_input, color = "gray") + 
#   geom_path() + 
#   theme(text = element_text(size = 15)) + 
#   xlim(0, 0.5)
# 
# # plot BH criterion at first point, and also get the 
# idx = 1
# BHDx.mlt = dplyr::filter(BHD.mlt, type == "BH", variable == paste0("X", idx))
# ggplot(BHDx.mlt, 
#        aes(x = candidate, y = value)) + 
#   geom_vline(xintercept = x_input, color = "gray", size = 1) + 
#   geom_vline(xintercept = dplyr::filter(BHDx.mlt, is.nan(value))$candidate, 
#              color = "orchid", linetype = 2, size = 1) + 
#   geom_path(size = 1) + 
#   theme(text = element_text(size = 15)) + 
#   xlim(0, 0.5)
# 
# # what do the KLs look like in the criterion? pick a specific new point
# ggdata.tmp = dplyr::filter(BHD.mlt, variable == paste0("X", idx))
# ggplot(ggdata.tmp, 
#        aes(x = candidate, y = log(value), color = type)) + 
#   facet_wrap(vars(type)) +
#   geom_vline(xintercept = x_input, color = "gray", size = 1) + 
#   geom_vline(xintercept = dplyr::filter(BHDx.mlt, is.nan(value))$candidate, 
#              color = "orchid", linetype = 2, size = 1) + 
#   geom_path(size = 1) + 
#   theme(text = element_text(size = 15)) + 
#   # theme_classic() + 
#   xlim(0, 0.5)
# # geom_point()
