
scenario = 4 # scenarios: 3,4,5,6
seq.type = 1 # 1 = fully sequential, 2 = stage-sequential 3x5

################################################################################
# Sources/Libraries
################################################################################
sims_dir = "gp_experiments/simulations_MSP"
output_dir = paste0(sims_dir, "/scenario_MSP/outputs")
data_dir = paste0(sims_dir, "/simulated_data")
functions_dir = "functions"

# for seqmed design
source(paste(functions_dir, "/SeqMEDgp.R", sep = ""))
source(paste(functions_dir, "/SeqMEDgp_batch.R", sep = ""))
source(paste(functions_dir, "/charge_function_q.R", sep = ""))
source(paste(functions_dir, "/covariance_functions.R", sep = ""))
source(paste(functions_dir, "/wasserstein_distance.R", sep = ""))
source(paste(functions_dir, "/gp_predictive.R", sep = ""))
source(paste(functions_dir, "/gp_plot.R", sep = ""))

# for box-hill design
source(paste(functions_dir, "/boxhill.R", sep = ""))
source(paste(functions_dir, "/boxhill_gp.R", sep = ""))
source(paste(functions_dir, "/kl_divergence.R", sep = ""))

library(mvtnorm)
rng.seed = 123

library(ggplot2)
library(reshape2)
library(ggpubr)
gg_color_hue = function(n) {
  hues = seq(15, 275, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
library(data.table)

################################################################################
# simulation settings, shared for both scenarios
################################################################################

# simulations settings
numSims = 25
Nin = 1
if(seq.type == 1){
  numSeq = 15
  seqN = 1
} else if(seq.type == 2){
  numSeq = 3
  seqN = 5
}
Nnew = numSeq * seqN
Nttl = Nin + Nnew
xmin = 0
xmax = 1
numx = 10^3 + 1
x_seq = seq(from = xmin, to = xmax, length.out = numx)
sigmasq_measuremt = 1e-10
sigmasq_signal = 1

# shared settings
nugget = sigmasq_measuremt
prior_probs = rep(1 / 2, 2)

################################################################################
# Scenario settings
################################################################################
if(scenario == 3){
  type01 = c("squaredexponential", "squaredexponential")
  typeT = "matern"
  l01= c(0.005, 0.01)
  lT = 0.01
  model0 = list(type = type01[1], l = l01[1], signal.var = sigmasq_signal, 
                measurement.var = nugget)
  model1 = list(type = type01[2], l = l01[2], signal.var = sigmasq_signal, 
                measurement.var = nugget)
} else if(scenario == 4){
  type01 = c("matern", "squaredexponential")
  typeT = "periodic"
  l01= c(0.01, 0.01)
  lT = 0.1
  pT = 0.05 # 0.05 or 0.1
  model0 = list(type = type01[1], l = l01[1], signal.var = sigmasq_signal, 
                measurement.var = nugget)
  model1 = list(type = type01[2], l = l01[2], signal.var = sigmasq_signal, 
                measurement.var = nugget)
} else if(scenario == 5){
  type01 = c("matern", "periodic")
  typeT = "squaredexponential"
  l01= c(0.01, 0.01)
  lT = 0.01
  p1 = 0.26
  model0 = list(type = type01[1], l = l01[1], signal.var = sigmasq_signal, 
                measurement.var = nugget)
  model1 = list(type = type01[2], l = l01[2], signal.var = sigmasq_signal, 
                measurement.var = nugget, p = p1)
} else if(scenario == 6){
  type01 = c("squaredexponential", "periodic")
  typeT = "matern"
  l01= c(0.01, 0.01)
  lT = 0.01
  p1 = 0.26
  model0 = list(type = type01[1], l = l01[1], signal.var = sigmasq_signal, 
                measurement.var = nugget)
  model1 = list(type = type01[2], l = l01[2], signal.var = sigmasq_signal, 
                measurement.var = nugget, p = p1)
} else{
  stop("invalid scenario number")
}

################################################################################
# import data
filename_append = ""
if(!is.null(sigmasq_measuremt)){
  filename_append = "_noise"
}
if(typeT == "periodic"){
  simulated_data_file = paste0(
    data_dir,
    "/", typeT,
    "_l", lT,
    "_p", pT,
    filename_append, 
    "_seed", rng.seed,
    ".rds")
} else{
  simulated_data_file = paste0(
    data_dir,
    "/", typeT,
    "_l", lT,
    filename_append, 
    "_seed", rng.seed,
    ".rds")
}
simulated.data = readRDS(simulated_data_file)
numSims = simulated.data$numSims
x_seq = simulated.data$x
numx = length(x_seq)
null_cov = simulated.data$null_cov
null_mean = simulated.data$null_mean
y_seq_mat = simulated.data$function_values_mat

################################################################################
# initial design
x_input_idx = ceiling(numx / 2)
x_input = x_seq[x_input_idx]

################################################################################
# read in the data

filename_append = ""
if(!is.null(sigmasq_measuremt)){
  filename_append = "_noise"
}
if(typeT == "periodic"){
  simulated_data_file = paste0(
    data_dir,
    "/", typeT,
    "_l", lT,
    "_p", pT,
    filename_append, 
    "_seed", rng.seed,
    ".rds")
} else{
  simulated_data_file = paste0(
    data_dir,
    "/", typeT,
    "_l", lT,
    filename_append, 
    "_seed", rng.seed,
    ".rds")
}
simulated.data = readRDS(simulated_data_file)
numSims = simulated.data$numSims
x_seq = simulated.data$x
numx = length(x_seq)
null_cov = simulated.data$null_cov
null_mean = simulated.data$null_mean
y_seq_mat = simulated.data$function_values_mat

################################################################################
# initial design
x_input_idx = ceiling(numx / 2)
x_input = x_seq[x_input_idx]

################################################################################
# read in the data

# filename_append.tmp for all methods alike
filename_append.tmp = paste0(
  filename_append, 
  "_seed", rng.seed,
  ".rds"
)
boxhill_sims = readRDS(paste0(
  output_dir,
  "/scenario", scenario, "_boxhill", 
  filename_append.tmp))
leaveout_sims = readRDS(paste0(
  output_dir,
  "/scenario", scenario, "_seqmed", 
  "_leaveout", 
  "_seq", seq.type,
  filename_append.tmp))
qcap_sims = readRDS(paste0(
  output_dir,
  "/scenario", scenario, "_seqmed", 
  "_cap",
  "_seq", seq.type,
  filename_append.tmp))
persist_sims = readRDS(paste0(
  output_dir,
  "/scenario", scenario, "_seqmed", 
  "_persist", 
  "_seq", seq.type,
  filename_append.tmp))
qcap_persist_sims = readRDS(paste0(
  output_dir,
  "/scenario", scenario, "_seqmed", 
  "_cap_persist",
  "_seq", seq.type,
  filename_append.tmp))

if(typeT == "periodic"){
  random_sims_file = paste0(
    sims_dir, 
    "/spacefilling_designs/outputs/random", 
    "_", typeT,
    "_l", lT,
    "_p", pT,
    filename_append.tmp)
  grid_sims_file = paste0(
    sims_dir,
    "/spacefilling_designs/outputs/grid", 
    "_", typeT,
    "_l", lT,
    "_p", pT,
    filename_append.tmp)
} else{
  random_sims_file = paste0(
    sims_dir, 
    "/spacefilling_designs/outputs/random", 
    "_", typeT,
    "_l", lT,
    filename_append.tmp)
  grid_sims_file = paste0(
    sims_dir,
    "/spacefilling_designs/outputs/grid", 
    "_", typeT,
    "_l", lT,
    filename_append.tmp)
}
random_sims = readRDS(random_sims_file)
grid_sims = readRDS(grid_sims_file)

################################################################################
# make posterior probability of H_\ell plots
################################################################################

# models
if(typeT == "periodic"){
  modelT = list(type = typeT, l = lT, signal.var = sigmasq_signal, 
                measurement.var = sigmasq_measuremt, p = pT)
} else{
  modelT = list(type = typeT, l = lT, signal.var = sigmasq_signal, 
                measurement.var = sigmasq_measuremt)
}

getPPHseq = function(design, model0, model1, modelT){
  PPH0_seq = rep(NA, length(as.vector(na.omit(design$y.new))))
  PPH1_seq = rep(NA, length(as.vector(na.omit(design$y.new))))
  PPHT_seq = rep(NA, length(as.vector(na.omit(design$y.new))))
  for(i in 1:length(as.vector(na.omit(design$y.new)))){
    y.tmp = c(design$y, as.vector(na.omit(design$y.new))[1:i])
    x.tmp = c(design$x, as.vector(na.omit(design$x.new))[1:i])
    PPHs.tmp = getHypothesesPosteriors(
      prior.probs = rep(1 / 3, 3), 
      evidences = c(
        Evidence_gp(y.tmp, x.tmp, model0),
        Evidence_gp(y.tmp, x.tmp, model1), 
        Evidence_gp(y.tmp, x.tmp, modelT)
      )
    )
    PPH0_seq[i] = PPHs.tmp[1]
    PPH1_seq[i] = PPHs.tmp[2]
    PPHT_seq[i] = PPHs.tmp[3]
  }
  if(length(PPH0_seq) < Nnew){
    PPH0_seq[(length(PPH0_seq) + 1):Nnew] = PPH0_seq[length(PPH0_seq)]
  }
  if(length(PPH1_seq) < Nnew){
    PPH1_seq[(length(PPH1_seq) + 1):Nnew] = PPH1_seq[length(PPH1_seq)]
  }
  if(length(PPHT_seq) < Nnew){
    PPHT_seq[(length(PPHT_seq) + 1):Nnew] = PPHT_seq[length(PPHT_seq)]
  }
  return(data.frame(
    index = 1:Nnew, 
    "H0" = PPH0_seq, 
    "H1" = PPH1_seq, 
    "HT" = PPHT_seq
  ))
}

PPH_seq = data.frame(
  PPH0 = numeric(), PPH1 = numeric(), PPHT = numeric(), 
  type = character(), sim = numeric())
for(j in 1:numSims){
  # designs at sim b
  bh = boxhill_sims[[j]]
  qc = qcap_sims[[j]]
  lo = leaveout_sims[[j]]
  qc2 = qcap_persist_sims[[j]]
  kq = persist_sims[[j]]
  r = random_sims[[j]]
  g = grid_sims[[j]]
  # sequence of PPHs for each design
  PPH_seq.bh = getPPHseq(bh, model0, model1, modelT)
  PPH_seq.qc = getPPHseq(qc, model0, model1, modelT)
  PPH_seq.lo = getPPHseq(lo, model0, model1, modelT)
  PPH_seq.qc2 = getPPHseq(qc2, model0, model1, modelT)
  PPH_seq.kq = getPPHseq(kq, model0, model1, modelT)
  PPH_seq.r = getPPHseq(r, model0, model1, modelT)
  PPH_seq.g = getPPHseq(g, model0, model1, modelT)
  # master data frame
  PPH_seq.bh$type = "boxhill"
  PPH_seq.qc$type = "qcap"
  PPH_seq.lo$type = "leaveout"
  PPH_seq.qc2$type = "keepq2"
  PPH_seq.kq$type = "keepq"
  PPH_seq.r$type = "random"
  PPH_seq.g$type = "grid"
  PPH_seq.tmp = rbind(
    PPH_seq.bh, PPH_seq.qc, PPH_seq.lo, PPH_seq.qc2, PPH_seq.kq, 
    PPH_seq.r, PPH_seq.g)
  PPH_seq.tmp$sim = j
  PPH_seq = rbind(PPH_seq, PPH_seq.tmp)
}

PPH0mean_seq = aggregate(PPH_seq$H0, by = list(PPH_seq$index, PPH_seq$type), 
                         FUN = function(x) mean(x, na.rm = TRUE))
names(PPH0mean_seq) = c("index", "type", "value")
PPH0mean_seq$Hypothesis = "H0"
PPH1mean_seq = aggregate(PPH_seq$H1, by = list(PPH_seq$index, PPH_seq$type), 
                         FUN = function(x) mean(x, na.rm = TRUE))
names(PPH1mean_seq) = c("index", "type", "value")
PPH1mean_seq$Hypothesis = "H1"
PPHTmean_seq = aggregate(PPH_seq$HT, by = list(PPH_seq$index, PPH_seq$type), 
                         FUN = function(x) mean(x, na.rm = TRUE))
names(PPHTmean_seq) = c("index", "type", "value")
PPHTmean_seq$Hypothesis = "HT"

# plot mean curves
PPHmean_seq = rbind(PPH0mean_seq, PPH1mean_seq, PPHTmean_seq)
epph.plt = ggplot(PPHmean_seq, aes(x = index, y = value, color = type, 
                                   linetype = type, shape = type)) + 
  facet_wrap(~Hypothesis) + 
  geom_path() + 
  geom_point() +
  theme_bw() +
  ylim(0, 1)
plot(epph.plt)

# plot the numSims individual posterior probability of HT curves
PPH_seq_bh = dplyr::filter(PPH_seq, type == "boxhill")
ggplot(PPH_seq_bh, aes(x = index, y = HT)) +
  facet_wrap(~sim) +
  geom_path() + 
  geom_point() +
  theme_bw() +
  ylim(0, 1) + 
  theme(panel.grid.minor = element_blank())

################################################################################
# make design plots
################################################################################

plt_list = list()
for(j in 1:numSims){
  print(boxhill_sims[[j]]$x.new)
  data.gg = data.frame(
    index = as.character(0:Nnew), 
    design = c(x_input, boxhill_sims[[j]]$x.new), 
    type = c("input", rep("boxhill", Nnew))
  )
  text.gg = data.frame(N = length(na.omit(boxhill_sims[[j]]$x.new)))
  plt_list[[j]] = ggplot() + 
    geom_point(data = data.gg, 
               mapping = aes(x = design, y = 0, color = type)) +
    geom_text(data = text.gg,
              aes(x = x_input, label = N, y = 0.01), size = 5) +
    xlim(xmin, xmax) + 
    ylim(0, 0.02) + 
    theme_classic() +
    theme(legend.position = "none", 
          axis.title.x = element_blank(), axis.text.x = element_blank(), 
          axis.title.y = element_blank(), axis.text.y = element_blank())
  # plot(plt_list[[j]])
}
# ggarrange(plotlist = plt_list, nrow = 5, ncol = 5)

################################################################################
# investigate one simulation only 
################################################################################

idx = 1
y_seq = y_seq_mat[, idx]
y_input = y_seq[x_input_idx]
x_new = boxhill_sims[[idx]]$x.new
x_new_idx = boxhill_sims[[idx]]$x.new.idx
# x_new = qcap_sims[[idx]]$x.new
# x_new_idx = qcap_sims[[idx]]$x.new.idx
y_new = y_seq[x_new_idx]

fit_new_idx = c(1)
x_fit = c(x_input, x_new[fit_new_idx])
y_fit = c(y_input, y_new[fit_new_idx])

# ggplot(data.frame(x = x_seq, y = y_seq), aes(x = x, y = y)) + 
#   geom_vline(xintercept = c(0.26 + 0:2 * 0.26), color = gg_color_hue(3)[1]) +
#   geom_path() + 
#   geom_point(data = data.frame(x = x_input, y = y_input), 
#              mapping = aes(x = x, y = y), 
#              color = gg_color_hue(3)[2], size = 3) + 
#   geom_point(data = data.frame(x = x_new, y = y_new), 
#              mapping = aes(x = x, y = y), 
#              color = gg_color_hue(3)[3], size = 3) + 
#   theme_classic() # + 
#   # xlim(min(c(x_input, x_new), na.rm = TRUE), 
#   #      max(c(x_input, x_new), na.rm = TRUE))


# plot with error bars

HT_predfn = getGPPredictive(x_seq, x_fit, y_fit, typeT, lT, 1, 
                            measurement.var = sigmasq_measuremt)
err = 2 * sqrt(diag(HT_predfn$pred_var))
ggdata = data.table(
  x = x_seq, 
  `True Function` = y_seq, 
  `HT Pred Mean` = HT_predfn$pred_mean, 
  lower = HT_predfn$pred_mean - err, 
  upper = HT_predfn$pred_mean + err
)
# ggdata.melted = melt(ggdata, id.vars = c("x"))

H0_predfn = getGPPredictive(x_seq, x_fit, y_fit, type01[1], l01[1], 1, 
                            measurement.var = sigmasq_measuremt)
err0 = 2 * sqrt(diag(H0_predfn$pred_var))
ggdata0 = data.table(
  x = x_seq, 
  `True Function` = y_seq, 
  `HT Pred Mean` = H0_predfn$pred_mean, 
  lower = H0_predfn$pred_mean - err0, 
  upper = H0_predfn$pred_mean + err0
)
# ggdata0.melted = melt(ggdata0, id.vars = c("x"))

H1_predfn = getGPPredictive(x_seq, x_fit, y_fit, type01[2], l01[2], 1, 
                            measurement.var = sigmasq_measuremt)
err1 = 2 * sqrt(diag(H1_predfn$pred_var))
ggdata1 = data.table(
  x = x_seq, 
  `True Function` = y_seq, 
  `HT Pred Mean` = H1_predfn$pred_mean, 
  lower = H1_predfn$pred_mean - err1, 
  upper = H1_predfn$pred_mean + err1
)
# ggdata1.melted = melt(ggdata1, id.vars = c("x"))

ggdata_pts = data.table(
  x = c(x_input, x_new), y = c(y_input, y_new), 
  color = c(rep(gg_color_hue(5)[2], length(x_input) + length(fit_new_idx)), 
            rep(gg_color_hue(5)[1], length(x_new) - length(fit_new_idx))), 
  shape = c(rep(8, length(x_input) + length(fit_new_idx)), 
            rep(16, length(x_new) - length(fit_new_idx)))
)
ggplot(ggdata) + 
  geom_path(aes(x = x, y = `True Function`)) + 
  geom_ribbon(aes(x = x, ymin = lower, ymax = upper), 
              color = gg_color_hue(5)[3], linetype = 2, 
              fill = gg_color_hue(5)[3], alpha = 0.1) +
  geom_path(aes(x = x, y = `HT Pred Mean`), color = gg_color_hue(5)[3]) +
  geom_ribbon(data = ggdata0, 
              aes(x = x, ymin = lower, ymax = upper), 
              color = gg_color_hue(5)[4], linetype = 2, 
              fill = gg_color_hue(5)[4], alpha = 0.1) +
  geom_path(data = ggdata0, 
            aes(x = x, y = `HT Pred Mean`), color = gg_color_hue(5)[4]) +
  geom_ribbon(data = ggdata1, 
              aes(x = x, ymin = lower, ymax = upper), 
              color = gg_color_hue(5)[5], linetype = 2, 
              fill = gg_color_hue(5)[5], alpha = 0.1) +
  geom_path(data = ggdata1, 
            aes(x = x, y = `HT Pred Mean`), color = gg_color_hue(5)[5]) +
  geom_point(data = ggdata_pts, mapping = aes(x = x, y = y), 
             inherit.aes = FALSE, color = ggdata_pts$color, 
             shape = ggdata_pts$shape, 
             size = 3) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  labs(y = "y", x = "x", fill = "Function", color = "Function") +
  xlim(x_input - pT, x_input + pT) + 
  geom_vline(xintercept = pT * (0:floor((xmax - xmin) / pT)), 
             color = "blue", alpha = 0.5)
   # xlim(min(c(x_input, x_new), na.rm = TRUE),
   # max(c(x_input, x_new), na.rm = TRUE))


# plot(plotGP(
#   x_seq, y_seq, x.data = c(x_input, x_new), y.data = c(y_input, y_new), 
#   kernel = typeT, length.scale = lT, xmin = 0, xmax = 1, signal.var = 1, 
#   measurement.var = nugget))
