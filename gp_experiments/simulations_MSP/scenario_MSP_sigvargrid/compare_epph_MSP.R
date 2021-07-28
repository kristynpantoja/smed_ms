
scenario = 4 # scenarios: 3,4,5,6
seq.type = 1 # 1 = fully sequential, 2 = stage-sequential 3x5
sigmasq_signal_seq = seq(0.5, 1.5, length.out = 101)

################################################################################
# Sources/Libraries
################################################################################
sims_dir = "gp_experiments/simulations_MSP/scenario_MSP_sigvargrid"
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

getPPHlast = function(design, model0, model1, modelT){
  len.new = length(as.vector(na.omit(design$y.new)))
    y.tmp = c(design$y, as.vector(na.omit(design$y.new))[1:len.new])
    x.tmp = c(design$x, as.vector(na.omit(design$x.new))[1:len.new])
    PPHs.tmp = getHypothesesPosteriors(
      prior.probs = rep(1 / 3, 3), 
      evidences = c(
        Evidence_gp(y.tmp, x.tmp, model0),
        Evidence_gp(y.tmp, x.tmp, model1), 
        Evidence_gp(y.tmp, x.tmp, modelT)
      )
    )
  return(data.frame(
    index = len.new, 
    "H0" = PPHs.tmp[1], 
    "H1" = PPHs.tmp[2], 
    "HT" = PPHs.tmp[3]
  ))
}

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

# shared settings
nugget = sigmasq_measuremt
prior_probs = rep(1 / 2, 2)

################################################################################
# simulations!

PPH_lasts_seq = data.table(
  type = character(), H0 = numeric(), H1 = numeric(), HT = numeric(), 
  sigvar = numeric(), index = numeric())
y_seq_sigvargrid = matrix(NA, nrow = numx, ncol = length(sigmasq_signal_seq))
for(i in 1:length(sigmasq_signal_seq)){
  ################################################################################
  # Scenario settings
  type01 = c("matern", "squaredexponential")
  typeT = "periodic"
  l01= c(0.01, 0.01)
  pT = 0.05
  lT = 0.5
  sigmasq_signal = sigmasq_signal_seq[i]
  model0 = list(type = type01[1], l = l01[1], signal.var = 1, 
                measurement.var = nugget)
  model1 = list(type = type01[2], l = l01[2], signal.var = 1, 
                measurement.var = nugget)
  
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
      "_sig", sigmasq_signal,
      filename_append, 
      "_seed", rng.seed,
      ".rds")
  } else{
    simulated_data_file = paste0(
      data_dir,
      "/", typeT,
      "_l", lT,
      "_sig", sigmasq_signal,
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
  
  y_seq_sigvargrid[, i] = y_seq_mat[, 1]
  
  ################################################################################
  # initial design
  x_input_idx = ceiling(numx / 2)
  x_input = x_seq[x_input_idx]
  
  ################################################################################
  # read in the designs
  # filename_append.tmp for all methods alike
  filename_append.tmp = paste0(
    filename_append, 
    "_seed", rng.seed,
    ".rds"
  )
  boxhill_sims = readRDS(paste0(
    output_dir,
    "/scenario", scenario, "_boxhill", 
    "_sig", sigmasq_signal,
    filename_append.tmp))
  leaveout_sims = readRDS(paste0(
    output_dir,
    "/scenario", scenario, "_seqmed", 
    "_leaveout", 
    "_seq", seq.type, 
    "_sig", sigmasq_signal,
    filename_append.tmp))
  qcap_sims = readRDS(paste0(
    output_dir,
    "/scenario", scenario, "_seqmed", 
    "_cap",
    "_seq", seq.type, 
    "_sig", sigmasq_signal,
    filename_append.tmp))
  persist_sims = readRDS(paste0(
    output_dir,
    "/scenario", scenario, "_seqmed", 
    "_persist", 
    "_seq", seq.type, 
    "_sig", sigmasq_signal,
    filename_append.tmp))
  qcap_persist_sims = readRDS(paste0(
    output_dir,
    "/scenario", scenario, "_seqmed", 
    "_cap_persist",
    "_seq", seq.type, 
    "_sig", sigmasq_signal,
    filename_append.tmp))
  
  if(typeT == "periodic"){
    random_sims_file = paste0(
      sims_dir, 
      "/spacefilling_designs/outputs/random", 
      "_", typeT,
      "_l", lT,
      "_p", pT,
      "_sig", sigmasq_signal,
      filename_append.tmp)
    grid_sims_file = paste0(
      sims_dir,
      "/spacefilling_designs/outputs/grid", 
      "_", typeT,
      "_l", lT,
      "_p", pT,
      "_sig", sigmasq_signal,
      filename_append.tmp)
  } else{
    random_sims_file = paste0(
      sims_dir, 
      "/spacefilling_designs/outputs/random", 
      "_", typeT,
      "_l", lT,
      "_sig", sigmasq_signal,
      filename_append.tmp)
    grid_sims_file = paste0(
      sims_dir,
      "/spacefilling_designs/outputs/grid", 
      "_", typeT,
      "_l", lT,
      "_sig", sigmasq_signal,
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
  PPHlast = data.frame(
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
    PPH.bh = getPPHlast(bh, model0, model1, modelT)
    PPH.qc = getPPHlast(qc, model0, model1, modelT)
    PPH.lo = getPPHlast(lo, model0, model1, modelT)
    PPH.qc2 = getPPHlast(qc2, model0, model1, modelT)
    PPH.kq = getPPHlast(kq, model0, model1, modelT)
    PPH.r = getPPHlast(r, model0, model1, modelT)
    PPH.g = getPPHlast(g, model0, model1, modelT)
    # master data frame
    PPH.bh$type = "boxhill"
    PPH.qc$type = "qcap"
    PPH.lo$type = "leaveout"
    PPH.qc2$type = "keepq2"
    PPH.kq$type = "keepq"
    PPH.r$type = "random"
    PPH.g$type = "grid"
    PPH.tmp = rbind(PPH.bh, PPH.qc, PPH.lo, PPH.qc2, PPH.kq, PPH.r, PPH.g)
    PPH.tmp$sim = j
    PPHlast = rbind(PPHlast, PPH.tmp)
  }
  PPHlast$index = NULL
  
  PPH0mean_last = aggregate(PPHlast$H0, by = list(PPHlast$type), 
                           FUN = function(x) mean(x, na.rm = TRUE))
  names(PPH0mean_last) = c("type", "value")
  PPH0mean_last$Hypothesis = "H0"
  PPH1mean_last = aggregate(PPHlast$H1, by = list(PPHlast$type), 
                           FUN = function(x) mean(x, na.rm = TRUE))
  names(PPH1mean_last) = c("type", "value")
  PPH1mean_last$Hypothesis = "H1"
  PPHTmean_last = aggregate(PPHlast$HT, by = list(PPHlast$type), 
                           FUN = function(x) mean(x, na.rm = TRUE))
  names(PPHTmean_last) = c("type", "value")
  PPHTmean_last$Hypothesis = "HT"
  
  # plot curves
  PPHmean_last = rbind(PPH0mean_last, PPH1mean_last, PPHTmean_last)
  PPHmean_last = data.table::dcast(data.table(PPHmean_last), type~Hypothesis)
  PPHmean_last[, sigvar := sigmasq_signal]
  PPHmean_last[, index := i]
  PPH_lasts_seq = rbind(PPH_lasts_seq, PPHmean_last)
}

PPH_lasts_seq_melt = data.table::melt(
  PPH_lasts_seq, id.vars = c("sigvar", "index", "type"))

ggplot(PPH_lasts_seq_melt, aes(x = sigvar, y = value, color = type)) + 
  facet_wrap(vars(variable)) + 
  geom_path()
  
plt_seq = c(1, 26, 51, 76, 101)
plt_list = list()
for(i in 1:length(plt_seq)){
  plt_list[[i]] = ggplot(
    data.frame(x = x_seq, y = y_seq_sigvargrid[, plt_seq[i]]), 
    aes(x = x, y = y)) + 
    geom_path() + 
    xlim(0.25, 0.75) + 
    ggtitle(sigmasq_signal_seq[plt_seq[i]])
}
ggpubr::ggarrange(plotlist = plt_list)

