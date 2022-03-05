rm(list = ls())

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
rng.seed = 1947 # 123, 345, 1947
registerDoRNG(rng.seed)

library(ggplot2)
library(reshape2)
library(ggpubr)
gg_color_hue = function(n) {
  hues = seq(15, 275, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
library(data.table)

################################################################################
# Simulation settings, shared for both scenarios
################################################################################

# simulations settings
numSims = 100
numSeq = 15
seqN = 1
Nttl = numSeq * seqN
xmin = -1
xmax = 1
numx = 10^3 + 1
x_seq = seq(from = xmin, to = xmax, length.out = numx)
sigmasq_measuremt = 1e-10
sigmasq_signal = 1

# shared settings
nugget = sigmasq_measuremt
prior_probs = rep(1 / 2, 2)

################################################################################
# Scenario 4 Settings: Matern vs. SE, f ~ Periodic
################################################################################

type01 = c("matern", "squaredexponential")
typeT = "periodic"
l01= c(0.01, 0.01)
lT = 0.1
pT = 0.05 # 0.05 or 0.1
model0 = list(type = type01[1], l = l01[1], signal.var = sigmasq_signal, 
              measurement.var = nugget)
model1 = list(type = type01[2], l = l01[2], signal.var = sigmasq_signal, 
              measurement.var = nugget)

################################################################################
# Generate data
################################################################################

# generate functions 
registerDoRNG(rng.seed)
null_cov = getCov(x_seq, x_seq, typeT, lT, pT, sigmasq_signal)
null_mean = rep(0, numx)

# the function values
y_seq_mat = t(rmvnorm(n = numSims, mean = null_mean, sigma = null_cov)) 
filename_append = ""
if(is.null(sigmasq_measuremt)){
  y_seq_mat = t(rmvnorm(n = numSims, mean = null_mean, sigma = null_cov))
  filename_append = ""
} else{
  y_seq_mat = t(rmvnorm(n = numSims, mean = null_mean, sigma = null_cov + 
                          sigmasq_measuremt * diag(numx))) 
  filename_append = "_noise"
}

sim_data = list(
  x = x_seq, 
  null_mean = null_mean, 
  null_cov = null_cov, 
  numSims = numSims, 
  function_values_mat = y_seq_mat
)

################################################################################
# Run SeqMED and Box & Hill design methods
################################################################################

# generate seqmeds 
registerDoRNG(rng.seed)
seqmeds = foreach(
  b = 1:numSims
) %dorng% {
  y_seq = y_seq_mat[ , b]
  SeqMEDgp(
    y.in = NULL, x.in = NULL, x.in.idx = NULL, 
    candidates = x_seq, function.values = y_seq, xmin = xmin, xmax = xmax,
    model0 = model0, model1 = model1, numSeq = numSeq, seqN = seqN)
}

# generate boxhills
registerDoRNG(rng.seed)
boxhills = foreach(
  i = 1:numSims
) %dorng% {
  y_seq = y_seq_mat[ , i]
  BHgp_m2(
    y.in = NULL, x.in = NULL, x.in.idx =  NULL, 
    prior.probs = prior_probs, model0 = model0, model1 = model1, n = Nttl, 
    candidates = x_seq, function.values = y_seq, seed = NULL)
}

boxhills2 = list()
for(i in 1:numSims){
  boxhills2.tmp = boxhills[[i]]
  x.last.idx = max(which(!is.na(boxhills2.tmp$x.new)))
  if(x.last.idx < (Nttl - 1)){
    boxhills2.tmp$x.new.idx[(x.last.idx + 1):(Nttl - 1)] = 
      boxhills2.tmp$x.new.idx[x.last.idx]
    boxhills2.tmp$x.new[(x.last.idx + 1):(Nttl - 1)] = 
      boxhills2.tmp$x.new[x.last.idx]
    boxhills2.tmp$y.new[(x.last.idx + 1):(Nttl - 1)] = 
      boxhills2.tmp$y.new[x.last.idx]
    # posterior probabilities
    postprobs.tmp = matrix(
      boxhills2.tmp$post.probs[1, ], nrow = Nttl, ncol = 2, byrow = TRUE)
    for(j in 1:(Nttl - 1)){
      x.cur = c(boxhills2.tmp$x.in, boxhills2.tmp$x.new[1:j])
      y.cur = c(boxhills2.tmp$y.in, boxhills2.tmp$y.new[1:j])
      post.probs.cur = postprobs.tmp[j, ]
      post.probs.cur = getHypothesesPosteriors(
        prior.probs = post.probs.cur, 
        evidences = c(
          Evidence_gp(y.cur, x.cur, model0),
          Evidence_gp(y.cur, x.cur, model1)
        )
      )
      postprobs.tmp[j + 1, ] = post.probs.cur
    }
    postprobs.last.idx = max(which(
      apply(boxhills2.tmp$post.probs, 1, function(row) any(!is.na(row)))))
    if(all.equal(
      postprobs.tmp[1:postprobs.last.idx, ], 
      boxhills2.tmp$post.probs[1:postprobs.last.idx, ])){
      boxhills2.tmp$post.probs = postprobs.tmp
    } else{
      stop("issue with postprobs.tmp")
    }
  }
  boxhills2[[i]] = boxhills2.tmp
}

################################################################################
# make sequential EPPH plots
################################################################################

# models
if(typeT == "periodic"){
  modelT = list(type = typeT, l = lT, signal.var = sigmasq_signal, 
                measurement.var = sigmasq_measuremt, p = pT)
} else{
  modelT = list(type = typeT, l = lT, signal.var = sigmasq_signal, 
                measurement.var = sigmasq_measuremt)
}

getPPHseq = function(
  design, model0, model1, modelT, n, randomize.order = FALSE){
  n.new = n - 1
  x.new.idx = design$x.new.idx
  x.new = design$x.new
  y.new = design$y.new
  if(n.new != length(y.new)) warning("getPPHseq: n argument does not match length of new data")
  len.tmp = length(as.vector(na.omit(y.new)))
  if(randomize.order){
    new.order = sample(1:len.tmp, len.tmp, replace = FALSE)
    x.new.idx = x.new.idx[new.order]
    x.new = x.new[new.order]
    y.new = y.new[new.order]
  }
  # calculate posterior probs for each new point
  PPH0_seq = rep(NA, len.tmp)
  PPH1_seq = rep(NA, len.tmp)
  PPHT_seq = rep(NA, len.tmp)
  for(i in 1:len.tmp){
    y.tmp = c(design$y.in, y.new[1:i])
    x.tmp = c(design$x.in, x.new[1:i])
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
  if(length(PPH0_seq) < n.new){
    PPH0_seq[(length(PPH0_seq) + 1):n.new] = PPH0_seq[length(PPH0_seq)]
  }
  if(length(PPH1_seq) < n.new){
    PPH1_seq[(length(PPH1_seq) + 1):n.new] = PPH1_seq[length(PPH1_seq)]
  }
  if(length(PPHT_seq) < n.new){
    PPHT_seq[(length(PPHT_seq) + 1):n.new] = PPHT_seq[length(PPHT_seq)]
  }
  # include posterior probs for initial point(s)
  len.init = length(design$y.in)
  PPHs.init = getHypothesesPosteriors(
    prior.probs = rep(1 / 3, 3), 
    evidences = c(
      Evidence_gp(design$y.in, design$x.in, model0),
      Evidence_gp(design$y.in, design$x.in, model1), 
      Evidence_gp(design$y.in, design$x.in, modelT)
    )
  )
  PPH0_seq = c(PPHs.init[1], PPH0_seq)
  PPH1_seq = c(PPHs.init[2], PPH1_seq)
  PPHT_seq = c(PPHs.init[3], PPHT_seq)
  return(data.frame(
    index = 1:n, 
    "H0" = PPH0_seq, 
    "H1" = PPH1_seq, 
    "HT" = PPHT_seq
  ))
}

PPH_seq = data.frame(
  PPH0 = numeric(), PPH1 = numeric(), PPHT = numeric(), 
  Design = character(), sim = numeric())
for(j in 1:numSims){
  # designs at sim b
  bh = boxhills[[j]]
  bh2 = boxhills2[[j]]
  sm = seqmeds[[j]]
  # sequence of PPHs for each design
  PPH_seq.bh = getPPHseq(bh, model0, model1, modelT, Nttl)
  PPH_seq.bh2 = getPPHseq(bh2, model0, model1, modelT, Nttl)
  PPH_seq.sm = getPPHseq(sm, model0, model1, modelT, Nttl)
  # master data frame
  PPH_seq.bh$Design = "BoxHill"
  PPH_seq.bh2$Design = "BoxHill2"
  PPH_seq.sm$Design = "SeqMED"
  PPH_seq.tmp = rbind(
    PPH_seq.bh, PPH_seq.bh2, PPH_seq.sm)
  PPH_seq.tmp$sim = j
  PPH_seq = rbind(PPH_seq, PPH_seq.tmp)
}

PPH0mean_seq = aggregate(PPH_seq$H0, by = list(PPH_seq$index, PPH_seq$Design), 
                         FUN = function(x) mean(x, na.rm = TRUE))
names(PPH0mean_seq) = c("index", "Design", "value")
PPH0mean_seq$Hypothesis = "H0"
PPH1mean_seq = aggregate(PPH_seq$H1, by = list(PPH_seq$index, PPH_seq$Design), 
                         FUN = function(x) mean(x, na.rm = TRUE))
names(PPH1mean_seq) = c("index", "Design", "value")
PPH1mean_seq$Hypothesis = "H1"
PPHTmean_seq = aggregate(PPH_seq$HT, by = list(PPH_seq$index, PPH_seq$Design), 
                         FUN = function(x) mean(x, na.rm = TRUE))
names(PPHTmean_seq) = c("index", "Design", "value")
PPHTmean_seq$Hypothesis = "HT"

PPHmean_seq = rbind(PPH0mean_seq, PPH1mean_seq, PPHTmean_seq)
epph.plt = ggplot(PPHmean_seq, aes(x = index, y = value, color = Design, 
                                   linetype = Design, shape = Design)) + 
  facet_wrap(~Hypothesis) + 
  geom_path() + 
  geom_point() +
  theme_bw() +
  ylim(0, 1) + 
  labs(x = "Stage Index", y = element_blank()) +
  scale_x_continuous(breaks = c(5, 10, 15))
plot(epph.plt)











