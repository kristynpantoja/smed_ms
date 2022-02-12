################################################################################
# last updated: 01/24/2022
# purpose: investigate seqmed with large alpha
# scenario 1:
#   linear vs. quadratic,
#   where the true function is quadratic
# scenario 2:
#   linear vs. quadratic,
#   where the true function is cubic
rm(list = ls())

scenario = 1 # 1, 2

################################################################################
# Sources/Libraries
################################################################################
output_dir = "lm_experiments/lm/outputs"
functions_dir = "functions"

# for seqmed design
source(paste(functions_dir, "/SeqMED.R", sep = ""))
source(paste(functions_dir, "/SeqMED_batch.R", sep = ""))
source(paste(functions_dir, "/charge_function_q.R", sep = ""))
source(paste(functions_dir, "/construct_design_matrix.R", sep = ""))
source(paste(functions_dir, "/wasserstein_distance.R", sep = ""))
source(paste(functions_dir, "/posterior_parameters.R", sep = ""))
source(paste(functions_dir, "/simulate_y.R", sep = ""))

# for generating initial data
source(paste(functions_dir, "/MMED.R", sep = ""))
source(paste(functions_dir, "/variance_marginal_y.R", sep = ""))

# for box-hill deisign
source(paste(functions_dir, "/boxhill.R", sep = ""))

# set up parallelization
library(foreach)
library(future)
library(doFuture)
library(parallel)
registerDoFuture()
nworkers = detectCores() / 2
plan(multisession, workers = nworkers)

library(rngtools)
library(doRNG)
rng.seed = 123 # 123, 2022
registerDoRNG(rng.seed)

################################################################################
# simulation settings, shared for both scenarios (linear vs. quadratic)
################################################################################

# simulations settings
numSims = 1 # 500 sims with N = 12, 1 sim with N = 100
numSeq = 100 # 12, 100
seqN = 1
Nttl = numSeq * seqN
xmin = -1
xmax = 1
numCandidates = 10^3 + 1
candidates = seq(from = xmin, to = xmax, length.out = numCandidates)
if(scenario == 1){
  if(numSeq == 100){
    sigmasq = 0.28
  } else if(numSeq == 12){
    sigmasq = 0.04
  }
} else if(scenario == 2){
  if(numSeq == 100){
    sigmasq = 0.21
  } else if(numSeq == 12){
    sigmasq = 0.038
  }
}

# shared settings
type01 = c(2, 3)
mu0 = c(0, 0)
mu1 = c(0, 0, 0)
sigmasq01 = 0.25
V0 = diag(rep(sigmasq01,length(mu0)))
V1 = diag(rep(sigmasq01,length(mu1)))
desX0 = function(x){
  n = length(x)
  return(cbind(rep(1, n), x))
}
desX1 = function(x){
  n = length(x)
  return(cbind(rep(1, n), x, x^2))
}
model0 = list(designMat = desX0, beta.mean = mu0, beta.var = V0)
model1 = list(designMat = desX1, beta.mean = mu1, beta.var = V1)

# boxhill settings
prior_probs = rep(1 / 2, 2)

################################################################################
# Scenarios
################################################################################
if(scenario == 1){
  betaT = c(-0.2, -0.4, 0.4)
  fT = function(x) betaT[1] + betaT[2] * x + betaT[3] * x^2
} else if(scenario == 2){
  betaT = c(0, -0.75, 0, 1)
  fT = function(x) betaT[1] + betaT[2] * x + betaT[3] * x^2 + betaT[4] * x^3
}

################################################################################
# run seqmed
################################################################################
seed = 111 
# 123: SeqMED10 is better
# 456: SeqMED100 is better
# 789: SeqMED 100 is slightly better
# 1947: SeqMED10 was better than everyone
# 777: SeqMED 100
# 555: SeqMED10 is better (even though SeqMED100 started off better) !
# 111: SeqMED is better, then worse, then better !!!!!!!!!!!

# seqmed10_sim = SeqMED(
#   y.in = NULL, x.in = NULL, true.function = fT,
#   model0 = model0, model1 = model1,
#   error.var = sigmasq, xmin = xmin, xmax = xmax,
#   candidates = candidates, numSeq = numSeq, seqN = seqN,
#   alpha_seq = 10, prints = TRUE, seed = seed)
# seqmed100_sim = SeqMED(
#   y.in = NULL, x.in = NULL, true.function = fT,
#   model0 = model0, model1 = model1,
#   error.var = sigmasq, xmin = xmin, xmax = xmax,
#   candidates = candidates, numSeq = numSeq, seqN = seqN,
#   alpha_seq = 100, prints = TRUE, seed = seed)
# boxhill_sim = BH_m2(NULL, NULL, prior_probs, model0, model1, Nttl,
#       candidates, fT, sigmasq, seed = seed)
# 
# saveRDS(
#   list(
#     seqmed10_sim = seqmed10_sim,
#     seqmed100_sim = seqmed100_sim,
#     boxhill_sim = boxhill_sim
#   ),
#   paste0(output_dir, "/alpha10vsalpha100_1sim", ".rds")
# )

saved_sims = readRDS(paste0(output_dir, "/alpha10vsalpha100_1sim", ".rds"))
seqmed10_sim = saved_sims$seqmed10_sim
seqmed100_sim = saved_sims$seqmed100_sim
boxhill_sim = saved_sims$boxhill_sim

################################################################################
# non-sequential designs
################################################################################
set.seed(rng.seed)

# for D-optimal design
library(AlgDesign)

# use AlgDesign to obtain Doptimal designs #

# consider linear model, y = b0 + b1 x
candidates_named = data.frame(x = candidates)
res_Fed_Doptlin = optFederov(
  ~1+x, data = candidates_named, approximate = TRUE, criterion = "D")
res_Fed_Doptlin$design

# consider quadratic model, y = b0 + b1 x + b2 x^2
res_Fed_Doptquad = optFederov(
  ~1+x+I(x^2), data = candidates_named, approximate = TRUE, criterion = "D")
res_Fed_Doptquad$design

# define designs #

# space-filling (grid)
space_filling = seq(from = xmin, to = xmax, length.out = Nttl)

# generate data
# Doptimal - linear
num_supportpts_Doptlin = nrow(res_Fed_Doptlin$design)
supportpt_assgnmt_Doptlin = cut(
  sample(1:Nttl, size = Nttl, replace = FALSE), # shuffle
  breaks = num_supportpts_Doptlin, labels = FALSE)
dopt_linear = rep(NA, Nttl)
for(i in 1:num_supportpts_Doptlin){
  dopt_linear[supportpt_assgnmt_Doptlin == i] = 
    res_Fed_Doptlin$design[i, "x"]
}
# # check:
# res_Fed_Doptlin$design
# table(dopt_linear) / Nttl

# Doptimal - quadratic
num_supportpts_Doptquad = nrow(res_Fed_Doptquad$design)
supportpts_Doptquad = res_Fed_Doptquad$design[, "x"]
supportpt_assgnmt_Doptquad = cut(
  sample(1:Nttl, size = Nttl, replace = FALSE), # shuffle
  breaks = num_supportpts_Doptquad, labels = FALSE)
dopt_quadratic = rep(NA, Nttl)
for(i in 1:num_supportpts_Doptquad){
  dopt_quadratic[supportpt_assgnmt_Doptquad == i] = 
    supportpts_Doptquad[i]
}
# # check:
# res_Fed_Doptquad$design
# table(dopt_quadratic) / Nttl

# half space-filling, half quadratic Doptimal, 
#   assumes Nttl is divisible by 2
supportpt_assgnmt_hybrid = cut(
  sample(1:(Nttl / 2), size = Nttl / 2, replace = FALSE), # shuffle
  breaks = num_supportpts_Doptquad, labels = FALSE)
hybrid_grid_doptq = rep(NA, Nttl / 2)
for(i in 1:num_supportpts_Doptquad){
  hybrid_grid_doptq[supportpt_assgnmt_hybrid == i] = 
    supportpts_Doptquad[i]
}
hybrid_grid_doptq[(Nttl / 2 + 1):Nttl] =c(
  seq(
    from = supportpts_Doptquad[1], to = supportpts_Doptquad[2],
    length.out = (Nttl / 4) + 2)[-c(1, (Nttl / 4) + 2)],
  seq(
    from = supportpts_Doptquad[2], to = supportpts_Doptquad[3],
    length.out = (Nttl / 4) + 2)[-c(1, (Nttl / 4) + 2)]
)
hybrid_grid_doptq = rev(hybrid_grid_doptq) # spacefilling -> doptimal
# simulations #
space_filling.tmp = sample(space_filling, replace = FALSE)
grid_sim= list(
  x = space_filling.tmp,
  y = sapply(space_filling.tmp, FUN = function(x) simulateY_fromfunction(
    x = x, true.function = fT, error.var = sigmasq)))
dopt_linear.tmp = sample(dopt_linear, replace = FALSE)
doptlin_sim = list(
  x = dopt_linear.tmp,
  y = sapply(dopt_linear.tmp, FUN = function(x) simulateY_fromfunction(
    x = x, true.function = fT, error.var = sigmasq)))
dopt_quadratic.tmp = sample(dopt_quadratic, replace = FALSE)
doptquad_sim = list(
  x = dopt_quadratic.tmp,
  y = sapply(dopt_quadratic.tmp, FUN = function(x) simulateY_fromfunction(
    x = x, true.function = fT, error.var = sigmasq)))
hybrid_grid_doptq.tmp = sample(hybrid_grid_doptq, replace = FALSE)
hybrid_sim = list(
  x = hybrid_grid_doptq.tmp,
  y = sapply(hybrid_grid_doptq.tmp, FUN = function(x) simulateY_fromfunction(
    x = x, true.function = fT, error.var = sigmasq)))

################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
# PLOTS #
################################################################################
################################################################################
################################################################################
################################################################################
 ################################################################################

# for plots
library(mvtnorm)
library(ggplot2)
library(reshape2)
library(ggpubr)
library(tidyverse)
gg_color_hue = function(n) {
  hues = seq(15, 275, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

include_hybrid = FALSE

################################################################################
# plot a seqmed with alpha = 10 or 100
sm = seqmed10_sim # seqmed100_sim
ggdata = data.frame(x = c(sm$x.in, sm$x.new), y = c(sm$y.in, sm$y.new))
plt0 = ggplot(ggdata) +
  geom_histogram(binwidth = 0.12, closed = "right",
                 aes(x = x, y = after_stat(density))) +
  theme_bw() + #base_size = 20) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
plt1 = ggplot(ggdata) +
  geom_point(aes(x, y), col = gg_color_hue(2)[1]) +
  stat_function(fun = fT) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# ggarrange(plt0, plt1)

# plot wasserstein distance
numseq = 1e2
x_seq = seq(from = xmin, to = xmax, length.out = numseq)
w_seq = sapply(x_seq, function(x) WNlm(
  x, sm$postmean0[, Nttl], sm$postmean1[, Nttl],
  diag(sm$postvar0[, Nttl]), diag(sm$postvar1[, Nttl]), 
  model0, model1, sigmasq))
f1est = function(x) sm$postmean1[1, Nttl] +
  sm$postmean1[2, Nttl] * x + sm$postmean1[3, Nttl] * x^2
f2est = function(x) sm$postmean0[1, Nttl] +
  sm$postmean0[2, Nttl] * x
f1est_seq = sapply(x_seq, f1est)
f2est_seq = sapply(x_seq, f2est)
fT_seq = sapply(x_seq, fT)

ggdata = data.table::data.table(
  x = x_seq,
  `Est. Quadr.` = f1est_seq,
  `Est. Line` = f2est_seq,
  `True Quadr.` = fT_seq,
  `Wasserstein` = w_seq
)
ggdata = data.table::melt(
  ggdata, id = c("x"), value.name = "y", variable.name = "Function")

ggdata_ribbon = data.table::data.table(
  x = x_seq,
  ymin = apply(cbind(f1est_seq, f2est_seq), 1, min),
  ymax = apply(cbind(f1est_seq, f2est_seq), 1, max)
)
pltw = ggplot(
  ggdata, aes(x = x, y = y, color = Function, linetype = Function)) +
  scale_linetype_manual(values = c(2, 2, 1, 1)) +
  ylim(-1.9, 1.9) +
  scale_color_manual(
    values = c(gg_color_hue(4)[c(3, 4)], "black", gg_color_hue(4)[2])) +
  geom_path() +
  geom_ribbon(
    data = ggdata_ribbon, mapping = aes(x = x, ymin = ymin, ymax = ymax),
    alpha = 0.2, inherit.aes = FALSE) +
  theme_bw() +
  theme(panel.grid.minor = element_blank())
pltw
ggarrange(plt0, plt1, pltw, nrow = 1, ncol = 3, widths = c(1, 1, 1.75))

################################################################################
# plot the posterior probabilities of the hypotheses
################################################################################
include_hybrid = FALSE

# functions ####################################################################
if(scenario == 1){
  models = list(model0, model1)
} else if(scenario == 2){
  desXT = function(x){
    n = length(x)
    return(cbind(rep(1, n), x, x^2, x^3))
  }
  modelT = list(
    designMat = desXT, 
    beta.mean = rep(0, 4), 
    beta.var = diag(rep(sigmasq01, 4)), 4)
  models = list(model0, model1, modelT)
}

getPPH = function(
  design, models, true.function, error.var, initial.data = FALSE, seed = NULL
){
  # format the data
  if(!is.null(seed)) set.seed(seed)
  if(initial.data){
    x = c(design$x.in, design$x.new)
    y = c(design$y.in, design$y.new)
  } else{
    x = design$x
    y = design$y
  }
  
  # get model evidences
  model.evidences = rep(NA, length(models))
  for(m in 1:length(models)){
    model.tmp = models[[m]]
    model.evidences[m] = Evidence_lm(
      y = y, x = x, model = model.tmp, error.var = error.var)
  }
  # get each hypotheses' posterior probability
  PPHs = getHypothesesPosteriors(
    prior.probs = rep(1 / length(models), length(models)), 
    evidences = model.evidences)
  PPHs = as.data.frame(t(matrix(PPHs)))
  row.names(PPHs) = NULL
  names(PPHs) = paste("H", 0:(length(models) - 1), sep = "")
  return(PPHs)
}

getPPHseq = function(
  design, models, n, true.function, error.var, initial.data = TRUE, 
  randomize.order = FALSE, randomize.halves.order = FALSE, seed = NULL
){
  # format the data and randomize if desired
  if(!is.null(seed)) set.seed(seed)
  if(initial.data){
    x.new = design$x.new
    if(randomize.order){
      new.order = sample(1:n, n, replace = FALSE)
      x.new = x.new[new.order]
    }
    if(randomize.halves.order){
      new.order = c(
        sample(1:(n / 2), n / 2, replace = FALSE), 
        sample(((n / 2) + 1):n, n / 2, replace = FALSE)
      )
      x.new = x.new[new.order]
    }
    x = c(design$x.in, x.new)
    y = c(design$y.in, design$y.new)
  } else{
    x = design$x
    y = design$y
    if(randomize.order){
      new.order = sample(1:n, n, replace = FALSE)
      x = x[new.order]
      y = y[new.order]
    }
    if(randomize.halves.order){
      new.order = c(
        sample(1:(n / 2), n / 2, replace = FALSE), 
        sample(((n / 2) + 1):n, n / 2, replace = FALSE)
      )
      x = x[new.order]
      y = y[new.order]
    }
  }
  
  # calculate posterior probabilities of hypothesized models
  PPH_mat = matrix(NA, nrow = n, ncol = length(models))
  for(i in 1:n){
    x.tmp = x[1:i]
    y.tmp = y[1:i]
    
    # get model evidences
    model.evidences.tmp = rep(NA, length(models))
    for(m in 1:length(models)){
      model.tmp = models[[m]]
      model.evidences.tmp[m] = Evidence_lm(
        y = y.tmp, x = x.tmp, model = model.tmp, error.var = error.var)
    }
    PPH_mat[i, ] = getHypothesesPosteriors(
      prior.probs = rep(1 / length(models), length(models)), 
      evidences = model.evidences.tmp)
  }
  colnames(PPH_mat) = paste("H", 0:(length(models) - 1), sep = "")
  PPH_mat = data.frame(PPH_mat)
  PPH_mat$index = 1:n
  return(PPH_mat)
}

#

################################################################################
# sequential & non-sequential designs' epph for the given alpha value

# non-sequential designs
# sequence of PPHs for each design
PPH_grid = getPPH(
  grid_sim, models, fT, sigmasq)
PPH_doptl = getPPH(
  doptlin_sim, models, fT, sigmasq)
PPH_doptq = getPPH(
  doptquad_sim, models, fT, sigmasq)
PPH_hybrid = getPPH(
  hybrid_sim, models, fT, sigmasq)
# master data frame
PPH_grid$Design = "Grid"
PPH_doptl$Design = "DOptLin."
PPH_doptq$Design = "DOptQuadr."
PPH_hybrid$Design = "Hybrid"
PPH_df = rbind(PPH_grid, PPH_doptl, PPH_doptq, PPH_hybrid)

PPHmean = aggregate(
  PPH_df[, names(PPH_df)[1:length(models)]],
  by = list(PPH_df[, "Design"]), FUN = function(x) mean(x, na.rm = TRUE))
names(PPHmean)[1] = "Design"
PPHmean$index = Nttl
# but we want a line, so allow interpolation by setting PPHmean$index = 0 too
PPHmean2 = PPHmean
PPHmean2$index = 0
PPHmean = rbind(PPHmean, PPHmean2)

# sequential designs
# sequence of PPHs for each design
PPH_seq.sm10 = getPPHseq(seqmed10_sim, models, Nttl, fT, sigmasq)
PPH_seq.sm100 = getPPHseq(seqmed100_sim, models, Nttl, fT, sigmasq)
PPH_seq.bh = getPPHseq(boxhill_sim, models, Nttl, fT, sigmasq)
# master data frame
PPH_seq.sm10$Design = "SeqMED 10"
PPH_seq.sm100$Design = "SeqMED 100"
PPH_seq.bh$Design = "BoxHill"
PPH_seq = rbind(PPH_seq.sm10, PPH_seq.sm100, PPH_seq.bh)

PPHmean_seq = aggregate(
  PPH_seq[, names(PPH_seq)[1:length(models)]],
  by = list(PPH_seq[, "Design"], PPH_seq[, "index"]),
  FUN = function(x) mean(x, na.rm = TRUE))
names(PPHmean_seq)[c(1, 2)] = c("Design", "index")

PPHmean_gg = rbind(PPHmean, PPHmean_seq)
PPHmean_gg = reshape2::melt(PPHmean_gg, id.vars = c("Design", "index"),
                            measure.vars = paste0("H", 0:(length(models) - 1), sep = ""),
                            variable.name = "hypothesis")
design_names = c(
  "Hybrid", "Grid", "DOptLin.", "DOptQuadr.", "BoxHill",
  "SeqMED 10", "SeqMED 100")
PPHmean_gg$Design = factor(PPHmean_gg$Design, levels = design_names)
if(scenario == 1){
  PPHmean_gg$hypothesis = factor(
    PPHmean_gg$hypothesis,
    levels = paste0("H", 0:(length(models) - 1), sep = ""),
    labels = paste0("Case ", scenario, ", H", 0:(length(models) - 1), sep = ""))
} else if(scenario == 2){
  PPHmean_gg$hypothesis = factor(
    PPHmean_gg$hypothesis,
    levels = paste0("H", 0:(length(models) - 1), sep = ""),
    labels = paste0("Case ", scenario, ", H", c(0, 1, "T"), sep = ""))
}
PPHmean_gg2 = PPHmean_gg[PPHmean_gg$index == Nttl, ]
PPHmean_gg2$Design = factor(PPHmean_gg2$Design, levels = design_names)
if(include_hybrid){
  num_fixed_designs = 4
} else{
  num_fixed_designs = 3
}
if(include_hybrid){
  epph.plt = ggplot(PPHmean_gg, aes(x = index, y = value, color = Design,
                                    linetype = Design, shape = Design)) +
    facet_wrap(~hypothesis) +
    geom_path() +
    scale_linetype_manual(values=c(
      rep("dashed", num_fixed_designs), rep("solid", 3))) +
    geom_point(data = PPHmean_gg2,
               mapping = aes(x = index, y = value, color = Design),
               inherit.aes = FALSE) +
    theme_bw() +
    ylim(0, 1) +
    labs(x = "Stage Index", y = element_blank())
  epph.plt
} else{
  PPHmean_gg_nohybrid = PPHmean_gg %>% filter(Design != "Hybrid")
  PPHmean_gg2_nohybrid = PPHmean_gg2 %>% filter(Design != "Hybrid")
  epph.plt = ggplot(
    PPHmean_gg_nohybrid, aes(x = index, y = value, color = Design,
                             linetype = Design, shape = Design)) +
    facet_wrap(~hypothesis) +
    geom_path() +
    scale_linetype_manual(
      values=c(rep("dashed", num_fixed_designs), rep("solid", 3))) +
    geom_point(
      data = PPHmean_gg2_nohybrid,
      mapping = aes(x = index, y = value, color = Design),
      inherit.aes = FALSE) +
    theme_bw() +
    ylim(0, 1) +
    labs(x = "Stage Index", y = element_blank())
  epph.plt
}

################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
# PLOTS -- one point at a time #
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
# uptoN = 11
for(uptoN in 2:Nttl){
  seqmed10_x = c(seqmed10_sim$x.in, seqmed10_sim$x.new)[1:uptoN]
  seqmed10_y = c(seqmed10_sim$y.in, seqmed10_sim$y.new)[1:uptoN]
  seqmed100_x = c(seqmed100_sim$x.in, seqmed100_sim$x.new)[1:uptoN]
  seqmed100_y = c(seqmed100_sim$y.in, seqmed100_sim$y.new)[1:uptoN]
  
  # E[P(H | y, X) | D]
  PPHmean_gg_uptoN = PPHmean_gg %>% filter(
    Design %in% c("BoxHill", "SeqMED 10", "SeqMED 100")) %>%
    filter(index %in% c(1:uptoN))
  PPHmean_gg2_uptoN = PPHmean_gg2 %>% filter(
    Design %in% c("BoxHill", "SeqMED 10", "SeqMED 100")) %>%
    filter(index %in% c(1:uptoN))
  epph.plt_uptoN = ggplot(
    PPHmean_gg_uptoN, aes(x = index, y = value, color = Design,
                          linetype = Design, shape = Design)) +
    facet_wrap(~hypothesis) +
    geom_path() +
    geom_point(
      data = PPHmean_gg2_uptoN,
      mapping = aes(x = index, y = value, color = Design),
      inherit.aes = FALSE) +
    theme_bw() +
    ylim(0, 1) +
    labs(x = "Stage Index", y = element_blank()) + 
    labs(title = paste0("N = ", uptoN))
  
  # Histogram
  seqmed_designpts = rbind(
    data.frame(
      x = seqmed10_x, 
      y = seqmed10_y, 
      Design = "SeqMED 10", 
      index = 1:uptoN, 
      alpha = seqmed10_sim$alpha_seq[1:uptoN], 
      specified_alpha = "10", 
      point_color = c(rep("Prev. Pts", uptoN - 1), "Newest Pt")), 
    data.frame(
      x = seqmed100_x, 
      y = seqmed100_y, 
      Design = "SeqMED 100", 
      index = 1:uptoN, 
      alpha = seqmed100_sim$alpha_seq[1:uptoN], 
      specified_alpha = "100",
      point_color = c(rep("Prev. Pts", uptoN - 1), "Newest Pt"))
  )
  hist_designpts = ggplot(seqmed_designpts) +
    geom_histogram(binwidth = 0.12, closed = "right",
                   aes(x = x, y = after_stat(count / max(count)))) +
    facet_wrap(vars(Design)) +
    theme_bw() + 
    geom_text(
      data = seqmed_designpts %>% filter(index == uptoN), 
      mapping = aes(x = x, y = 0, label = round(x, 5)), color = 4) +
    geom_text(
      data = seqmed_designpts %>% filter(index == uptoN), 
      mapping = aes(x = 0, y = 1, label = paste0("alpha=", alpha)), color = 2) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  # ggarrange(epph.plt_uptoN, hist_designpts, nrow = 2)
  
  # Model Fits & Wasserstein Distance
  numseq = 1e2
  x_seq = seq(from = xmin, to = xmax, length.out = numseq)
  w_seq = sapply(x_seq, function(x) WNlm(
    x, sm$postmean0[, uptoN], sm$postmean1[, uptoN],
    diag(sm$postvar0[, uptoN]), diag(sm$postvar1[, uptoN]), 
    model0, model1, sigmasq))
  f0est_seq = as.numeric(model0$designMat(x_seq) %*% sm$postmean0[, uptoN])
  f1est_seq = as.numeric(model1$designMat(x_seq) %*% sm$postmean1[, uptoN])
  fT_seq = sapply(x_seq, fT)
  
  ggdata = data.table::data.table(
    x = x_seq,
    `Est. Quadr.` = f1est_seq,
    `Est. Line` = f0est_seq,
    `True Quadr.` = fT_seq,
    `Wasserstein` = w_seq
  )
  ggdata = data.table::melt(
    ggdata, id = c("x"), value.name = "y", variable.name = "Function")
  
  ggdata_ribbon = data.table::data.table(
    x = x_seq,
    ymin = apply(cbind(f0est_seq, f1est_seq), 1, min),
    ymax = apply(cbind(f0est_seq, f1est_seq), 1, max)
  )
  ggdata2 = rbind(
    ggdata %>% mutate(specified_alpha = "10"), 
    ggdata %>% mutate(specified_alpha = "100")
  )
  pltw = ggplot(
    ggdata2, aes(x = x, y = y, color = Function, linetype = Function)) +
    facet_wrap(vars(specified_alpha)) +
    scale_linetype_manual(values = c(2, 2, 1, 1), guide = "none") +
    # ylim(-1.9, 1.9) +
    scale_color_manual(
      values = c(gg_color_hue(4)[c(3, 4)], "red", "orange", "gray50", gg_color_hue(4)[2])) + #gg_color_hue(6)[c(3, 4)])) +
    geom_path() +
    geom_ribbon(
      data = ggdata_ribbon, mapping = aes(x = x, ymin = ymin, ymax = ymax),
      alpha = 0.2, inherit.aes = FALSE) +
    # geom_text(
    #   data = seqmed_designpts %>% filter(index == uptoN), 
    #   mapping = aes(x = x, y = 0, label = round(x, 5)), color = 4, 
    #   inherit.aes = FALSE) +
    # geom_text(
    #   data = seqmed_designpts %>% filter(index == uptoN), 
    #   mapping = aes(x = 0, y = 1, label = paste0("alpha=", alpha)), color = 2, 
    #   inherit.aes = FALSE) +
    geom_point(
      data = seqmed_designpts %>% mutate(point_color = factor(point_color)), 
      mapping = aes(x = x, y = y, color = point_color), alpha = 0.5,
      inherit.aes = FALSE) +
    theme_bw() +
    theme(panel.grid.minor = element_blank())
  ggarrange(epph.plt_uptoN, pltw, nrow = 2)
  
  ggsave(
    filename = paste0(
      "lm", "_scen", scenario,
      "_epph_fits_alpha10vs100_N", uptoN, ".pdf"),
    plot = last_plot(),
    width = 8, height = 6, units = c("in")
  )
}











