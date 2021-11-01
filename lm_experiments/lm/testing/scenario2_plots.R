################################################################################
# last updated: 10/29/21
# purpose: to create a list of seqmed simulations for scenario 2:
#   linear vs. quadratic,
#   where the true function is cubic

scenario = 2

beta_setting = 0 # 0, 4
sigmasq = 0.1 # 0.1, 0.05
discontinuity = 0.01
height = 1
numSeq = 12 # 12

################################################################################
# Sources/Libraries
################################################################################
output_dir = "lm_experiments/lm/testing"
functions_home = "functions"

# for seqmed design
source(paste(functions_home, "/SeqMED.R", sep = ""))
source(paste(functions_home, "/SeqMED_batch.R", sep = ""))
source(paste(functions_home, "/charge_function_q.R", sep = ""))
source(paste(functions_home, "/construct_design_matrix.R", sep = ""))
source(paste(functions_home, "/wasserstein_distance.R", sep = ""))
source(paste(functions_home, "/posterior_parameters.R", sep = ""))
source(paste(functions_home, "/simulate_y.R", sep = ""))

# for generating initial data
source(paste(functions_home, "/MMED.R", sep = ""))
source(paste(functions_home, "/variance_marginal_y.R", sep = ""))

# for box-hill deisign
source(paste(functions_home, "/boxhill.R", sep = ""))

# for D-optimal design
library(AlgDesign)

# set up parallelization
library(doFuture)
library(parallel)
registerDoFuture()
nworkers = detectCores()
plan(multisession, workers = nworkers)

library(mvtnorm)
library(ggplot2)
library(reshape2)
library(ggpubr)
gg_color_hue = function(n) {
  hues = seq(15, 275, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

library(doRNG)
registerDoRNG(1995)

################################################################################
# simulation settings, shared for both scenarios (linear vs. quadratic)
################################################################################

# simulations settings
numSims = 100 #100
# numSeq = 12 #100, 36
seqN = 1
Nttl = numSeq * seqN
xmin = -1
xmax = 1
numCandidates = 10^3 + 1
candidates = seq(from = xmin, to = xmax, length.out = numCandidates)
# sigmasq = 0.1 # 0.1

# shared settings
type01 = c(2, 3)
mu0 = c(0, 0)
mu1 = c(0, 0, 0)
sigmasq01 = 0.25
V0 = diag(rep(sigmasq01,length(mu0)))
V1 = diag(rep(sigmasq01,length(mu1)))
f0 = function(x) mu0[1] + mu0[2] * x
f1 = function(x) mu1[1] + mu1[2] * x + mu1[3] * x^2
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
# Scenario 2: True function is cubic
################################################################################
if(beta_setting == 0){
  betaT = c(0, -0.75, 0, 1)
  fT = function(x) betaT[1] + betaT[2] * x + betaT[3] * x^2 + betaT[4] * x^3
} else if(beta_setting == 1){
  betaT = c(0, 0, 0, 1)
  fT = function(x) betaT[1] + betaT[2] * x + betaT[3] * x^2 + betaT[4] * x^3
} else if(beta_setting == 2){
  betaT = c(0, 0, 1, 1)
  fT = function(x) betaT[1] + betaT[2] * x + betaT[3] * x^2 + betaT[4] * x^3
} else if(beta_setting == 3){
  betaT = c(0, 0, 1)
  fx_outer = function(x) betaT[1] + betaT[2] * x + betaT[3] * x^2
  # discontinuity = 0.01
  # height = 1
  fx_innerneg = function(x) height - ((fx_outer(discontinuity) - height) / discontinuity) * x
  fx_innerpos = function(x) height + ((fx_outer(discontinuity) - height) / discontinuity) * x
  fT = function(x){
    fx = rep(NA, length(x))
    for(i in 1:length(x)){
      if(abs(x[i]) >= discontinuity){
        fx[i] = fx_outer(x[i])
      } else if(x[i] > -discontinuity & x[i] <= 0){
        fx[i] = fx_innerneg(x[i])
      } else{
        fx[i] = fx_innerpos(x[i])
      }
    }
    return(fx)
  }
} else if(beta_setting == 4){
  betaT = c(0, 0, 1)
  fx_outer = function(x) betaT[1] + betaT[2] * x + betaT[3] * x^2
  # discontinuity = 0.01
  # height = 1
  fx_inner = function(x) height
  fT = function(x){
    fx = rep(NA, length(x))
    for(i in 1:length(x)){
      if(abs(x[i]) >= discontinuity){
        fx[i] = fx_outer(x[i])
      } else{
        fx[i] = fx_inner(x[i])
      }
    }
    return(fx)
  }
}
ggplot() + 
  xlim(xmin, xmax) + 
  geom_function(fun = fT) + 
  theme_bw()
# ggsave(
#   filename = paste0(
#     "lm_", "_scen", scenario,
#     "_beta", beta_setting,
#     "_curve", ".pdf"),
#   plot = last_plot(),
#   width = 6.5, height = 3.5, units = c("in")
# )

################################################################################
# import sequential designs
################################################################################

if(beta_setting == 0){
  seqmed_sims = readRDS(file = paste0(
    output_dir,
    "/sm", "_scen", scenario, 
    "_beta", beta_setting, 
    "_N", Nttl, 
    "_sigmasq", sigmasq,
    "_numSims", numSims, 
    ".rds"
  ))
}
if(beta_setting == 4){
  seqmed_sims = readRDS(file = paste0(
    output_dir,
    "/sm", "_scen", scenario, 
    "_beta", beta_setting, 
    "_height", height, 
    "_discontinuity", discontinuity,
    "_N", Nttl, 
    "_sigmasq", sigmasq,
    "_numSims", numSims, 
    ".rds"
  ))
}

################################################################################
# non-sequential designs
################################################################################

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
supportpt_assgnmt_Doptquad = cut(
  sample(1:Nttl, size = Nttl, replace = FALSE), # shuffle
  breaks = num_supportpts_Doptquad, labels = FALSE)
dopt_quadratic = rep(NA, Nttl)
for(i in 1:num_supportpts_Doptquad){
  dopt_quadratic[supportpt_assgnmt_Doptquad == i] = 
    res_Fed_Doptquad$design[i, "x"]
}
# # check:
# res_Fed_Doptquad$design
# table(dopt_quadratic) / Nttl

# half space-filling, half quadratic Doptimal
supportpt_assgnmt_hybrid = cut(
  sample(1:(Nttl / 2), size = Nttl / 2, replace = FALSE), # shuffle
  breaks = num_supportpts_Doptquad, labels = FALSE)
hybrid_grid_doptq = rep(NA, Nttl / 2)
for(i in 1:num_supportpts_Doptquad){
  hybrid_grid_doptq[supportpt_assgnmt_hybrid == i] = 
    res_Fed_Doptquad$design[i, "x"]
}
hybrid_grid_doptq[(Nttl / 2 + 1):Nttl] = seq(
  from = xmin, to = xmax, length.out = Nttl / 2)

grid_sims = list()
doptlin_sims = list()
doptquad_sims = list()
hybrid_sims = list()
for(i in 1:numSims){
  grid_sims[[i]] = list(
    x = space_filling,
    y = sapply(space_filling, FUN = function(x) simulateY_fromfunction(
      x = x, true.function = fT, error.var = sigmasq)))
  doptlin_sims[[i]] = list(
    x = dopt_linear,
    y = sapply(dopt_linear, FUN = function(x) simulateY_fromfunction(
      x = x, true.function = fT, error.var = sigmasq)))
  doptquad_sims[[i]] = list(
    x = dopt_quadratic,
    y = sapply(dopt_quadratic, FUN = function(x) simulateY_fromfunction(
      x = x, true.function = fT, error.var = sigmasq)))
  hybrid_sims[[i]] = list(
    x = hybrid_grid_doptq,
    y = sapply(hybrid_grid_doptq, FUN = function(x) simulateY_fromfunction(
      x = x, true.function = fT, error.var = sigmasq)))
}



################################################################################
################################################################################
################################################################################
# plots!!!
################################################################################
################################################################################
################################################################################



################################################################################
# plot the designs
################################################################################
sim.idx = 1

# plot a seqmed
sm = seqmed_sims[[sim.idx]]
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

################################################################################
# plot the wasserstein distance when scenario == 1
# plot the posterior mean curve and the true curve when scenario == 2
################################################################################
library(data.table)

if(scenario == 1){
  numseq = 1e2
  x_seq = seq(from = xmin, to = xmax, length.out = numseq)
  w_seq = sapply(x_seq, function(x) WNlm(
    x, sm$postmean0, sm$postmean1,
    diag(sm$postvar0), diag(sm$postvar1), model0, model1, sigmasq))
  f1est = function(x) sm$postmean1[1, ] +
    sm$postmean1[2, ] * x + sm$postmean1[3, ] * x^2
  f2est = function(x) sm$postmean0[1, ] +
    sm$postmean0[2, ] * x
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
    ylim(-1.1, 2.1) +
    scale_color_manual(
      values = c(gg_color_hue(4)[c(3, 4)], "black", gg_color_hue(4)[2])) +
    geom_path() +
    geom_ribbon(
      data = ggdata_ribbon, mapping = aes(x = x, ymin = ymin, ymax = ymax), 
      alpha = 0.2, inherit.aes = FALSE) +
    theme_bw() +
    theme(panel.grid.minor = element_blank())
  # ggarrange(plt0, plt1, pltw, nrow = 1, ncol = 3, widths = c(1, 1, 1.75))
  
  # # manuscript plot
  # ggsave(
  #   filename = paste0("scen", scenario, "_seqmedexamplewasserstein.pdf"),
  #   plot = last_plot(),
  #   width = 6.5, height = 1.75, units = c("in")
  # )
} else if(scenario == 2){
  # plot ...
}

################################################################################
# plot the posterior probabilities of the hypotheses
################################################################################
library(data.table)

if(scenario == 1){
  models = list(model0, model1)
} else if(scenario == 2){
  if(beta_setting %in% 0:2){
    desXT = function(x){
      n = length(x)
      return(cbind(rep(1, n), x, x^2, x^3))
    }
    modelT = list(
      designMat = desXT, 
      beta.mean = rep(0, 4), 
      beta.var = diag(rep(sigmasq01, 4)), 4)
    models = list(model0, model1, modelT)
  } else if(beta_setting == 3){
    xi1 = -discontinuity
    xi2 = 0
    xi3 = discontinuity
    h1 = function(x) 1
    h2 = function(x) x
    h3 = function(x) x^2
    h4 = function(x) ifelse(x - xi1 > 0, x - xi1, 0)
    h5 = function(x) ifelse(x - xi2 > 0, x - xi2, 0)
    h6 = function(x) ifelse(x - xi3 > 0, x - xi3, 0)
    desXT = function(x){
      cbind(h1(x), h2(x), h3(x), h4(x), h5(x), h6(x))
    }
    modelT = list(
      designMat = desXT, 
      beta.mean = rep(0, 6), 
      beta.var = diag(rep(sigmasq01, 6)), 6)
    models = list(model0, model1, modelT)
  } else if(beta_setting == 4){
    xi1 = -discontinuity
    xi2 = discontinuity
    h1 = function(x) 1
    h2 = function(x) x
    h3 = function(x) x^2
    h4 = function(x) ifelse(x < xi1, 1, 0)
    h5 = function(x) ifelse(x < xi2, 1, 0)
    desXT = function(x){
      cbind(h1(x), h2(x), h3(x), h4(x), h5(x))
    }
    modelT = list(
      designMat = desXT, 
      beta.mean = rep(0, 5), 
      beta.var = diag(rep(sigmasq01, 5)), 5)
    models = list(model0, model1, modelT)
  }
}

getPPH = function(
  design, models, true.function, error.var, initial.data = FALSE, seed = NULL
){
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
  randomize.order = FALSE, seed = NULL
){
  if(!is.null(seed)) set.seed(seed)
  if(initial.data){
    x.new = design$x.new
    if(randomize.order){
      new.order = sample(1:n, n, replace = FALSE)
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

################################################################################

# non-sequential designs
PPH_df = data.frame()
for(j in 1:numSims){
  # sequence of PPHs for each design
  PPH_grid = getPPH(
    grid_sims[[j]], models, fT, sigmasq)
  PPH_doptl = getPPH(
    doptlin_sims[[j]], models, fT, sigmasq)
  PPH_doptq = getPPH(
    doptquad_sims[[j]], models, fT, sigmasq)
  PPH_hybrid = getPPH(
    hybrid_sims[[j]], models, fT, sigmasq)
  # master data frame
  PPH_grid$Design = "Grid"
  PPH_doptl$Design = "DOptLin."
  PPH_doptq$Design = "DOptQuadr."
  PPH_hybrid$Design = "Hybrid"
  PPH.tmp = rbind(PPH_grid, PPH_doptl, PPH_doptq, PPH_hybrid)
  PPH.tmp$sim = j
  PPH_df = rbind(PPH_df, PPH.tmp)
}
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
PPH_seq = data.frame()
for(j in 1:numSims){
  # sequence of PPHs for each design
  PPH_seq.sm = getPPHseq(seqmed_sims[[j]], models, Nttl, fT, sigmasq)
  # master data frame
  PPH_seq.sm$Design = "SeqMED"
  PPH_seq.tmp =PPH_seq.sm
  PPH_seq.tmp$sim = j
  PPH_seq = rbind(PPH_seq, PPH_seq.tmp)
}

PPHmean_seq = aggregate(
  PPH_seq[, names(PPH_seq)[1:length(models)]],
  by = list(PPH_seq[, "Design"], PPH_seq[, "index"]),
  FUN = function(x) mean(x, na.rm = TRUE))
names(PPHmean_seq)[c(1, 2)] = c("Design", "index")

PPHmean_gg = rbind(PPHmean, PPHmean_seq)
PPHmean_gg = melt(PPHmean_gg, id.vars = c("Design", "index"),
                  measure.vars = paste0("H", 0:(length(models) - 1), sep = ""),
                  variable.name = "hypothesis")
design_names = rev(c("SeqMED", "DOptLin.", "DOptQuadr.", "Grid", "Hybrid"))
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
PPHmean_gg = setorder(PPHmean_gg, cols = "Design")
PPHmean_gg2 = PPHmean_gg[PPHmean_gg$index == Nttl, ]
PPHmean_gg2$Design = factor(PPHmean_gg2$Design, levels = design_names)
PPHmean_gg2 = setorder(PPHmean_gg2, cols = "Design")
epph.plt = ggplot(PPHmean_gg, aes(x = index, y = value, color = Design,
                                  linetype = Design, shape = Design)) +
  facet_wrap(~hypothesis) +
  geom_path() +
  # scale_linetype_manual(values=c(rep("dashed", 4), rep("solid", 2))) +
  geom_point(data = PPHmean_gg2,
             mapping = aes(x = index, y = value, color = Design),
             inherit.aes = FALSE) +
  theme_bw() +
  ylim(0, 1) +
  labs(x = "Stage Index", y = element_blank())
plot(epph.plt)

if(beta_setting == 0){
  ggsave(
    filename = paste0(
      "lm_", "_scen", scenario,
      "_beta", beta_setting,
      "_N", Nttl,
      "_sigmasq", sigmasq,
      "_numSims", numSims,
      "_epphs", ".pdf"),
    plot = last_plot(),
    width = 6.5, height = 3.5, units = c("in")
  )
}
if(beta_setting == 4){
  ggsave(
    filename = paste0(
      "lm_", "_scen", scenario,
      "_beta", beta_setting,
      "_height", height, 
      "_discontinuity", discontinuity,
      "_N", Nttl,
      "_sigmasq", sigmasq,
      "_numSims", numSims,
      "_epphs", ".pdf"),
    plot = last_plot(),
    width = 6.5, height = 3.5, units = c("in")
  )
}

################################################################################

# sequential designs
PPH_seq = data.frame()
for(j in 1:numSims){
  # sequence of PPHs for each design
  PPH_seq.g = getPPHseq(grid_sims[[j]], models, Nttl, fT, sigmasq, 
                        initial.data = FALSE, randomize.order = TRUE) 
  PPH_seq.dl = getPPHseq(doptlin_sims[[j]], models, Nttl, fT, sigmasq, 
                         initial.data = FALSE, randomize.order = TRUE) 
  PPH_seq.dq = getPPHseq(doptquad_sims[[j]], models, Nttl, fT, sigmasq, 
                         initial.data = FALSE, randomize.order = TRUE) 
  PPH_seq.h = getPPHseq(hybrid_sims[[j]], models, Nttl, fT, sigmasq, 
                        initial.data = FALSE, randomize.order = TRUE) 
  PPH_seq.sm = getPPHseq(seqmed_sims[[j]], models, Nttl, fT, sigmasq) 
  # master data frame
  PPH_seq.g$Design = "Grid"
  PPH_seq.dl$Design = "DOptLin."
  PPH_seq.dq$Design = "DOptQuadr."
  PPH_seq.h$Design = "Hybrid"
  PPH_seq.sm$Design = "SeqMED"
  PPH_seq.tmp = rbind(PPH_seq.g, PPH_seq.dl, PPH_seq.dq, PPH_seq.h, PPH_seq.sm)
  PPH_seq.tmp$sim = j
  PPH_seq = rbind(PPH_seq, PPH_seq.tmp)
}

PPHmean_seq = aggregate(
  PPH_seq[, names(PPH_seq)[1:length(models)]], 
  by = list(PPH_seq[, "Design"], PPH_seq[, "index"]), 
  FUN = function(x) mean(x, na.rm = TRUE))
names(PPHmean_seq)[c(1, 2)] = c("Design", "index")

PPHmean_gg = PPHmean_seq
PPHmean_gg = melt(PPHmean_gg, id.vars = c("Design", "index"), 
                  measure.vars = paste0("H", 0:(length(models) - 1), sep = ""), 
                  variable.name = "hypothesis")
design_names = rev(c("SeqMED", "DOptLin.", "DOptQuadr.", "Grid", "Hybrid"))
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
PPHmean_gg = setorder(PPHmean_gg, cols = "Design")
PPHmean_gg2 = PPHmean_gg[PPHmean_gg$index == Nttl, ]
PPHmean_gg2$Design = factor(PPHmean_gg2$Design, levels = design_names)
PPHmean_gg2 = setorder(PPHmean_gg2, cols = "Design")
epph.plt2 = ggplot(PPHmean_gg, aes(x = index, y = value, color = Design,
                                  linetype = Design, shape = Design)) +
  facet_wrap(~hypothesis) +
  geom_path() +
  # scale_linetype_manual(values=c(rep("dashed", 4), rep("solid", 2))) +
  geom_point(data = PPHmean_gg2, 
             mapping = aes(x = index, y = value, color = Design),
             inherit.aes = FALSE) +
  theme_bw() +
  ylim(0, 1) + 
  labs(x = "Stage Index", y = element_blank())
plot(epph.plt2)

if(beta_setting == 0){
  ggsave(
    filename = paste0(
      "lm_", "_scen", scenario,
      "_beta", beta_setting,
      "_N", Nttl,
      "_sigmasq", sigmasq,
      "_numSims", numSims,
      "_epphs_seq", ".pdf"),
    plot = last_plot(),
    width = 6.5, height = 3.5, units = c("in")
  )
}
if(beta_setting == 4){
  ggsave(
    filename = paste0(
      "lm_", "_scen", scenario,
      "_beta", beta_setting,
      "_height", height, 
      "_discontinuity", discontinuity,
      "_N", Nttl,
      "_sigmasq", sigmasq,
      "_numSims", numSims,
      "_epphs_seq", ".pdf"),
    plot = last_plot(),
    width = 6.5, height = 3.5, units = c("in")
  )
}

################################################################################
# plot the estimated true curve?
################################################################################
sim.idx = 1
design.name = "seqmed"

if(design.name == "seqmed"){
  design_sim = seqmed_sims[[sim.idx]]
  x = c(design_sim$x.in, design_sim$x.new)
  y = c(design_sim$y.in, design_sim$y.new)
} else if(design.name == "hybrid"){
  design_sim = hybrid_sims[[sim.idx]]
  x = design_sim$x
  y = design_sim$y
} else if(design.name == "doptlinear"){
  design_sim = doptlin_sims[[sim.idx]]
  x = design_sim$x
  y = design_sim$y
} else if(design.name == "doptquadratic"){
  design_sim = doptquad_sims[[sim.idx]]
  x = design_sim$x
  y = design_sim$y
} else if(design.name == "grid"){
  design_sim = grid_sims[[sim.idx]]
  x = design_sim$x
  y = design_sim$y
}

posteriorT = getBetaPosterior(
  y, 
  modelT$designMat(x), 
  modelT$beta.mean, 
  modelT$beta.var, 
  sigmasq
)
posteriorT$mean
predictiveT = modelT$designMat(candidates) %*% posteriorT$mean

plot(candidates, predictiveT, type = "l")


