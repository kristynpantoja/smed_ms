################################################################################
# last updated: 2/25/22
# purpose: to create a list of seqmed simulations
# scenario 1:
#   linear vs. quadratic,
#   where the true function is quadratic
# scenario 2:
#   linear vs. quadratic,
#   where the true function is cubic
rm(list = ls())

scenario = 1 # 1, 2
given_Dinit = FALSE

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

# for box-hill design
source(paste(functions_dir, "/boxhill.R", sep = ""))

# for D-optimal design
library(AlgDesign)

# set up parallelization
rng.seed = 123 # 123, 345

library(mvtnorm)
library(ggplot2)
library(reshape2)
library(ggpubr)
library(tidyverse)
gg_color_hue = function(n) {
  hues = seq(15, 275, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

################################################################################
# simulation settings, shared for both scenarios (linear vs. quadratic)
################################################################################

# simulations settings
numSims = 100 # numSims = 500 & numSeq = 12 OR numSims = 100 & numSeq = 100
numSeq = 100 # 12, 100
seqN = 1
Nttl = numSeq * seqN
xmin = -1
xmax = 1
numCandidates = 10^3 + 1
candidates = seq(from = xmin, to = xmax, length.out = numCandidates)
if(scenario == 1){
  if(Nttl == 100){
    sigmasq = 0.28
  } else if(Nttl == 12){
    sigmasq = 0.04
  }
} else if(scenario == 2){
  if(Nttl == 100){
    sigmasq = 0.21
  } else if(Nttl == 12){
    sigmasq = 0.038
  }
}
alpha = 1

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
# curve(fT, from = xmin, to = xmax)

################################################################################
# import box & hill and seqmed simulations
################################################################################


if(given_Dinit){
  N0 = 3
  Dinit_label = paste0("_Dinit", N0)
}  else{
  Dinit_label = ""
}
seqmed_sims = readRDS(paste0(
  output_dir,
  "/scenario", scenario, 
  "_seqmed",
  Dinit_label, 
  "_N", Nttl, 
  "_sigmasq", sigmasq,
  "_alpha", alpha,
  "_numSims", numSims,
  "_seed", rng.seed,
  ".rds"
))
boxhill_sims = readRDS(paste0(
  output_dir, "/scenario", scenario, 
  "_boxhill", 
  "_N", Nttl, 
  "_sigmasq", sigmasq,
  "_numSims", numSims,
  "_seed", rng.seed,
  ".rds"))

################################################################################
# non-sequential designs
################################################################################

# use AlgDesign to obtain Doptimal designs #

# consider linear model, y = b0 + b1 x
candidates_named = data.frame(x = candidates)
res_Fed_Doptlin = optFederov(
  ~1+I(x), data = candidates_named, approximate = TRUE, criterion = "D")
res_Fed_Doptlin$design

# consider quadratic model, y = b0 + b1 x + b2 x^2
res_Fed_Doptquad = optFederov(
  ~1+x+I(x^2), data = candidates_named, approximate = TRUE, criterion = "D")
res_Fed_Doptquad$design

# define designs #

# space-filling (grid)
space_filling = seq(from = xmin, to = xmax, length.out = Nttl)

grid_sims = list()
doptlin_sims = list()
doptquad_sims = list()
hybrid_sims = list()
set.seed(rng.seed)
for(j in 1:numSims){
  # Doptimal - linear
  num_DoptL = nrow(res_Fed_Doptlin$design) # number of support points
  pts_DoptL = res_Fed_Doptlin$design[, "x"] # the actual support points
  numeachL = floor(Nttl / num_DoptL) # smallest number of times can recycle spt pts
  dopt_linear0 = c(sapply(pts_DoptL, function(x) rep(x, numeachL)))
  dopt_linear = sample(c(
    dopt_linear0, 
    sample(pts_DoptL, size = Nttl - numeachL * num_DoptL, replace = FALSE)
  ), size = Nttl, replace = FALSE)
  # # check:
  # res_Fed_Doptlin$design
  # table(dopt_linear) / Nttl
  
  # Doptimal - quadratic
  num_DoptQ = nrow(res_Fed_Doptquad$design) # number of support points
  pts_DoptQ = res_Fed_Doptquad$design[, "x"] # the actual support points
  numeachQ = floor(Nttl / num_DoptQ) # smallest number of times can recycle spt pts
  dopt_quadratic0 = c(sapply(pts_DoptQ, function(x) rep(x, numeachQ)))
  dopt_quadratic = sample(c(
    dopt_quadratic0, 
    sample(pts_DoptQ, size = Nttl - numeachQ * num_DoptQ, replace = FALSE)
  ), size = Nttl, replace = FALSE)
  # # check:
  # res_Fed_Doptquad$design
  # table(dopt_quadratic) / Nttl
  
  # # half space-filling, half quadratic Doptimal, 
  # #   assumes Nttl is divisible by 2
  # supportpt_assgnmt_hybrid = cut(
  #   sample(1:(Nttl / 2), size = Nttl / 2, replace = FALSE), # shuffle
  #   breaks = num_supportpts_Doptquad, labels = FALSE)
  # hybrid_grid_doptq = rep(NA, Nttl / 2)
  # for(i in 1:num_supportpts_Doptquad){
  #   hybrid_grid_doptq[supportpt_assgnmt_hybrid == i] = 
  #     supportpts_Doptquad[i]
  # }
  # hybrid_grid_doptq[(Nttl / 2 + 1):Nttl] =c(
  #   seq(
  #     from = supportpts_Doptquad[1], to = supportpts_Doptquad[2],
  #     length.out = (Nttl / 4) + 2)[-c(1, (Nttl / 4) + 2)],
  #   seq(
  #     from = supportpts_Doptquad[2], to = supportpts_Doptquad[3],
  #     length.out = (Nttl / 4) + 2)[-c(1, (Nttl / 4) + 2)]
  # )
  # hybrid_grid_doptq = rev(hybrid_grid_doptq) # spacefilling -> doptimal
  
  # simulations #
  
  space_filling.tmp = sample(space_filling, replace = FALSE)
  grid_sims[[j]] = list(
    x = space_filling.tmp,
    y = sapply(space_filling.tmp, FUN = function(x) simulateY_fromfunction(
      x = x, true.function = fT, error.var = sigmasq)))
  dopt_linear.tmp = sample(dopt_linear, replace = FALSE)
  doptlin_sims[[j]] = list(
    x = dopt_linear.tmp,
    y = sapply(dopt_linear.tmp, FUN = function(x) simulateY_fromfunction(
      x = x, true.function = fT, error.var = sigmasq)))
  dopt_quadratic.tmp = sample(dopt_quadratic, replace = FALSE)
  doptquad_sims[[j]] = list(
    x = dopt_quadratic.tmp,
    y = sapply(dopt_quadratic.tmp, FUN = function(x) simulateY_fromfunction(
      x = x, true.function = fT, error.var = sigmasq)))
  # hybrid_grid_doptq.tmp = sample(hybrid_grid_doptq, replace = FALSE)
  # hybrid_sims[[j]] = list(
  #   x = hybrid_grid_doptq.tmp,
  #   y = sapply(hybrid_grid_doptq.tmp, FUN = function(x) simulateY_fromfunction(
  #     x = x, true.function = fT, error.var = sigmasq)))
}


################################################################################
################################################################################
################################################################################
# plots!!!
################################################################################
################################################################################
################################################################################


################################################################################
# plot the posterior probabilities of the hypotheses
################################################################################
include_hybrid = FALSE

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

# ################################################################################
# # sequential & non-sequential designs' epph for the given alpha value
# 
# # non-sequential designs
# PPH_df = data.frame()
# for(j in 1:numSims){
#   # sequence of PPHs for each design
#   PPH_grid = getPPH(
#     grid_sims[[j]], models, fT, sigmasq)
#   PPH_doptl = getPPH(
#     doptlin_sims[[j]], models, fT, sigmasq)
#   PPH_doptq = getPPH(
#     doptquad_sims[[j]], models, fT, sigmasq)
#   # PPH_hybrid = getPPH(
#   #   hybrid_sims[[j]], models, fT, sigmasq)
#   # master data frame
#   PPH_grid$Design = "Grid"
#   PPH_doptl$Design = "DOptLin."
#   PPH_doptq$Design = "DOptQuadr."
#   # PPH_hybrid$Design = "Hybrid"
#   PPH.tmp = rbind(PPH_grid, PPH_doptl, PPH_doptq) #, PPH_hybrid)
#   PPH.tmp$sim = j
#   PPH_df = rbind(PPH_df, PPH.tmp)
# }
# PPHmean = aggregate(
#   PPH_df[, names(PPH_df)[1:length(models)]],
#   by = list(PPH_df[, "Design"]), FUN = function(x) mean(x, na.rm = TRUE))
# names(PPHmean)[1] = "Design"
# PPHmean$index = Nttl
# # but we want a line, so allow interpolation by setting PPHmean$index = 0 too
# PPHmean2 = PPHmean
# PPHmean2$index = 0
# PPHmean = rbind(PPHmean, PPHmean2)
# 
# # sequential designs
# PPH_seq = data.frame()
# for(j in 1:numSims){
#   # sequence of PPHs for each design
#   PPH_seq.sm = getPPHseq(seqmed_sims[[j]], models, Nttl, fT, sigmasq)
#   PPH_seq.bh = getPPHseq(boxhill_sims[[j]], models, Nttl, fT, sigmasq)
#   # master data frame
#   PPH_seq.sm$Design = "SeqMED"
#   PPH_seq.bh$Design = "BoxHill"
#   PPH_seq.tmp = rbind(PPH_seq.sm, PPH_seq.bh)
#   PPH_seq.tmp$sim = j
#   PPH_seq = rbind(PPH_seq, PPH_seq.tmp)
# }
# 
# PPHmean_seq = aggregate(
#   PPH_seq[, names(PPH_seq)[1:length(models)]],
#   by = list(PPH_seq[, "Design"], PPH_seq[, "index"]),
#   FUN = function(x) mean(x, na.rm = TRUE))
# names(PPHmean_seq)[c(1, 2)] = c("Design", "index")
# 
# PPHmean_gg = rbind(PPHmean, PPHmean_seq)
# PPHmean_gg = reshape2::melt(
#   PPHmean_gg, id.vars = c("Design", "index"),
#   measure.vars = paste0("H", 0:(length(models) - 1), sep = ""),
#   variable.name = "hypothesis")
# design_names = rev(c(
#   "SeqMED", "BoxHill", "DOptLin.", "DOptQuadr.", "Grid")) #, "Hybrid"))
# PPHmean_gg$Design = factor(PPHmean_gg$Design, levels = design_names)
# if(scenario == 1){
#   PPHmean_gg$hypothesis = factor(
#     PPHmean_gg$hypothesis,
#     levels = paste0("H", 0:(length(models) - 1), sep = ""),
#     labels = paste0("Case ", scenario, ", H", 0:(length(models) - 1), sep = ""))
# } else if(scenario == 2){
#   PPHmean_gg$hypothesis = factor(
#     PPHmean_gg$hypothesis,
#     levels = paste0("H", 0:(length(models) - 1), sep = ""),
#     labels = paste0("Case ", scenario, ", H", c(0, 1, "T"), sep = ""))
# }
# PPHmean_gg = setorder(PPHmean_gg, cols = "Design")
# PPHmean_gg2 = PPHmean_gg[PPHmean_gg$index == Nttl, ]
# PPHmean_gg2$Design = factor(PPHmean_gg2$Design, levels = design_names)
# PPHmean_gg2 = setorder(PPHmean_gg2, cols = "Design")
# # if(include_hybrid){
# #   num_fixed_designs = 4
# # } else{
#   num_fixed_designs = 3
# # }
# # if(include_hybrid){
# #   epph.plt = ggplot(PPHmean_gg, aes(x = index, y = value, color = Design,
# #                                     linetype = Design, shape = Design)) +
# #     facet_wrap(~hypothesis) +
# #     geom_path() +
# #     scale_linetype_manual(values=c(
# #       rep("dashed", num_fixed_designs), rep("solid", 2))) +
# #     geom_point(data = PPHmean_gg2,
# #                mapping = aes(x = index, y = value, color = Design),
# #                inherit.aes = FALSE) +
# #     theme_bw() +
# #     ylim(0, 1) +
# #     labs(x = "Stage Index", y = element_blank())
# #   epph.plt
# # } else{
#   PPHmean_gg_nohybrid = PPHmean_gg %>% filter(Design != "Hybrid")
#   PPHmean_gg2_nohybrid = PPHmean_gg2 %>% filter(Design != "Hybrid")
#   epph.plt = ggplot(
#     PPHmean_gg_nohybrid, aes(x = index, y = value, color = Design,
#                              linetype = Design, shape = Design)) +
#     facet_wrap(~hypothesis) +
#     geom_path() +
#     scale_linetype_manual(
#       values=c(rep("dashed", num_fixed_designs), rep("solid", 2))) +
#     geom_point(
#       data = PPHmean_gg2_nohybrid,
#       mapping = aes(x = index, y = value, color = Design),
#       inherit.aes = FALSE) +
#     theme_bw() +
#     ylim(0, 1) +
#     labs(x = "Stage Index", y = element_blank())
#   epph.plt
# # }

# if(scenario == 1){
#   epph_scen1 = epph.plt
# } else if(scenario == 2){
#   epph_scen2 = epph.plt + theme(legend.position = "none")
#   if(Nttl == 12){
#     epph_scen2 = epph_scen2 + scale_x_continuous(breaks = seq(0, 12, by = 3))
#   }
# }
# if(!is.null(epph_scen1) && !is.null(epph_scen2)){
#   ggarrange(epph_scen1, epph_scen2, nrow = 2, ncol = 1, widths = 0.9)
# 
#   # manuscript plot
#   ggsave(
#     filename = paste0(
#       "lm", Dinit_label, "_hybrid", include_hybrid,
#       "_epphs", ".pdf"),
#     plot = last_plot(),
#     width = 6.5, height = 5.5, units = c("in")
#   )
# }



################################################################################
# all the alphas
################################################################################

# alphas = c(0, 1, 10, 25, 50, 100)
# alphas = c(0, 1, 10, 100)
alphas = c(0, 1, 10)

seqmed_sims_alphas = list()
min_alpha_used = matrix(NA, length(alphas), numSims)
for(i in 1:length(alphas)){
  seqmed_sims_alphas[[i]] = readRDS(paste0(
    output_dir,
    "/scenario", scenario, 
    "_seqmed",
    Dinit_label, 
    "_N", Nttl, 
    "_sigmasq", sigmasq,
    "_alpha", alphas[i],
    "_numSims", numSims,
    "_seed", rng.seed,
    ".rds"
  ))
  min_alpha_used[i, ] = sapply(
    seqmed_sims_alphas[[i]], function(sim) min(sim$alpha_seq))
}
avg_min_alpha_used = apply(min_alpha_used, 1, mean)


################################################################################
# sequential & non-sequential designs' epph

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
  # PPH_hybrid = getPPH(
  #   hybrid_sims[[j]], models, fT, sigmasq)
  # master data frame
  PPH_grid$Design = paste0("(vii) ", "Grid")
  PPH_doptl$Design = paste0("(v) ", "DOptLin.")
  PPH_doptq$Design = paste0("(vi) ", "DOptQuadr.")
  # PPH_hybrid$Design = "Hybrid"
  PPH.tmp = rbind(PPH_grid, PPH_doptl, PPH_doptq) #, PPH_hybrid)
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
  PPH_seq.tmp = data.frame()
  for(m in 1:length(alphas)){
    # sequence of PPHs for each design
    PPH_seq.sm.m = getPPHseq(
      seqmed_sims_alphas[[m]][[j]], models, Nttl, fT, sigmasq)
    # master data frame
    PPH_seq.sm.m$Design = paste0("(", paste(rep("i", m), collapse = ""), ")", " SeqMED ", alphas[m])
    PPH_seq.tmp = rbind(PPH_seq.tmp, PPH_seq.sm.m)
  }
  PPH_seq.bh = getPPHseq(boxhill_sims[[j]], models, Nttl, fT, sigmasq)
  PPH_seq.bh$Design = paste0("(iv) ", "BoxHill")
  PPH_seq.tmp = rbind(PPH_seq.tmp, PPH_seq.bh)
  PPH_seq.tmp$sim = j
  PPH_seq = rbind(PPH_seq, PPH_seq.tmp)
}

PPHmean_seq = aggregate(
  PPH_seq[, names(PPH_seq)[1:length(models)]],
  by = list(PPH_seq[, "Design"], PPH_seq[, "index"]),
  FUN = function(x) mean(x, na.rm = TRUE))
names(PPHmean_seq)[c(1, 2)] = c("Design", "index")

PPHmean_gg = rbind(PPHmean, PPHmean_seq)
PPHmean_gg = reshape2::melt(
  PPHmean_gg, id.vars = c("Design", "index"),
  measure.vars = paste0("H", 0:(length(models) - 1), sep = ""),
  variable.name = "hypothesis")
# design_names = c(
#   "Grid", "DOptLin.", "DOptQuadr.", "BoxHill",
#   paste("SeqMED", alphas, sep = " "))
# PPHmean_gg$Design = factor(PPHmean_gg$Design, levels = design_names)
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
# PPHmean_gg = setorder(PPHmean_gg, cols = "Design")
PPHmean_gg2 = PPHmean_gg[PPHmean_gg$index == Nttl, ]
# PPHmean_gg2$Design = factor(PPHmean_gg2$Design, levels = design_names)
# PPHmean_gg2 = setorder(PPHmean_gg2, cols = "Design")
# if(include_hybrid){
#   epph.plt3 = ggplot(PPHmean_gg, aes(x = index, y = value, color = Design,
#                                      linetype = Design)) +
#     facet_wrap(~hypothesis) +
#     geom_path() +
#     scale_linetype_manual(
#       values=c(
#         rep("dashed", num_fixed_designs), rep("solid", length(alphas) + 1))) +
#     geom_point(data = PPHmean_gg2,
#                mapping = aes(x = index, y = value, color = Design),
#                inherit.aes = FALSE) +
#     theme_bw() +
#     ylim(0, 1) +
#     labs(x = "Stage Index", y = element_blank())
#   if(Nttl == 12){
#     epph.plt3 = epph.plt3 + scale_x_continuous(breaks = seq(0, 12, by = 3))
#   }
#   epph.plt3
# } else{
  # PPHmean_gg_nohybrid = PPHmean_gg %>% filter(Design != "Hybrid")
  # PPHmean_gg2_nohybrid = PPHmean_gg2 %>% filter(Design != "Hybrid")
num_fixed_designs = 3
  epph.plt3 = ggplot(
    PPHmean_gg, aes(x = index, y = value, color = Design,
                             linetype = Design)) +
    facet_wrap(~hypothesis) +
    geom_path() +
    scale_linetype_manual(values=c(
      rep("solid", length(alphas) + 1), rep("dashed", num_fixed_designs))) +
    geom_point(
      data = PPHmean_gg2,
      mapping = aes(x = index, y = value, color = Design),
      inherit.aes = FALSE) +
    theme_bw() +
    ylim(0, 1) +
    labs(x = "Index", y = element_blank())
  if(Nttl == 12){
    epph.plt3 = epph.plt3 + scale_x_continuous(breaks = seq(0, 12, by = 3))
  }
  epph.plt3
# }


# # manuscript plot
# ggsave(
#   filename = paste0(
#     "lm", "_scen", scenario, Dinit_label, "_hybrid", include_hybrid,
#     "_epph_alphas_", paste(alphas, collapse = "_"), ".pdf"),
#   plot = last_plot(),
#   width = 6.5, height = 3.5, units = c("in")
# )

if(scenario == 1){
  epph_scen1 = epph.plt3
} else if(scenario == 2){
  epph_scen2 = epph.plt3 + theme(legend.position = "none")
}
if(!is.null(epph_scen1) && !is.null(epph_scen2) && numSeq == 100){
  ggarrange(epph_scen1, epph_scen2, nrow = 2, ncol = 1)

  # manuscript plot
  ggsave(
    filename = paste0(
      "lm", Dinit_label, "_hybrid", include_hybrid, 
      "_epphs_alphas_", paste(alphas, collapse = "_"), ".pdf"),
    plot = last_plot(),
    width = 6.5, height = 5.5, units = c("in")
  )
}

# # all sequential epph plot #####################################################
# 
# PPH_seq_alpha = data.frame()
# for(j in 1:numSims){
#   # sequence of PPHs for each design
#   PPH_seq.g = getPPHseq(grid_sims[[j]], models, Nttl, fT, sigmasq, 
#                         initial.data = FALSE, randomize.order = TRUE) 
#   PPH_seq.dl = getPPHseq(doptlin_sims[[j]], models, Nttl, fT, sigmasq, 
#                          initial.data = FALSE, randomize.order = TRUE) 
#   PPH_seq.dq = getPPHseq(doptquad_sims[[j]], models, Nttl, fT, sigmasq, 
#                          initial.data = FALSE, randomize.order = TRUE) 
#   PPH_seq.h = getPPHseq(hybrid_sims[[j]], models, Nttl, fT, sigmasq, 
#                         initial.data = FALSE, randomize.halves.order = TRUE)
#     PPH_seq.sm0 = getPPHseq(
#       seqmed_sims_alphas[[1]][[j]], models, Nttl, fT, sigmasq) 
#     PPH_seq.sm1 = getPPHseq(
#       seqmed_sims_alphas[[2]][[j]], models, Nttl, fT, sigmasq) 
#     PPH_seq.sm5 = getPPHseq(
#       seqmed_sims_alphas[[3]][[j]], models, Nttl, fT, sigmasq) 
#     PPH_seq.sm10 = getPPHseq(
#       seqmed_sims_alphas[[4]][[j]], models, Nttl, fT, sigmasq) 
#     PPH_seq.bh = getPPHseq(
#       boxhill_sims[[j]], models, Nttl, fT, sigmasq) 
#   # master data frame
#   PPH_seq.g$Design = "Grid"
#   PPH_seq.dl$Design = "DOptLin."
#   PPH_seq.dq$Design = "DOptQuadr."
#   PPH_seq.h$Design = "Hybrid"
#   PPH_seq.sm0$Design = "SeqMED 0"
#   PPH_seq.sm1$Design = "SeqMED 1"
#   PPH_seq.sm5$Design = "SeqMED 5"
#   PPH_seq.sm10$Design = "SeqMED 10"
#   PPH_seq.bh$Design = "BoxHill"
#   PPH_seq.tmp = rbind(
#     PPH_seq.g, PPH_seq.dl, PPH_seq.dq, PPH_seq.h, 
#     PPH_seq.sm0, PPH_seq.sm1, PPH_seq.sm5, PPH_seq.sm10, PPH_seq.bh)
#   PPH_seq.tmp$sim = j
#   PPH_seq_alpha = rbind(PPH_seq_alpha, PPH_seq.tmp)
# }
# 
# PPHmean_seq = aggregate(
#   PPH_seq_alpha[, names(PPH_seq_alpha)[1:length(models)]], 
#   by = list(PPH_seq_alpha[, "Design"], PPH_seq_alpha[, "index"]), 
#   FUN = function(x) mean(x, na.rm = TRUE))
# names(PPHmean_seq)[c(1, 2)] = c("Design", "index")
# 
# PPHmean_gg = PPHmean_seq
# PPHmean_gg = reshape2::melt(
#   PPHmean_gg, id.vars = c("Design", "index"), 
#   measure.vars = paste0("H", 0:(length(models) - 1), sep = ""), 
#   variable.name = "hypothesis")
# design_names = rev(c(
#   "SeqMED 10", "SeqMED 5", "SeqMED 1", "SeqMED 0", 
#   "BoxHill",
#   "DOptLin.", "DOptQuadr.", "Grid", "Hybrid"))
# PPHmean_gg$Design = factor(PPHmean_gg$Design, levels = design_names)
# if(scenario == 1){
#   PPHmean_gg$hypothesis = factor(
#     PPHmean_gg$hypothesis, 
#     levels = paste0("H", 0:(length(models) - 1), sep = ""), 
#     labels = paste0("Case ", scenario, ", H", 0:(length(models) - 1), sep = ""))
# } else if(scenario == 2){
#   PPHmean_gg$hypothesis = factor(
#     PPHmean_gg$hypothesis, 
#     levels = paste0("H", 0:(length(models) - 1), sep = ""), 
#     labels = paste0("Case ", scenario, ", H", c(0, 1, "T"), sep = ""))
# }
# PPHmean_gg = setorder(PPHmean_gg, cols = "Design")
# PPHmean_gg2 = PPHmean_gg[PPHmean_gg$index == Nttl, ]
# PPHmean_gg2$Design = factor(PPHmean_gg2$Design, levels = design_names)
# PPHmean_gg2 = setorder(PPHmean_gg2, cols = "Design")
# epph.plt4 = ggplot(PPHmean_gg, aes(x = index, y = value, color = Design,
#                                    linetype = Design)) +
#   facet_wrap(~hypothesis) +
#   geom_path() +
#   scale_linetype_manual(values=c(rep("dashed", 3), rep("solid", 5))) +
#   geom_point(data = PPHmean_gg2, 
#              mapping = aes(x = index, y = value, color = Design), 
#              inherit.aes = FALSE) +
#   theme_bw() +
#   ylim(0, 1) + 
#   labs(x = "Stage Index", y = element_blank())
# epph.plt4



