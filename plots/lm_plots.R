
scenario = 1
# scenario = 1 : linear vs quadratic, f is quadratic
# scenario = 2 : linear vs quadratic, f is cubic
scenario1_illustration = TRUE
# scenario1_illustration = FALSE : for Figure 1
# scenario1_illustration = TRUE : for other figures in linear models examples
text_size = 12

################################################################################
# Sources/Libraries
################################################################################
output_dir = "plots"
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
library(foreach)
library(future)
library(doFuture)
library(parallel)
registerDoFuture()
nworkers = detectCores() / 2
plan(multisession, workers = nworkers)
library(rngtools)
library(doRNG)
rng.seed = 123 # 123, 345
registerDoRNG(rng.seed)

# other libraries
library(expm)
library(matrixStats)
library(MASS)
library(mvtnorm)
library(ggplot2)
library(reshape2)
library(tidyverse)
gg_color_hue = function(n) {
  hues = seq(15, 275, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

################################################################################
# simulation settings, shared for both scenarios (linear vs. quadratic)
################################################################################

# simulations settings
numSims = 100
numSeq = 100
seqN = 1
Nttl = numSeq * seqN
xmin = -1
xmax = 1
numCandidates = 10^3 + 1
candidates = seq(from = xmin, to = xmax, length.out = numCandidates)
if(scenario == 1){
  if(Nttl == 100){
    if(scenario1_illustration){
      sigmasq = 0.1
    } else{
      sigmasq = 0.28
    }
  } 
} else if(scenario == 2){
  if(Nttl == 100){
    sigmasq = 0.21
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
# run simulations
################################################################################

# generate seqmeds
alphas = c(0, 1, 10)
for(l in 1:length(alphas)){
  registerDoRNG(rng.seed)
  seqmed_list = foreach(i = 1:numSims) %dopar% {
    print(paste0("starting simulation ", i, " out of ", numSims))
    SeqMED(
      y.in = NULL, x.in = NULL, true.function = fT,
      model0 = model0, model1 = model1, 
      error.var = sigmasq, xmin = xmin, xmax = xmax,
      candidates = candidates, numSeq = numSeq, seqN = seqN, 
      alpha_seq = alphas[l])
  }
  saveRDS(seqmed_list, paste0(
    output_dir, "/scenario", scenario, 
    "_seqmed", 
    "_N", Nttl, 
    "_sigmasq", sigmasq,
    "_alpha", alphas[l],
    "_numSims", numSims,
    "_seed", rng.seed,
    ".rds"))
}

# generate boxhills
registerDoRNG(rng.seed)
bh_list = foreach(i = 1:numSims) %dorng% {
  print(paste0("starting simulation ", i, " out of ", numSims))
  BH_m2(NULL, NULL, prior_probs, model0, model1, Nttl, 
        candidates, fT, sigmasq)
}
saveRDS(bh_list, paste0(
  output_dir, "/scenario", scenario, 
  "_boxhill", 
  "_N", Nttl, 
  "_sigmasq", sigmasq,
  "_numSims", numSims,
  "_seed", rng.seed,
  ".rds"))

# import them
if(!scenario1_illustration){
  seqmed_sims_alphas = list()
  for(i in 1:length(alphas)){
    seqmed_sims_alphas[[i]] = readRDS(paste0(
      output_dir,
      "/scenario", scenario, 
      "_seqmed",
      "_N", Nttl, 
      "_sigmasq", sigmasq,
      "_alpha", alphas[i],
      "_numSims", numSims,
      "_seed", rng.seed,
      ".rds"
    ))
  }
  boxhill_sims = readRDS(paste0(
    output_dir, "/scenario", scenario, 
    "_boxhill", 
    "_N", Nttl, 
    "_sigmasq", sigmasq,
    "_numSims", numSims,
    "_seed", rng.seed,
    ".rds"))
}
seqmed_sims = readRDS(paste0(
  output_dir,
  "/scenario", scenario, 
  "_seqmed",
  "_N", Nttl, 
  "_sigmasq", sigmasq,
  "_alpha", 1,
  "_numSims", numSims,
  "_seed", rng.seed,
  ".rds"
))

################################################################################
# non-sequential designs
################################################################################

if(!scenario1_illustration){
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
  
  # consider D-optimal design for b2
  res_Fed_Ds = list()
  res_Fed_Ds$design = res_Fed_Doptquad$design
  res_Fed_Ds$design[, "Proportion"] = c(0.25, 0.5, 0.25)
  
  # define designs #
  
  # space-filling (grid)
  space_filling = seq(from = xmin, to = xmax, length.out = Nttl)
  
  grid_sims = list()
  doptlin_sims = list()
  doptquad_sims = list()
  doptequalinterest_sims = list()
  set.seed(rng.seed)
  for(j in 1:numSims){
    # Doptimal - linear
    num_DoptL = nrow(res_Fed_Doptlin$design) # number of support points
    pts_DoptL = res_Fed_Doptlin$design[, "x"] # the actual support points
    wts_DoptL = res_Fed_Doptlin$design[, "Proportion"] # their weights
    n_DoptL = floor(Nttl * wts_DoptL) # n for each support point
    dopt_linear0 = c()
    for(i in 1:num_DoptL){
      dopt_linear0 = c(dopt_linear0, rep(pts_DoptL[i], n_DoptL[i]))
    }
    dopt_linear = sample(c( # shuffle
      dopt_linear0, 
      sample(pts_DoptL, size = Nttl - sum(n_DoptL), replace = FALSE)
    ), size = Nttl, replace = FALSE)
    
    # Doptimal - quadratic
    num_DoptQ = nrow(res_Fed_Doptquad$design) # number of support points
    pts_DoptQ = res_Fed_Doptquad$design[, "x"] # the actual support points
    wts_DoptQ = res_Fed_Doptquad$design[, "Proportion"] # their weights
    n_DoptQ = floor(Nttl * wts_DoptQ) # n for each support point
    dopt_quadratic0 = c()
    for(i in 1:num_DoptQ){
      dopt_quadratic0 = c(dopt_quadratic0, rep(pts_DoptQ[i], n_DoptQ[i]))
    }
    dopt_quadratic = sample(c( # shuffle
      dopt_quadratic0, 
      sample(pts_DoptQ, size = Nttl - sum(n_DoptQ), replace = FALSE)
    ), size = Nttl, replace = FALSE)
    
    # Ds optimal - linear nested in quadratic
    num_Ds = nrow(res_Fed_Ds$design) # number of support points
    pts_Ds = res_Fed_Ds$design[, "x"] # the actual support points
    wts_Ds = res_Fed_Ds$design[, "Proportion"] # their weights
    n_Ds = floor(Nttl * wts_Ds) # n for each support point
    dopt_equalinterest0 = c()
    for(i in 1:num_Ds){
      dopt_equalinterest0 = c(dopt_equalinterest0, rep(pts_Ds[i], n_Ds[i]))
    }
    dopt_equalinterest = sample(c( # shuffle
      dopt_equalinterest0, 
      sample(pts_Ds, size = Nttl - sum(n_Ds), replace = FALSE)
    ), size = Nttl, replace = FALSE)
    
    
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
    dopt_equalinterest.tmp = sample(dopt_equalinterest, replace = FALSE)
    doptequalinterest_sims[[j]] = list(
      x = dopt_equalinterest.tmp,
      y = sapply(dopt_equalinterest.tmp, FUN = function(x) simulateY_fromfunction(
        x = x, true.function = fT, error.var = sigmasq)))
  }
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

# scenario-specific models
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

# useful functions #############################################################

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

# design model probs ###########################################################

if(!scenario1_illustration){
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
    PPH_doptei = getPPH(
      doptequalinterest_sims[[j]], models, fT, sigmasq)
    # master data frame
    PPH_grid$Design = paste0("(vii) ", "Grid")
    PPH_doptl$Design = paste0("(v) ", "DOptLin.")
    PPH_doptq$Design = paste0("(vi) ", "DOptQuadr.")
    PPH_doptei$Design = paste0("(viii) ", "DOptEqInt.")
    PPH.tmp = rbind(PPH_grid, PPH_doptl, PPH_doptq, PPH_doptei)
    PPH.tmp$sim = j
    PPH_df = rbind(PPH_df, PPH.tmp)
  }
  PPHmean = aggregate(
    PPH_df[, names(PPH_df)[1:length(models)]],
    by = list(PPH_df[, "Design"]), FUN = function(x) mean(x, na.rm = TRUE))
  names(PPHmean)[1] = "Design"
  # add points
  PPHmean = PPHmean[rep(
    seq_len(nrow(PPHmean)), each = length(c(1, 25, 50, 75, 100))), ]
  PPHmean$index = rep(c(1, 25, 50, 75, 100), 4)
  
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
  PPHmean_gg2 = PPHmean_gg[PPHmean_gg$index %in% c(1, 25, 50, 75, 100), ]
  num_fixed_designs = 3
  PPHmean_gg = PPHmean_gg %>% 
    filter(Design != paste0("(viii) ", "DOptEqInt."))
  PPHmean_gg2 = PPHmean_gg2 %>% 
    filter(Design != paste0("(viii) ", "DOptEqInt."))
  if(scenario == 1){
    PPHmean_gg = PPHmean_gg %>% 
      filter(hypothesis != "Case 1, H0")
    PPHmean_gg2 = PPHmean_gg2 %>% 
      filter(hypothesis != "Case 1, H0")
  }
  epph.plt3 = ggplot(
    PPHmean_gg, aes(x = index, y = value, color = Design,
                    linetype = Design)) +
    facet_wrap(~hypothesis) +
    geom_path() +
    scale_linetype_manual(values=c(
      rep("solid", length(alphas) + 1), rep("dashed", num_fixed_designs))) +
    geom_point(
      data = PPHmean_gg2,
      mapping = aes(x = index, y = value, color = Design, shape = Design),
      size = 2, inherit.aes = FALSE) +
    ylim(0, 1) +
    labs(x = element_blank(), y = element_blank()) + 
    theme_bw() + 
    theme(text = element_text(size = text_size), 
          panel.background = element_rect(fill = "transparent"), 
          plot.background = element_rect(fill = "transparent", color = NA), 
          legend.background = element_rect(fill = "transparent"), 
          legend.box.background = element_rect(
            fill = "transparent", color = "transparent"))
  
  epph.plt3
  if(scenario == 1){
    ggsave(
      filename = paste0(
        "scen1_lm_epphs_alphas_", paste(alphas, collapse = "_"), ".pdf"),
      plot = epph.plt3,
      width = 4, height = 2.5, units = c("in")
    )
  } else if(scenario == 2){
    ggsave(
      filename = paste0(
        "scen2_lm_epphs_alphas_", paste(alphas, collapse = "_"), ".pdf"),
      plot = epph.plt3 + theme(legend.position = "none"),
      width = 5.5, height = 2.5, units = c("in")
    )
  }
}



################################################################################
# plot the designs
################################################################################

sim.idx = 5

# plot a seqmed
sm = seqmed_sims[[sim.idx]]
ggdata0 = data.frame(x = c(sm$x.in, sm$x.new), y = c(sm$y.in, sm$y.new))
plt0 = ggplot(ggdata0) +
  geom_histogram(binwidth = 0.12, closed = "right",
                 aes(x = x, y = after_stat(density))) +
  scale_x_continuous(breaks = c(-1, 0, 1)) +
  labs(x = element_blank(), y = element_blank()) + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),  
        text = element_text(size = text_size), 
        panel.background = element_rect(fill = "transparent"), 
        plot.background = element_rect(fill = "transparent", color = NA), 
        legend.background = element_rect(fill = "transparent"), 
        legend.box.background = element_rect(
          fill = "transparent", color = "transparent"))
plt1 = ggplot(ggdata0) +
  geom_point(aes(x, y), col = gg_color_hue(2)[1]) +
  stat_function(fun = fT) +
  scale_x_continuous(breaks = c(-1, 0, 1)) + 
  labs(x = element_blank(), y = element_blank()) + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),  
        text = element_text(size = text_size), 
        panel.background = element_rect(fill = "transparent"), 
        plot.background = element_rect(fill = "transparent", color = NA), 
        legend.background = element_rect(fill = "transparent"), 
        legend.box.background = element_rect(
          fill = "transparent", color = "transparent"))

if(!scenario1_illustration){
  # plot all of them #############################################################
  sim.idx2 = 1
  sm_x = c(seqmed_sims[[sim.idx2]]$x.in, seqmed_sims[[sim.idx2]]$x.new)
  bh_x = c(boxhill_sims[[sim.idx2]]$x.in, boxhill_sims[[sim.idx2]]$x.new)
  grid_x = grid_sims[[sim.idx2]]$x
  dl_x = doptlin_sims[[sim.idx2]]$x
  dq_x = doptquad_sims[[sim.idx2]]$x
  ggdata3 = rbind(
    data.frame(x = sm_x, Design = "SeqMED"), 
    data.frame(x = bh_x, Design = "BoxHill"), 
    data.frame(x = grid_x, Design = "Grid"), 
    data.frame(x = dl_x, Design = "DOptLin"), 
    data.frame(x = dq_x, Design = "DOptQuad")
  )
  ggdata3$Design = factor(
    ggdata3$Design, 
    levels = c("SeqMED", "Grid", "BoxHill", "DOptLin", "DOptQuad"))
  plt4 = ggplot(ggdata3, aes(x = x)) +
    facet_wrap(~Design, nrow = 1) +
    geom_histogram(binwidth = 0.12, 
                   aes(x = x)) +
    scale_x_continuous(breaks = c(-1, 0, 1)) +
    theme_bw() + 
    theme(
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
      panel.background = element_rect(fill = "transparent"), 
      plot.background = element_rect(fill = "transparent", color = NA), 
      legend.background = element_rect(fill = "transparent"), 
      legend.box.background = element_rect(
        fill = "transparent", color = "transparent"))
  
  plt4
  ggsave(
    filename = paste0("scen", scenario, "_seqmedex_designs.pdf"),
    plot = plt4,
    width = 6, height = 2.5, units = c("in"), bg = "transparent"
  )
}

# plot the wasserstein distance ################################################

library(data.table)

if(scenario == 1 & scenario1_illustration){
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
    scale_x_continuous(breaks = c(-1, 0, 1)) +
    scale_color_manual(
      values = c(gg_color_hue(4)[c(3, 4)], "black", gg_color_hue(4)[2])) +
    geom_path() +
    geom_ribbon(
      data = ggdata_ribbon, mapping = aes(x = x, ymin = ymin, ymax = ymax),
      alpha = 0.2, inherit.aes = FALSE) + 
    labs(x = element_blank(), y = element_blank()) + 
    theme_bw() + 
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),  
          text = element_text(size = text_size), 
          panel.background = element_rect(fill = "transparent"), 
          plot.background = element_rect(fill = "transparent", color = NA), 
          legend.background = element_rect(fill = "transparent"), 
          legend.box.background = element_rect(
            fill = "transparent", color = "transparent"))
  
  pltw
  ggsave(
    filename = paste0("scen", scenario, "_seqmedex_hist.pdf"),
    plot = plt0,
    width = 1.5, height = 2, units = c("in")
  )
  ggsave(
    filename = paste0("scen", scenario, "_seqmedex_plot.pdf"),
    plot = plt1,
    width = 1.5, height = 2, units = c("in")
  )
  ggsave(
    filename = paste0("scen", scenario, "_seqmedex_wass.pdf"),
    plot = pltw,
    width = 2.75, height = 2, units = c("in")
  )
} 

# plot the MSE of beta-hat #####################################################

if(!scenario1_illustration){
  source(paste(functions_dir, "/posterior_mean_mse.R", sep = ""))
  
  sim.idx = 1
  sm = seqmed_sims[[sim.idx]]
  bh = boxhill_sims[[sim.idx]]
  
  if(scenario == 1){
    muT = mu1
    VT = V1
    typeT = 3
  } else{
    muT = rep(0, 4)
    VT = diag(rep(sigmasq01, length(muT)))
    typeT = 4
  }
  MSEbetahat_doptlin = getMSEBeta(
    dopt_linear, Nttl, betaT, muT, VT, sigmasq, typeT)$MSE_postmean
  MSEbetahat_doptquad = getMSEBeta(
    dopt_quadratic, Nttl, betaT, muT, VT, sigmasq, typeT)$MSE_postmean
  MSEbetahat_space = getMSEBeta(
    space_filling, Nttl, betaT, muT, VT, sigmasq, typeT)$MSE_postmean
  MSEbetahat_seqmed = getMSEBeta(
    c(sm$x.in, sm$x.new), Nttl, betaT, muT, VT, sigmasq, typeT)$MSE_postmean
  MSEbetahat_bh = getMSEBeta(
    c(bh$x.in, bh$x.new), Nttl, betaT, muT, VT, sigmasq, typeT)$MSE_postmean
  
  b0 = c(MSEbetahat_doptlin[1], MSEbetahat_doptquad[1], MSEbetahat_space[1],
         MSEbetahat_seqmed[1], MSEbetahat_bh[1])
  b1 = c(MSEbetahat_doptlin[2], MSEbetahat_doptquad[2], MSEbetahat_space[2],
         MSEbetahat_seqmed[2], MSEbetahat_bh[2])
  b2 = c(MSEbetahat_doptlin[3], MSEbetahat_doptquad[3], MSEbetahat_space[3],
         MSEbetahat_seqmed[3], MSEbetahat_bh[3])
  
  if(scenario == 1){
    ggdata = data.frame(
      Designs = rep(c(
        "(v) DOptLin.", "(vi) DOptQuadr.", "(vii) Grid", "(ii) SeqMED", "(iv) BoxHill"), 3),
      MSE = c(b0, b1, b2), beta = rep(c("B0", "B1", "B2"), each = length(b0)))
    mseb.plt = ggplot(ggdata, aes(x = Designs, y = MSE)) +
      geom_bar(stat = "identity") +
      facet_wrap(vars(beta), scales = "free_y") +
      theme_bw() +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            axis.text.x = element_text(angle = 90, vjust = 0.5)) +
      labs(y = NULL)
  } else if(scenario == 2){
    b3 = c(MSEbetahat_doptlin[4], MSEbetahat_doptquad[4], MSEbetahat_space[4],
           MSEbetahat_seqmed[4], MSEbetahat_bh[4])
    
    ggdata = data.frame(
      Designs = rep(c(
        "(v) DOptLin.", "(vi) DOptQuadr.", "(vii) Grid", "(ii) SeqMED", "(iv) BoxHill"), 4),
      MSE = c(b0, b1, b2, b3), beta = rep(c("B0", "B1", "B2", "B3"), each = length(b0)))
    mseb.plt = ggplot(ggdata, aes(x = Designs, y = MSE)) +
      geom_bar(stat = "identity") +
      facet_wrap(vars(beta), scales = "free_y", ncol = 4) +
      theme_bw() +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            axis.text.x = element_text(angle = 90, vjust = 0.5),  
            text = element_text(size = text_size)) +
      labs(y = NULL)
  }
  
  mseb.plt
  ggsave(
    filename = paste0("scen", scenario, "_mseb_sim", sim.idx, ".pdf"),
    plot = mseb.plt,
    width = 6, height = 2.5, units = c("in")
  )
}

# plot the MSE of y-hat ########################################################

if(!scenario1_illustration){
  source(paste(functions_dir, "/predictive_yhat_mse.R", sep = ""))
  
  x_seq2 = seq(from = -1.25, to = 1.25, length.out = 1e4)
  
  yhatmse_space = getMSEYhat_seq(
    x_seq2, space_filling, Nttl, betaT, typeT, muT, VT, sigmasq, typeT)
  yhatmse_doptquad = getMSEYhat_seq(
    x_seq2, dopt_quadratic, Nttl, betaT, typeT, muT, VT, sigmasq, typeT)
  yhatmse_doptlin = getMSEYhat_seq(
    x_seq2, dopt_linear, Nttl, betaT, typeT, muT, VT, sigmasq, typeT)
  yhatmse_seqmed = getMSEYhat_seq(
    x_seq2, c(sm$x.in, sm$x.new),
    Nttl, betaT, typeT, muT, VT, sigmasq, typeT)
  yhatmse_bh = getMSEYhat_seq(
    x_seq2, c(bh$x.in, bh$x.new),
    Nttl, betaT, typeT, muT, VT, sigmasq, typeT)
  
  if(scenario == 1){
    ylimarg = range(
      0, yhatmse_space$MSEyhat, yhatmse_doptquad$MSEyhat, yhatmse_seqmed$MSEyhat,
      yhatmse_bh$MSEyhat)
  } else{
    ylimarg = c(0, 0.125)
  }
  
  ggdata = data.table(
    x = x_seq2,
    `(v) DOptLin.` = yhatmse_doptlin$MSEyhat,
    `(vi) DOptQuadr.` = yhatmse_doptquad$MSEyhat,
    `(vii) Grid` = yhatmse_space$MSEyhat,
    `(ii) SeqMED` = yhatmse_seqmed$MSEyhat,
    `(iv) BoxHill` = yhatmse_bh$MSEyhat
  )
  ggdata = melt(ggdata, id = c("x"), value.name = "yhatmse", variable.name = "Design")
  ggdata$Design = factor(
    ggdata$Design, 
    levels = c("(ii) SeqMED", "(iv) BoxHill", "(v) DOptLin.", "(vi) DOptQuadr.", "(vii) Grid"))
  ggdata2 = ggdata %>% filter(x %in% x_seq2[c(1, (1:10) * 1000)])
  msey.plt = ggplot(ggdata, aes(x = x, y = yhatmse, color = Design, linetype = Design)) +
    coord_cartesian(ylim = ylimarg, xlim = c(-1, 1)) +
    scale_x_continuous(breaks = c(-1, 0, 1)) +
    geom_path() +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),  
          text = element_text(size = text_size)) +
    labs(y = element_blank(), x = element_blank()) + 
    geom_point(data = ggdata2, 
               aes(x = x, y = yhatmse, color = Design, shape = Design), 
               size = 2)
  
  msey.plt
  if(scenario == 1){
    msey_scen1 = msey.plt + theme(legend.position = "none")
    ggsave(
      filename = paste0("scen", scenario, "_mseyhat_sim", sim.idx, ".pdf"),
      plot = msey_scen1,
      width = 2.2, height = 2.5, units = c("in")
    )
  } else if(scenario == 2){
    msey_scen2 = msey.plt 
    ggsave(
      filename = paste0("scen", scenario, "_mseyhat_sim", sim.idx, ".pdf"),
      plot = msey_scen2,
      width = 3.8, height = 2.5, units = c("in")
    )
  }
}







