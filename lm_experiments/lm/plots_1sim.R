################################################################################
# last updated: 2/28/22
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
gg_color_hue = function(n) {
  hues = seq(15, 275, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

################################################################################
# simulation settings, shared for both scenarios (linear vs. quadratic)
################################################################################

# simulations settings
numSims = 100 # 1 simulation with 100 design points
numSeq = 100 # 100 design points
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
# res_Fed_Doptlin$design

# consider quadratic model, y = b0 + b1 x + b2 x^2
res_Fed_Doptquad = optFederov(
  ~1+x+I(x^2), data = candidates_named, approximate = TRUE, criterion = "D")
# res_Fed_Doptquad$design

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
# plot the designs
################################################################################
# better with Nttl = 100
sim.idx = 5

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

# manuscript plot
# if(Nttl == 100){
#   ggsave(
#     filename = paste0("scen", scenario, "_seqmedexample.pdf"),
#     plot = last_plot(),
#     width = 4.5, height = 2, units = c("in")
#   )
# }

# # plot a boxhill
bh = boxhill_sims[[sim.idx]]
# ggdata2 = data.frame(x = c(bh$x.in, bh$x.new), y = c(bh$y.in, bh$y.new))
# plt2 = ggplot(ggdata2) +
#   geom_histogram(binwidth = 0.12, closed = "right",
#                  aes(x = x, y = after_stat(density))) +
#   theme_bw() + #base_size = 20) +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# plt3 = ggplot(ggdata2) +
#   geom_point(aes(x, y)) +
#   stat_function(fun = fT) +
#   theme_bw() +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# ggarrange(plt2, plt3)

################################################################################
# plot the wasserstein distance when scenario == 1
# plot the posterior mean curve and the true curve when scenario == 2
################################################################################
library(data.table)

if(scenario == 1){
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

  # # manuscript plot
  if(Nttl == 100){
    # ggsave(
    #   filename = paste0("scen", scenario, "_seqmedexamplewasserstein.pdf"),
    #   plot = last_plot(),
    #   width = 6.5, height = 1.75, units = c("in")
    # )
  }
} else if(scenario == 2){
  # plot ...
}

################################################################################
# plot the MSE of beta-hat (posterior mean) of the hypotheses
################################################################################
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
          axis.text.x = element_text(angle = 45, vjust = 0.5)) +
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
          axis.text.x = element_text(angle = 45, vjust = 0.5)) +
    labs(y = NULL)
}
mseb.plt

# # manuscript plot
# ggsave(
#   filename = paste0("scen", scenario, "_mseb_sim", sim.idx, ".pdf"),
#   plot = mseb.plt,
#   width = 6.5, height = 2.5, units = c("in")
# )

################################################################################
# plot the MSE of y-hat (posterior mean) of the hypotheses
################################################################################
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
msey.plt = ggplot(ggdata, aes(x = x, y = yhatmse, color = Design)) +
  coord_cartesian(ylim = ylimarg, xlim = c(-1, 1)) +
  geom_path() +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  labs(y = "", x = "x")
msey.plt
# # manuscript plot
# ggsave(
#   filename = paste0("scen", scenario, "_mseyhat_sim", sim.idx, ".pdf"),
#   plot = msey.plt,
#   width = 6.5, height = 2.5, units = c("in")
# )

if(scenario == 1){
  msey_scen1 = msey.plt + theme(legend.position = "none")
} else if(scenario == 2){
  msey_scen2 = msey.plt 
}
if(!is.null(msey_scen1) && !is.null(msey_scen2) && numSeq == 100){
  ggarrange(msey_scen1, msey_scen2, nrow = 1, ncol = 2, widths = c(1, 1.5))
  
  # manuscript plot
  ggsave(
    filename = paste0(
      "lm",
      "_mseyhats_sim", sim.idx, ".pdf"),
    plot = last_plot(),
    width = 6.5, height = 2, units = c("in")
  )
}


################################################################################
# all the alphas: 0, 0.5, 1, 2, 3, 4
################################################################################

# when N=12, alpha=50 leads to something like Box & Hill
alphas = c(0, 1, 10, 25, 50, 100)

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
# plot the designs

sim.idx = 1
seqmed_designs_alphas = data.frame(alpha = c(), x = c(), y = c())
for(i in 1:length(seqmed_sims_alphas)){
  sm.tmp = seqmed_sims_alphas[[i]][[sim.idx]]
  sm_design_alpha.tmp = data.frame(
    alpha = alphas[i],
    avg_min_alpha_used = paste0(
      min_alpha_used[i, sim.idx], " (mean = ", avg_min_alpha_used[i], ")"),
    x = c(sm.tmp$x.in, sm.tmp$x.new),
    y = c(sm.tmp$y.in, sm.tmp$y.new)
  )
  seqmed_designs_alphas = rbind(seqmed_designs_alphas, sm_design_alpha.tmp)
}
seqmed_designs_alphas$alpha = factor(seqmed_designs_alphas$alpha)

plt_alphas2 = ggplot(seqmed_designs_alphas) +
  facet_wrap(vars(alpha, avg_min_alpha_used)) +
  geom_histogram(binwidth = 0.12, closed = "right",
                 aes(x = x)) +#, y = after_stat(density))) +
  theme_bw() + #base_size = 20) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
plt_alphas2

# manuscript plot
ggsave(
  filename = paste0(
    "lm", "_scen", scenario, Dinit_label,
    "_designs_alphas", ".pdf"),
  plot = last_plot(),
  width = 6.5, height = 3.5, units = c("in")
)

