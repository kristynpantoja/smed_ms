rm(list=ls())
################################################################################
# last updated: 11/02/21
# purpose: to create a list of seqmed simulations
# dimT = 2:
#   dimensions (1, 2) vs dimensions (1, 2, 3)
#   where the true dimensions are (1, 2)
# dimT = 3:
#   dimensions (1, 2) vs dimensions (1, 2, 3)
#   where the true dimensions are (1, 2, 3)

dimT = 3 # 2, 3
sigmasq = 0.3 # 0.1, 0.3

################################################################################
# Sources/Libraries
################################################################################
output_dir = "lm_experiments/lmvs/outputs"
functions_dir = "functions"

# for seqmed design
source(paste(functions_dir, "/SeqMEDvs.R", sep = ""))
source(paste(functions_dir, "/SeqMEDvs_batch.R", sep = ""))
source(paste(functions_dir, "/charge_function_q.R", sep = ""))
source(paste(functions_dir, "/construct_design_matrix.R", sep = ""))
source(paste(functions_dir, "/wasserstein_distance.R", sep = ""))
source(paste(functions_dir, "/posterior_parameters.R", sep = ""))
source(paste(functions_dir, "/simulate_y.R", sep = ""))

# for generating initial data
source(paste(functions_dir, "/variance_marginal_y.R", sep = ""))

# for box-hill design
source(paste(functions_dir, "/boxhill.R", sep = ""))
source(paste(functions_dir, "/kl_divergence.R", sep = ""))

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
# simulation settings, shared for both scenarios
################################################################################

# simulations settings
numSims = 100
Nin = 1 # 1, 5
numSeq = 27
seqN = 1
Nnew = numSeq * seqN
Nttl = Nin + Nnew 
xmin = -1
xmax = 1
dimX = 3
numCandidates = 5000
x_seq = seq(from = xmin, to = xmax, length.out = floor((numCandidates)^(1 / 3)))
candidates = as.matrix(expand.grid(x_seq, x_seq, x_seq))

# hypothesis settings
mu_full = c(0.5, 0.5, 0.5)
indices0 = c(1, 2)
indices1 = 1:length(mu_full)
mu0 = rep(0, length(indices0))
mu1 = rep(0, length(indices1))
sigmasq01 = 0.25 # 0.25 - cannot distinguish dopt and 3fact when Nin = 5
V0 = diag(rep(sigmasq01,length(mu0)))
V1 = diag(rep(sigmasq01,length(mu1)))
model0 = list(
  indices = indices0, beta.mean = mu0, beta.var = V0)
model1 = list(
  indices = indices1, beta.mean = mu1, beta.var = V1)

# seqmed settings
p = dimX
k = 4 * p

# boxhill settings
prior_probs = rep(1 / 2, 2)

################################################################################
# Scenarios
################################################################################
if(dimT == 2){
  betaT = mu_full[indices0]
  indicesT = indices0
  fT = function(x) x[, indicesT, drop = FALSE] %*% betaT
} else if(dimT == 3){
  betaT = mu_full[indices1]
  indicesT = indices1
  fT = function(x) x[, indicesT, drop = FALSE] %*% betaT
}

################################################################################
# import seqmed simulations
################################################################################

seqmed_sims = readRDS(paste0(
  output_dir, "/3dim", 
  "_dim", dimT, 
  "_seqmed", 
  "_Nttl", Nttl,
  "_Nin", Nin,
  "_numSeq", numSeq,
  "_seqN", seqN,
  "_seed", rng.seed,
  "_noise", strsplit(as.character(sigmasq), split = "\\.")[[1]][2],
  ".rds"))

################################################################################
# non-sequential designs
################################################################################

# use AlgDesign to obtain Doptimal designs #

# consider linear model y = b1 x1 + b2 x2 + b3 x3
candidates_named = candidates
colnames(candidates_named) = paste("x", 1:dimX, sep = "")
res_Fed_Doptlin = optFederov(
  ~-1+x1+x2+x3, data = candidates_named, approximate = TRUE, criterion = "D")
res_Fed_Doptlin$design

# # consider linear model with 1st order interactions 
#   y = x1+x2+x3+x1x2+x1x3+x2x3
# res_Fed_Dopt1int = optFederov(
#   ~-1+x1+x2+x3+x1*x2+x1*x3+x2*x3, data = candidates, approximate = TRUE,
#   criterion = "D")
# # this is exactly like a 3-level factorial design

# define designs #

# random design
# see sims

# doptimal - linear
num_supportpts_Doptlin = nrow(res_Fed_Doptlin$design)
supportpt_assgnmt_Doptlin = cut(
  sample(1:Nnew, size = Nnew, replace = FALSE), # shuffle
  breaks = num_supportpts_Doptlin, labels = FALSE)
x_doptimal = matrix(NA, nrow = Nnew, ncol = dimX)
for(i in 1:Nnew){
  supportpt_tmp = as.numeric(res_Fed_Doptlin$design[
    supportpt_assgnmt_Doptlin[i], c("x1", "x2", "x3")])
  x_doptimal[i, ] = supportpt_tmp
}
# # check: ???
# res_Fed_Doptlin$design
# # table(x_doptimal) / Nttl # ???

# 2 level factorial design
fact2_pts = as.matrix(expand.grid(c(-1, 1), c(-1, 1), c(-1, 1)))
fact2_pts = fact2_pts[sample(1:nrow(fact2_pts), replace = FALSE), ] # shuffle
x_2fact = matrix(NA, nrow = Nnew, ncol = dimX)
for(i in 0:(Nnew - 1)){
  x_2fact[ i + 1, ] = as.numeric(fact2_pts[ 1 + (i %% nrow(fact2_pts)), ])
}

# 3 level factorial design
fact3_pts = as.matrix(expand.grid(c(-1, 0, 1), c(-1, 0, 1), c(-1, 0, 1)))
fact3_pts = fact3_pts[sample(1:nrow(fact3_pts), replace = FALSE), ] # shuffle
x_3fact = fact3_pts

random_sims = list()
dopt_sims = list()
fact2_sims = list()
fact3_sims = list()
for(i in 1:numSims){
  # input data
  x.in.tmp = matrix(runif(
    n = dimX * Nin, min = xmin, max = xmax), nrow = Nin, ncol = dimX)
  y.in.tmp = simulateY_frommultivarfunction(
    x = x.in.tmp, true.function = fT, error.var = sigmasq)
  
  # random design sims
  x_random.tmp = matrix(runif(n = dimX * Nnew, min = xmin, max = xmax), 
                    nrow = Nnew, ncol = dimX)
  random_sims[[i]] = list(
    x.in = x.in.tmp, y.in = y.in.tmp,
    x.new = x_random.tmp,
    y.new = as.vector(simulateY_frommultivarfunction(
      x = x_random.tmp, true.function = fT, error.var = sigmasq))
    )
  
  # Dopt-linear sims
  dopt_sims[[i]] = list(
    x.in = x.in.tmp, y.in = y.in.tmp,
    x.new = x_doptimal,
    y.new = as.vector(simulateY_frommultivarfunction(
      x = x_doptimal, true.function = fT, error.var = sigmasq))
    )
  
  # 2-level factorial sims
  fact2_sims[[i]] = list(
    x.in = x.in.tmp, y.in = y.in.tmp,
    x.new = x_2fact,
    y.new = as.vector(simulateY_frommultivarfunction(
      x = x_2fact, true.function = fT, error.var = sigmasq))
    )
  
  # 3-level factorial sims
  fact3_sims[[i]] = list(
    x.in = x.in.tmp, y.in = y.in.tmp,
    x.new = x_3fact,
    y.new = as.vector(simulateY_frommultivarfunction(
      x = x_3fact, true.function = fT, error.var = sigmasq))
  )
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
name_of_design_to_plot = "2factorial"

numBreaks = 8
sm = seqmed_sims[[sim.idx]]
ra = random_sims[[sim.idx]]
dl = dopt_sims[[sim.idx]]
f2 = fact2_sims[[sim.idx]]
f3 = fact3_sims[[sim.idx]]

if(name_of_design_to_plot == "seqmed"){
  design.tmp = sm
} else if(name_of_design_to_plot == "random"){
  design.tmp = ra
} else if(name_of_design_to_plot == "doptimal"){
  design.tmp = dl
} else if(name_of_design_to_plot == "2factorial"){
  design.tmp = f2
} else if(name_of_design_to_plot == "3factorial"){
  design.tmp = f3
}

# maxcounts = rep(NA, length(mu_full))
# for(i in 1:length(mu_full)){
#   marginal = i
#   h = hist(design.tmp$x.new[ , marginal], breaks = numBreaks, plot = FALSE)
#   maxcounts[i] = max(h$counts)
# }

# plot marginals
marginals = matrix(NA, nrow = Nnew, ncol = dimX)
for(i in 1:ncol(marginals)) {
  marginals[, i] = design.tmp$x.new[ , i]
}
colnames(marginals) = paste("Variable", 1:3, sep = " ")
library(data.table)
marginals = as.data.table(marginals)
marginals.tall = melt(marginals, measure.vars = 1:3)
ggplot(marginals.tall, aes(x = value)) + 
  facet_wrap(vars(variable)) +
  geom_histogram(binwidth = 0.12, closed = "right") + #, 
                 # aes(y = after_stat(density))) + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) + 
  labs(x = "x")

# # plot marginals 9 points at a time
# Nnew.tmp = 9
# marginals.tmp = matrix(
#   NA, nrow = Nnew.tmp, ncol = dimX)
# for(i in 1:ncol(marginals.tmp)) {
#   marginals.tmp[, i] = design.tmp$x.new[1:Nnew.tmp , i]
# }
# colnames(marginals.tmp) = paste("Variable", 1:3, sep = " ")
# marginals.tmp = as.data.table(marginals.tmp)
# marginals.tall.tmp = melt(marginals.tmp, measure.vars = 1:3)
# ggplot(marginals.tall.tmp, aes(x = value)) + 
#   facet_wrap(vars(variable)) +
#   geom_histogram(binwidth = 0.12, closed = "right") +
#   theme_bw() + 
#   theme(panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank()) + 
#   labs(x = "x")

# # 3d scatterplot?
# library(scatterplot3d)
# scatterplot3d(marginals)

################################################################################
# plot the posterior probabilities of the hypotheses
################################################################################
library(data.table)

models = list(model0, model1)

# boxhill_vs methods that don't exist yet #
# get evidence
Evidence_lmvs = function(
  y, 
  x, 
  model, 
  error.var
){
  # make sure y is a vector
  if("matrix" %in% class(y)) y = as.vector(y)
  # get mean and variance of marginal density of y, which is n-dim multivariate normal pdf
  marg = getLMMarginal(
    x[, model$indices, drop = FALSE], model$beta.mean, model$beta.var, 
    error.var)
  evidence = dmvnorm(y, mean = marg$mean, sigma = marg$var)
  return(evidence)
}

# pph methods #
getPPH = function(
  design, models, true.function, error.var, initial.data = FALSE, seed = NULL
){
  if(!is.null(seed)) set.seed(seed)
  if(initial.data){
    x = rbind(design$x.in, design$x.new)
    y = c(design$y.in, design$y.new)
  } else{
    x = design$x
    y = design$y
  }
  
  # get model evidences
  model.evidences = rep(NA, length(models))
  for(m in 1:length(models)){
    model.tmp = models[[m]]
    model.evidences[m] = Evidence_lmvs(
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
  # define x and y, for which posterior probs will be sequentially computed
  if(initial.data){ # if there is initial data, just do seq calc for new data
    # in this case, n = # of new points
    x = design$x.new
    y = design$y.new
    if(length(y) != n) stop("n does not equal the number of new points")
    if(randomize.order){
      new.order = sample(1:n, n, replace = FALSE)
      x = x.new[new.order, ]
      y = y.new[new.order]
    }
  } else{ # otherwise, do seq calc for all data
    # in this case, n = total # of points
    x = design$x
    y = design$y
    if(length(y) != n) stop("n does not equal the length of y, i.e. total number of points")
    if(randomize.order){
      new.order = sample(1:n, n, replace = FALSE)
      x = x[new.order, ]
      y = y[new.order]
    }
  }
  
  # calculate sequential posterior probabilities of hypothesized models
  PPH_mat = matrix(NA, nrow = n, ncol = length(models))
  for(i in 1:n){
    x.tmp = x[1:i, , drop = FALSE]
    y.tmp = y[1:i]
    
    # get model evidences
    model.evidences.tmp = rep(NA, length(models))
    for(m in 1:length(models)){
      model.tmp = models[[m]]
      model.evidences.tmp[m] = Evidence_lmvs(
        y = y.tmp, x = x.tmp, model = model.tmp, error.var = error.var)
    }
    PPH_mat[i, ] = getHypothesesPosteriors(
      prior.probs = rep(1 / length(models), length(models)), 
      evidences = model.evidences.tmp)
  }
  # if there is initial data, add that to the matrix of PP calculations
  if(initial.data){
    x.tmp = design$x.in
    y.tmp = design$y.in
    
    # get model evidences
    model.evidences.tmp = rep(NA, length(models))
    for(m in 1:length(models)){
      model.tmp = models[[m]]
      model.evidences.tmp[m] = Evidence_lmvs(
        y = y.tmp, x = x.tmp, model = model.tmp, error.var = error.var)
    }
    PPH_mat = rbind(
      getHypothesesPosteriors(
        prior.probs = rep(1 / length(models), length(models)), 
        evidences = model.evidences.tmp), 
      PPH_mat)
  }
  # turn the matrix into a data frame
  colnames(PPH_mat) = paste("H", 0:(length(models) - 1), sep = "")
  PPH_mat = data.frame(PPH_mat)
  if(initial.data){
    PPH_mat$index = 0:n
  } else{
    PPH_mat$index = 1:n
  }
  return(PPH_mat)
}

#

# non-sequential designs
PPH_df = data.frame()
for(j in 1:numSims){
  # sequence of PPHs for each design
  PPH_random = getPPH(
    random_sims[[j]], models, true.function, sigmasq, TRUE)
  PPH_dopt = getPPH(
    dopt_sims[[j]], models, true.function, sigmasq, TRUE)
  PPH_fact2 = getPPH(
    fact2_sims[[j]], models, true.function, sigmasq, TRUE)
  PPH_fact3 = getPPH(
    fact3_sims[[j]], models, true.function, sigmasq, TRUE)
  # master data frame
  PPH_random$Design = "Random"
  PPH_dopt$Design = "DOptimal"
  PPH_fact2$Design = "2Factorial"
  PPH_fact3$Design = "3Factorial"
  PPH.tmp = rbind(PPH_random, PPH_dopt, PPH_fact2, PPH_fact3)
  PPH.tmp$sim = j
  PPH_df = rbind(PPH_df, PPH.tmp)
}
PPHmean = aggregate(
  PPH_df[, names(PPH_df)[1:length(models)]], 
  by = list(PPH_df[, "Design"]), FUN = function(x) mean(x, na.rm = TRUE))
names(PPHmean)[1] = "Design"
PPHmean$index = Nnew #Nttl
# but we want a line, so allow interpolation by setting PPHmean$index = 0 too
PPHmean2 = PPHmean
PPHmean2$index = 0
PPHmean = rbind(PPHmean, PPHmean2)

# sequential designs
PPH_seq = data.frame()
for(j in 1:numSims){
  # sequence of PPHs for each design
  PPH_seq.sm = getPPHseq(
    seqmed_sims[[j]], models, Nnew, fT, sigmasq, TRUE, FALSE) 
  # master data frame
  PPH_seq.sm$Design = "SeqMED"
  PPH_seq.tmp = rbind(PPH_seq.sm)
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
design_names = rev(c("SeqMED", "DOptimal", "3Factorial", "2Factorial", "Random"))
PPHmean_gg$Design = factor(PPHmean_gg$Design, levels = design_names)
PPHmean_gg = setorder(PPHmean_gg, cols = "Design")
PPHmean_gg2 = PPHmean_gg[PPHmean_gg$index == Nnew, ]
PPHmean_gg2$Design = factor(PPHmean_gg2$Design, levels = design_names)
PPHmean_gg2 = setorder(PPHmean_gg2, cols = "Design")
epph.plt = ggplot(PPHmean_gg, aes(x = index, y = value, color = Design,
                                  linetype = Design, shape = Design)) +
  facet_wrap(~hypothesis) +
  geom_path() +
  scale_linetype_manual(values=c(rep("dashed", 4), rep("solid", 2))) +
  geom_point(data = PPHmean_gg2, 
             mapping = aes(x = index, y = value, color = Design),
             inherit.aes = FALSE) +
  theme_bw() +
  ylim(0, 1) +
  labs(x = "Stage Index", y = element_blank())
plot(epph.plt)

# plot
ggsave(
  filename = paste0(
    "3dim", 
    "_dim", dimT, 
    "_noise", strsplit(as.character(sigmasq), split = "\\.")[[1]][2],
    "_epph.pdf"), 
  plot = epph.plt, 
  width = 6, height = 3, units = c("in")
)

################################################################################
# plot the MSE of beta-hat (posterior mean) of the hypotheses
################################################################################
# source(paste(functions_dir, "/posterior_mean_mse.R", sep = ""))
# 
# if(dimT == 2){
#   hyp_mu = mu0
#   hyp_V = V0
#   hyp_ind = indices0
# } else{
#   hyp_mu = mu1
#   hyp_V = V1
#   hyp_ind = indices1
# }
# mseBn_smmed = getMSEBeta(
#   rbind(sm$x.in, sm$x.new), 
#   Ntot, betaT, hyp_mu, hyp_V, sigmasq, NULL, hyp_ind)$MSE_postmean
# mseBn_rand = getMSEBeta(
#   rbind(ra$x.in, ra$x.new), 
#   Ntot, betaT, hyp_mu, hyp_V, sigmasq, NULL, hyp_ind)$MSE_postmean
# mseBn_dopt = getMSEBeta(
#   rbind(dl$x.in, dl$x.new), 
#   Ntot, betaT, hyp_mu, hyp_V, sigmasq, NULL, hyp_ind)$MSE_postmean
# mseBn_fact = getMSEBeta(
#   rbind(fa$x.in, fa$x.new), 
#   Ntot, betaT, hyp_mu, hyp_V, sigmasq, NULL, hyp_ind)$MSE_postmean
# 
# # plot
# b1 = c(mseBn_smmed[1], mseBn_rand[1], mseBn_dopt[1], mseBn_fact[1])
# b2 = c(mseBn_smmed[2], mseBn_rand[2], mseBn_dopt[2], mseBn_fact[2])
# ggdesigns = c("SeqMED", "Random", "Doptimal", "Factorial3")
# ggdata = data.frame(Designs = factor(rep(ggdesigns, 2), 
#                                      levels = ggdesigns[c(2, 1, 4, 3)]), 
#                     MSE = c(b1, b2), beta = rep(c("B1", "B2"), each = length(b1)))
# mse.plt = ggplot(ggdata, aes(x = Designs, y = MSE)) + 
#   geom_bar(stat = "identity") +
#   facet_wrap(vars(beta)) +
#   theme_bw() + 
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
#         axis.text.x = element_text(angle = 45, vjust = 0.5)) +
#   labs(y = NULL)
# mse.plt
# 
# # manuscript plot
# ggsave(
#   filename = paste0("dim", dimT, "_mseb.pdf"), 
#   plot = mse.plt, 
#   width = 4.5, height = 2, units = c("in")
# )
