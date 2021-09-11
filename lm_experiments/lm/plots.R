################################################################################
# last updated: 09/02/21
# purpose: to create a list of seqmed simulations
# scenario 1:
#   linear vs. quadratic,
#   where the true function is quadratic
# scenario 2:
#   linear vs. quadratic,
#   where the true function is cubic

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

# for box-hill design
source(paste(functions_dir, "/boxhill.R", sep = ""))

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
numSims = 100
numSeq = 100
seqN = 1
Nttl = numSeq * seqN
xmin = -1
xmax = 1
numCandidates = 10^3 + 1
candidates = seq(from = xmin, to = xmax, length.out = numCandidates)
sigmasq = 0.1

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
# import box & hill and seqmed simulations
################################################################################

seqmeds = readRDS(paste0(
  output_dir, "/scenario", scenario, 
  "_seqmed", 
  "_N", Nttl, 
  "_seed", rng.seed,
  ".rds"))
boxhills = readRDS(paste0(
  output_dir, "/scenario", scenario, 
  "_boxhill", 
  "_N", Nttl, 
  "_seed", rng.seed,
  ".rds"))

################################################################################
# non-sequential designs
################################################################################

space_filling =  seq(from = xmin, to = xmax, length.out = Nttl)
dopt_linear = c(rep(1, floor(Nttl / 2)), rep(-1, Nttl - floor(Nttl / 2)))
dopt_quadratic = c(rep(1, floor(Nttl / 3)), 
                   rep(0, ceiling(Nttl / 3)), 
                   rep(-1, Nttl - floor(Nttl / 3) - ceiling(Nttl / 3)))



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
sm = seqmeds[[sim.idx]]
ggdata = data.frame(x = c(sm$x.in, sm$x.new), y = c(sm$y.in, sm$y.new))
plt0 = ggplot(ggdata) + 
  geom_histogram(binwidth = 0.12, closed = "right", 
                 aes(x = x, y = after_stat(density))) + 
  theme_bw() + #base_size = 20) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
plt1 = ggplot(ggdata) + 
  geom_point(aes(x, y)) +
  stat_function(fun = fT) + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggarrange(plt0, plt1)

# plot a boxhill
bh = boxhills[[sim.idx]]
ggdata2 = data.frame(x = c(bh$x.in, bh$x.new), y = c(bh$y.in, bh$y.new))
plt2 = ggplot(ggdata2) + 
  geom_histogram(binwidth = 0.12, closed = "right", 
                 aes(x = x, y = after_stat(density))) + 
  theme_bw() + #base_size = 20) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
plt3 = ggplot(ggdata2) + 
  geom_point(aes(x, y)) +
  stat_function(fun = fT) + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# ggarrange(plt2, plt3)

# side-by-sides of seqmed and bh
# ggarrange(
#   plt0, plt1, # seqmed plots
#   plt2, plt3, # boxhill plots
#   nrow = 2, ncol = 2)

################################################################################
# plot the posterior probabilities of the hypotheses
################################################################################
source(paste(functions_dir, "/postprob_hypotheses.R", sep = ""))
library(data.table)

#######################################################################################################
if(scenario == 1){
  models = list("H0" = list(mu0, V0, 2),
                "H1" = list(mu1, V1, 3))
  typeT = 3
} else if(scenario == 2){
  models = list("H0" = list(mu0, V0, 2),
                "H1" = list(mu1, V1, 3),
                "H2" = list(betaT, diag(rep(sigmasq01, 4)), 4))
  typeT = 4
}

# non-sequential methods
epphs_space = calcEPPH(
  space_filling, Nttl, betaT, typeT, models, sigmasq, numSims = 100, seed = 123)
epphs_dopt1 = calcEPPH(
  dopt_linear, Nttl, betaT, typeT, models, sigmasq, numSims = 100, seed = 123)
epphs_dopt2 = calcEPPH(
  dopt_quadratic, Nttl, betaT, typeT, models, sigmasq, numSims = 100, seed = 123)

numModels = length(models)

# seqmed
pphs_seqmed= array(NA, dim = c(numModels, numSeq, numSims))
for(k in 1:numSims){
  smmed_data_k = seqmeds[[k]]
  for(i in 1:numSeq){
    if(i == 1){
      pphs_seqmed[ , i, k] = calcEPPHdata(
        smmed_data_k$y.in, 
        smmed_data_k$x.in, 
        N = seqN, models, sigmasq)
    } else{
      pphs_seqmed[ , i, k] = calcEPPHdata(
        c(smmed_data_k$y.in, smmed_data_k$y.new[1:(seqN * (i - 1))]), 
        c(smmed_data_k$x.in, smmed_data_k$x.new[1:(seqN * (i - 1))]), 
        N = seqN + seqN * (i - 1), models, sigmasq)
    }
  }
}
epphs_seqmed = apply(pphs_seqmed, c(1,2), mean)

# box-hill
pphs_bh = array(NA, dim = c(numModels, numSeq, numSims))
for(k in 1:numSims){
  bh_data_k = boxhills[[k]]
  for(i in 1:numSeq){
    if(i == 1){
      pphs_bh[ , i, k] = calcEPPHdata(
        bh_data_k$y.in, 
        bh_data_k$x.in, 
        N = seqN, models, sigmasq)
    } else{
      pphs_bh[ , i, k] = calcEPPHdata(
        c(bh_data_k$y.in, bh_data_k$y.new[1:(seqN * (i - 1))]), 
        c(bh_data_k$x.in, bh_data_k$x.new[1:(seqN * (i - 1))]), 
        N = seqN + seqN * (i - 1), models, sigmasq)
    }
  }
}
epphs_bh = apply(pphs_bh, c(1,2), mean)

# plot 1
ggdata0 = data.table(
  x = 1:numSeq, 
  Dlinear = rep(epphs_dopt1[1], numSeq), 
  Dquadratic = rep(epphs_dopt2[1], numSeq), 
  SpaceFilling = rep(epphs_space[1], numSeq), 
  SeqMED = epphs_seqmed[ 1, ], 
  BoxHill = epphs_bh[ 1, ], 
  Hypothesis = rep("H0", numSeq)
)
ggdata1 = data.table(
  x = 1:numSeq, 
  Dlinear = rep(epphs_dopt1[2], numSeq), 
  Dquadratic = rep(epphs_dopt2[2], numSeq), 
  SpaceFilling = rep(epphs_space[2], numSeq), 
  SeqMED = epphs_seqmed[ 2, ], 
  BoxHill = epphs_bh[ 2, ], 
  Hypothesis = rep("H1", numSeq)
)
if(scenario == 1){
  ggdata = rbind(ggdata0, ggdata1)
} else{
  ggdataT = data.table(
    x = 1:numSeq, 
    Dlinear = rep(epphs_dopt1[3], numSeq), 
    Dquadratic = rep(epphs_dopt2[3], numSeq), 
    SpaceFilling = rep(epphs_space[3], numSeq), 
    SeqMED = epphs_seqmed[ 3, ], 
    BoxHill = epphs_bh[ 3, ], 
    Hypothesis = rep("HT", numSeq)
  )
  ggdata = rbind(ggdata0, ggdata1, ggdataT)
}
ggdata.melted = melt(ggdata, id = c("x", "Hypothesis"), value.name = "epph", 
                     variable.name = "Design")
epph.plt = ggplot(ggdata.melted, aes(x = x, y = epph, color = Design, linetype = Design)) +
  facet_wrap(vars(Hypothesis)) + 
  geom_path(size = 1) + 
  scale_linetype_manual(values=c(rep("dashed", 3), rep("solid", 2))) + 
  geom_point(data = ggdata.melted[ggdata.melted$x == numSeq, ], 
             aes(x = x, y = epph), size = 2) + 
  theme_bw() +#base_size = 20) + 
  theme(panel.grid.minor = element_blank()) + 
  labs(y = "", x = "Stages") + 
  ylim(0, 1)
epph.plt
#######################################################################################################

models2 = list(model0, model1)

getPPH = function(
  design, models, n, true.function, error.var, randomize.order = FALSE, 
  seed = NULL
){
  ysims = sapply(design$x, FUN = simulateY_fromfunction(
    x = x, true.function = true.function, error.var = error.var))
  
  model.postprobs = matrix(NA, length(models), numSims)
  rownames(model.postprobs) = paste("model", 1:length(models), sep = "")
  colnames(model.postprobs) = paste("sim", 1:length(models), sep = "")
  for(j in 1:numSims){
    y.tmp = ysims[ , j]
    
    # calculate posterior probabilities for each model
    model.evidences.tmp = rep(NA, length(models))
    # get model evidences
    for(m in 1:length(models)){
      model.tmp = models[[m]]
      model.evidences.tmp[m] = Evidence_lm(
        y = y.tmp, x = design$x, model = model.tmp, error.var = error.var)
    }
    # get each hypotheses' posterior probability
    model.postprobs[, j] = getHypothesesPosteriors(
      prior.probs = rep(1 / length(models), length(models)), 
      evidences = model.evidences.tmp)
  }
  return(model.postprobs)
}

getPPHseq = function(
  design, models, n, true.function, error.var, randomize.order = FALSE, 
  seed = NULL
){
  x.new.idx = design$x.new.idx
  x.new = design$x.new
  y.new = design$y.new
  if(n != length(y.new)) warning("getPPHseq: n argument does not match length of new data")
  len.tmp = length(as.vector(na.omit(y.new)))
  if(randomize.order){
    new.order = sample(1:len.tmp, len.tmp, replace = FALSE)
    x.new.idx = x.new.idx[new.order]
    x.new = x.new[new.order]
    y.new = y.new[new.order]
  }
  ############### stopped here.... #####################################################
  for(j in 1:numSims){
    for(i in 1:n){
      
    }
  }
  
  
  
  ysims = sapply(design$x, FUN = simulateY_fromfunction(
    x = x, true.function = true.function, error.var = error.var))
  
  model.postprobs = matrix(NA, length(models), numSims)
  rownames(model.postprobs) = paste("model", 1:length(models), sep = "")
  colnames(model.postprobs) = paste("sim", 1:length(models), sep = "")
  for(j in 1:numSims){
    y.tmp = ysims[ , j]
    
    # calculate posterior probabilities for each model
    model.evidences.tmp = rep(NA, length(models))
    # get model evidences
    for(m in 1:length(models)){
      model.tmp = models[[m]]
      model.evidences.tmp[m] = Evidence_lm(
        y = y.tmp, x = design$x, model = model.tmp, error.var = error.var)
    }
    # get each hypotheses' posterior probability
    model.postprobs[, j] = getHypothesesPosteriors(
      prior.probs = rep(1 / length(models), length(models)), 
      evidences = model.evidences.tmp)
  }
  return(model.postprobs)
}


# seqmed
pphs_seqmed= array(NA, dim = c(numModels, numSeq, numSims))
for(k in 1:numSims){
  smmed_data_k = seqmeds[[k]]
  for(i in 1:numSeq){
    if(i == 1){
      
      pphs_seqmed[ , i, k] = calcEPPHdata(
        smmed_data_k$y.in, 
        smmed_data_k$x.in, 
        N = seqN, models, sigmasq)
    } else{
      pphs_seqmed[ , i, k] = calcEPPHdata(
        c(smmed_data_k$y.in, smmed_data_k$y.new[1:(seqN * (i - 1))]), 
        c(smmed_data_k$x.in, smmed_data_k$x.new[1:(seqN * (i - 1))]), 
        N = seqN + seqN * (i - 1), models, sigmasq)
    }
  }
}
epphs_seqmed = apply(pphs_seqmed, c(1,2), mean)

# box-hill
pphs_bh = array(NA, dim = c(numModels, numSeq, numSims))
for(k in 1:numSims){
  bh_data_k = boxhills[[k]]
  for(i in 1:numSeq){
    if(i == 1){
      pphs_bh[ , i, k] = calcEPPHdata(
        bh_data_k$y.in, 
        bh_data_k$x.in, 
        N = seqN, models, sigmasq)
    } else{
      pphs_bh[ , i, k] = calcEPPHdata(
        c(bh_data_k$y.in, bh_data_k$y.new[1:(seqN * (i - 1))]), 
        c(bh_data_k$x.in, bh_data_k$x.new[1:(seqN * (i - 1))]), 
        N = seqN + seqN * (i - 1), models, sigmasq)
    }
  }
}
epphs_bh = apply(pphs_bh, c(1,2), mean)

# plot 1
ggdata0 = data.table(
  x = 1:numSeq, 
  Dlinear = rep(epphs_dopt1[1], numSeq), 
  Dquadratic = rep(epphs_dopt2[1], numSeq), 
  SpaceFilling = rep(epphs_space[1], numSeq), 
  SeqMED = epphs_seqmed[ 1, ], 
  BoxHill = epphs_bh[ 1, ], 
  Hypothesis = rep("H0", numSeq)
)
ggdata1 = data.table(
  x = 1:numSeq, 
  Dlinear = rep(epphs_dopt1[2], numSeq), 
  Dquadratic = rep(epphs_dopt2[2], numSeq), 
  SpaceFilling = rep(epphs_space[2], numSeq), 
  SeqMED = epphs_seqmed[ 2, ], 
  BoxHill = epphs_bh[ 2, ], 
  Hypothesis = rep("H1", numSeq)
)
if(scenario == 1){
  ggdata = rbind(ggdata0, ggdata1)
} else{
  ggdataT = data.table(
    x = 1:numSeq, 
    Dlinear = rep(epphs_dopt1[3], numSeq), 
    Dquadratic = rep(epphs_dopt2[3], numSeq), 
    SpaceFilling = rep(epphs_space[3], numSeq), 
    SeqMED = epphs_seqmed[ 3, ], 
    BoxHill = epphs_bh[ 3, ], 
    Hypothesis = rep("HT", numSeq)
  )
  ggdata = rbind(ggdata0, ggdata1, ggdataT)
}
ggdata.melted = melt(ggdata, id = c("x", "Hypothesis"), value.name = "epph", 
                     variable.name = "Design")
epph.plt = ggplot(ggdata.melted, aes(x = x, y = epph, color = Design, linetype = Design)) +
  facet_wrap(vars(Hypothesis)) + 
  geom_path(size = 1) + 
  scale_linetype_manual(values=c(rep("dashed", 3), rep("solid", 2))) + 
  geom_point(data = ggdata.melted[ggdata.melted$x == numSeq, ], 
             aes(x = x, y = epph), size = 2) + 
  theme_bw() +#base_size = 20) + 
  theme(panel.grid.minor = element_blank()) + 
  labs(y = "", x = "Stages") + 
  ylim(0, 1)
epph.plt



################################################################################
# plot the posterior mean curve and the true curve when scenario == 2
################################################################################

if(scenario == 2){
  
}







################################################################################
# plot the MSE of beta-hat (posterior mean) of the hypotheses
################################################################################
source(paste(functions_dir, "/posterior_mean_mse.R", sep = ""))

if(scenario == 1){
  muT = mu1
  VT = V1
} else{
  muT = rep(0, 4)
  VT = diag(rep(sigmasq01, length(muT)))
}
MSEbetahat_doptlin = getMSEBeta(
  dopt_linear, Nttl, betaT, muT, VT, sigmasq, typeT)$MSE_postmean
MSEbetahat_doptquad = getMSEBeta(
  dopt_quadratic, Nttl, betaT, muT, VT, sigmasq, typeT)$MSE_postmean
MSEbetahat_space = getMSEBeta(
  space_filling, Nttl, betaT, muT, VT, sigmasq, typeT)$MSE_postmean
MSEbetahat_seqmed = getMSEBeta(
  c(seqmeds[[sim.idx]]$x.in, seqmeds[[sim.idx]]$x.new), 
  Nttl, betaT, muT, VT, sigmasq, typeT)$MSE_postmean
MSEbetahat_bh = getMSEBeta(
  c(boxhills[[sim.idx]]$x.in, boxhills[[sim.idx]]$x.new), 
  Nttl, betaT, muT, VT, sigmasq, typeT)$MSE_postmean

b0 = c(MSEbetahat_doptlin[1], MSEbetahat_doptquad[1], MSEbetahat_space[1], 
       MSEbetahat_seqmed[1], MSEbetahat_bh[1])
b1 = c(MSEbetahat_doptlin[2], MSEbetahat_doptquad[2], MSEbetahat_space[2], 
       MSEbetahat_seqmed[2], MSEbetahat_bh[2])
b2 = c(MSEbetahat_doptlin[3], MSEbetahat_doptquad[3], MSEbetahat_space[3], 
       MSEbetahat_seqmed[3], MSEbetahat_bh[3])

ggdata = data.frame(
  Designs = rep(c("Dlinear", "Dquadratic", "SpaceFilling", "SeqMED", "BoxHill"), 3), 
  MSE = c(b0, b1, b2), beta = rep(c("B0", "B1", "B2"), each = length(b0)))
mseb.plt = ggplot(ggdata, aes(x = Designs, y = MSE)) + 
  geom_bar(stat = "identity") +
  facet_wrap(vars(beta)) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.text.x = element_text(angle = 45, vjust = 0.5)) +
  labs(y = NULL)
mseb.plt

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
  x_seq2, c(seqmeds[[sim.idx]]$x.in, seqmeds[[sim.idx]]$x.new), 
  Nttl, betaT, typeT, muT, VT, sigmasq, typeT)
yhatmse_bh = getMSEYhat_seq(
  x_seq2, c(boxhills[[sim.idx]]$x.in, boxhills[[sim.idx]]$x.new), 
  Nttl, betaT, typeT, muT, VT, sigmasq, typeT)

if(scenario == 1){
  ylimarg = range(
    0, yhatmse_space$MSEyhat, yhatmse_doptquad$MSEyhat, yhatmse_seqmed$MSEyhat, 
    yhatmse_bh$MSEyhat)
} else{
  ylimarg = c(0, 0.15)
}

ggdata = data.table(
  x = x_seq2, 
  Dlinear = yhatmse_doptlin$MSEyhat, 
  Dquadratic = yhatmse_doptquad$MSEyhat, 
  SpaceFilling = yhatmse_space$MSEyhat, 
  SeqMED = yhatmse_seqmed$MSEyhat,
  BH = yhatmse_bh$MSEyhat
)
ggdata = melt(ggdata, id = c("x"), value.name = "yhatmse", variable.name = "Design")
msey.plt = ggplot(ggdata, aes(x = x, y = yhatmse, color = Design)) +
  coord_cartesian(ylim = ylimarg, xlim = c(-1, 1)) +
  geom_path() + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  labs(y = "", x = "x")
msey.plt
