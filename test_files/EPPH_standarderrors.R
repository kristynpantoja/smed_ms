################################################################################
# last updated: 12/17/20
# purpose: to calculate and plot metrics for scenario 1 simulations



################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
# LM Plots
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################



################################################################################
# Sources/Libraries
################################################################################
# tamu cluster
output_home = "run_designs/updated_simulations/lm"
functions_home = "functions"

# for seqmed design
source(paste(functions_home, "/SeqMED.R", sep = ""))
source(paste(functions_home, "/SeqMED_batch.R", sep = ""))
source(paste(functions_home, "/MMED.R", sep = ""))
source(paste(functions_home, "/variance_marginal_y.R", sep = ""))
source(paste(functions_home, "/charge_function_q.R", sep = ""))
source(paste(functions_home, "/construct_design_matrix.R", sep = ""))
source(paste(functions_home, "/wasserstein_distance.R", sep = ""))
source(paste(functions_home, "/posterior_parameters.R", sep = ""))
source(paste(functions_home, "/simulate_y.R", sep = ""))

# for evaluating designs
source(paste(functions_home, "/simulate_y.R", sep = ""))
source(paste(functions_home, "/postprob_hypotheses.R", sep = ""))
source(paste(functions_home, "/posterior_mean_mse.R", sep = ""))
source(paste(functions_home, "/predictive_yhat_mse.R", sep = ""))

library(expm)
library(matrixStats)
library(MASS)
library(mvtnorm)
library(knitr)

# for plots
library(ggplot2)
library(ggpubr)
library(reshape2)
library(data.table)
gg_color_hue = function(n) {
  hues = seq(15, 275, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
image_path = "plots"

################################################################################
# simulation settings, shared for both scenarios (linear vs. quadratic)
################################################################################

# simulations settings
numSims = 500
isParallelized = TRUE

# simulation settings
numSeq = 100
seqN = 1
N = numSeq * seqN
xmin = -1
xmax = 1
numCandidates = 10^3 + 1
candidates = seq(from = xmin, to = xmax, length.out = numCandidates)

# SeqMED settings
type01 = c(2, 3)
sigmasq = 0.1
mu0 = c(0, 0)
mu1 = c(0, 0, 0)
sigmasq01 = 0.25
V0 = diag(rep(sigmasq01,length(mu0)))
V1 = diag(rep(sigmasq01,length(mu1)))
f0 = function(x) mu0[1] + mu0[2] * x
f1 = function(x) mu1[1] + mu1[2] * x + mu1[3] * x^2

# boxhill settings
MMEDinputdata = FALSE
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
prior_probs = rep(1 / 2, 2)

################################################################################
# non-sequential designs
################################################################################

space_filling =  seq(from = xmin, to = xmax, length.out = N)
dopt_linear = c(rep(1, floor(N/2)), rep(-1, N - floor(N / 2)))
dopt_quadratic = c(rep(1, floor(N / 3)), 
                   rep(0, ceiling(N / 3)), 
                   rep(-1, N - floor(N / 3) - ceiling(N / 3)))

################################################################################
# Scenario 1: True function is quadratic
################################################################################
betaT = c(-0.2, -0.4, 0.4)
fT = function(x) betaT[1] + betaT[2] * x + betaT[3] * x^2

# seqmed settings
typeT = 3

# import the simulations
seqmed_file = paste0(output_home, "/seqmed/scenario1_seqmed_simulations", 
                     "_numSeq", numSeq, 
                     "_seqN", seqN,
                     "_numSims", numSims, 
                     ".rds")
bh_file = paste0(output_home, "/boxhill/scenario1_boxhill_simulations", 
                 "_N", N, 
                 "_MMEDinput", as.numeric(MMEDinputdata),
                 "_numSims", numSims, 
                 ".rds")
seqmeds = readRDS(seqmed_file)
bhs = readRDS(bh_file)

if(MMEDinputdata){
  # check if preliminary data are the same
  seqmed_x_input = c()
  seqmed_y_input = c()
  bh_x_input = c()
  bh_y_input = c()
  for(k in 1:numSims){
    seqmed.temp = seqmeds[[k]]
    bh_input.temp = bhs_inputs[[k]]
    seqmed_x_input = c(seqmed_x_input, seqmed.temp$D[1:seqN])
    seqmed_y_input = c(seqmed_y_input, seqmed.temp$y[1:seqN])
    bh_x_input = c(bh_x_input, bh_input.temp$x_input)
    bh_y_input = c(bh_y_input, bh_input.temp$y_input)
  }
  all.equal(seqmed_x_input, bh_x_input) # true
  all.equal(seqmed_y_input, bh_y_input) # true
}

################################################################################
# plot the designs
################################################################################
which.sim = 25

# plot a seqmed
seqmed1 = seqmeds[[which.sim]]

################################################################################
# plot the wasserstein distance and model fits
################################################################################
numseq = 1e2
x_seq = seq(from = xmin, to = xmax, length.out = numseq)

w_seq = sapply(x_seq, function(x) WNlm(
  x, seqmed1$postmean0[, numSeq], seqmed1$postmean1[, numSeq],
  diag(seqmed1$postvar0[, numSeq]), diag(seqmed1$postvar1[, numSeq]),
  sigmasq, type01))
f1est = function(x) seqmed1$postmean1[1, numSeq] +
  seqmed1$postmean1[2, numSeq] * x + seqmed1$postmean1[3, numSeq] * x^2
f2est = function(x) seqmed1$postmean0[1, numSeq] +
  seqmed1$postmean0[2, numSeq] * x
f1est_seq = sapply(x_seq, f1est)
f2est_seq = sapply(x_seq, f2est)
fT_seq = sapply(x_seq, fT)

################################################################################
# plot the posterior probabilities of the hypotheses
################################################################################

# first, calculate the posterior probabilities

# for the non-sequential methods;
epphs_space = calcExpPostProbH(
  space_filling, N, betaT, mu0, mu1, V0, V1, sigmasq, numSims, 
  typeT, type01, seed = 123, saveSims = T)
epphs_dopt1 = calcExpPostProbH(
  dopt_linear, N, betaT, mu0, mu1, V0, V1, sigmasq, numSims, 
  typeT, type01, seed = 123, saveSims = T)
epphs_dopt2 = calcExpPostProbH(
  dopt_quadratic, N, betaT, mu0, mu1, V0, V1, sigmasq, numSims, 
  typeT, type01, seed = 123, saveSims = T)
# save H0sims and H1sims
epphsH0sims = as.data.frame(
  cbind(
    SpaceFill = epphs_space$H0sims, 
    Dlinear = epphs_dopt1$H0sims, 
    Dquadratic = epphs_dopt2$H0sims
  )
)
epphsH1sims = as.data.frame(
  cbind(
    SpaceFill = epphs_space$H1sims, 
    Dlinear = epphs_dopt1$H1sims, 
    Dquadratic = epphs_dopt2$H1sims
  )
)

# seqmed posterior probabilities
postprobs0seq = matrix(NA, numSeq, numSims)
postprobs1seq = matrix(NA, numSeq, numSims)
BF01seq = matrix(NA, numSeq, numSims)
for(k in 1:numSims){
  smmed_data_k = seqmeds[[k]]
  for(i in 1:numSeq){
    changing_postprobs = calcExpPostProbH_data(
      smmed_data_k$y[1:(seqN * i)], smmed_data_k$D[1:(seqN * i)], 
      N = seqN * i, mu0, V0, mu1, V1, sigmasq, type01)
    postprobs0seq[i, k] = changing_postprobs[1]
    postprobs1seq[i, k] = changing_postprobs[2]
    BF01seq[i, k] = changing_postprobs[3]
  }
}
# save H0sims and H1sims (at last stage)
epphsH0sims$SeqMED = postprobs0seq[numSeq, ]
epphsH1sims$SeqMED = postprobs1seq[numSeq, ]
# get expected value (average)
epph0seq_seqmed = apply(postprobs0seq, 1, mean)
epph1seq_seqmed = apply(postprobs1seq, 1, mean)
BF01seq_seqmed = apply(BF01seq, 1, mean)

# box-hill posterior probabilities
# first, put it in the format
bhs.format = list()
for(k in 1:numSims){
  bh.temp = bhs[[k]]
  bhs.format[[k]] = list(
    D = c(bh.temp$x, bh.temp$x.new),
    y = c(bh.temp$y, bh.temp$y.new),
    post.probs = bh.temp$post.probs
  )
}

postprobs0seq = matrix(NA, N, numSims)
postprobs1seq = matrix(NA, N, numSims)
BF01seq = matrix(NA, N, numSims)
for(k in 1:numSims){
  bh_data_k = bhs.format[[k]]
  for(i in 1:N){
    changing_postprobs = suppressWarnings(calcExpPostProbH_data(
      bh_data_k$y[1:i], bh_data_k$D[1:i], 
      N = i, mu0, V0, mu1, V1, sigmasq, type01))
    postprobs0seq[i, k] = changing_postprobs[1]
    postprobs1seq[i, k] = changing_postprobs[2]
    BF01seq[i, k] = changing_postprobs[3]
  }
}
# save H0sims and H1sims (at last stage)
epphsH0sims$BoxHill = postprobs0seq[numSeq, ]
epphsH1sims$BoxHill = postprobs1seq[numSeq, ]
# get expected value (average)
epph0seq_bh2 = apply(postprobs0seq, 1, mean)
epph1seq_bh2 = apply(postprobs1seq, 1, mean)
BF01seq_bh2 = apply(BF01seq, 1, mean)
ggdata00 = data.table(
  x = 1:numSeq, 
  Dlinear = rep(epphs_dopt1[[1]], numSeq), 
  Dquadratic = rep(epphs_dopt2[[1]], numSeq), 
  SpaceFill = rep(epphs_space[[1]], numSeq), 
  SeqMED = epph0seq_seqmed,
  Hypothesis = rep("H0", numSeq), 
  key = "x"
)
ggdata01 = data.table(
  x = seq(1, numSeq, length.out = N), 
  BoxHill = epph0seq_bh2, 
  Hypothesis = rep("H0", N), 
  key = "x"
)
ggdata10 = data.table(
  x = 1:numSeq, 
  Dlinear = rep(epphs_dopt1[[2]], numSeq), 
  Dquadratic = rep(epphs_dopt2[[2]], numSeq), 
  SpaceFill = rep(epphs_space[[2]], numSeq), 
  SeqMED = epph1seq_seqmed,
  Hypothesis = rep("H1", numSeq), 
  key = "x"
)
ggdata11 = data.table(
  x = seq(1, numSeq, length.out = N), 
  BoxHill = epph1seq_bh2, 
  Hypothesis = rep("H1", N), 
  key = c("x")
)
ggdata.woBH = merge(ggdata00, ggdata10, all = TRUE, 
                    by = names(ggdata00))
ggdata.BH = merge(ggdata01, ggdata11, all = TRUE, 
                  by = names(ggdata01))
ggdata.tog = merge(ggdata.woBH, ggdata.BH, all = TRUE, 
                   by = c("x", "Hypothesis"))
ggdata.m2 = melt(ggdata.tog, id = c("x", "Hypothesis"), value.name = "epph", 
                 variable.name = "Design")
ggdata0.m2 = melt(ggdata.woBH, id = c("x", "Hypothesis"), value.name = "epph", 
                  variable.name = "Design")
ggdata1.m2 = melt(ggdata.BH, id = c("x", "Hypothesis"), value.name = "epph", 
                  variable.name = "Design")
ggplot(ggdata0.m2, aes(x = x, y = epph, color = Design, linetype = Design)) +
  facet_wrap(vars(Hypothesis)) + 
  geom_path() + 
  scale_linetype_manual(values=c("solid", "dashed", "dashed", "solid", "dashed")) +
  geom_point(data = ggdata0.m2[x == numSeq], aes(x = x, y = epph)) + 
  theme_bw() + 
  scale_x_continuous(breaks = 1:10) +
  theme(panel.grid.minor = element_blank(), axis.text.x = element_text(size = 5)) + 
  labs(y = "", x = "Stages") + 
  ylim(0, 1) + 
  geom_path(data = ggdata1.m2, mapping = aes(x = x, y = epph, color = Design, 
                                             linetype = Design), 
            inherit.aes = FALSE) + 
  geom_point(data = ggdata1.m2[x == numSeq], aes(x = x, y = epph))

# calculate standard errors
dim(epphsH0sims)
dim(epphsH1sims)
H0stats = data.frame(
  mean = apply(epphsH0sims, MARGIN = 2, FUN = mean),
  sd = apply(epphsH0sims, MARGIN = 2, FUN = sd),
  se = apply(epphsH0sims, MARGIN = 2, FUN = sd) / numSims
)
H0stats$lower = H0stats$mean - 2 * H0stats$se
H0stats$upper = H0stats$mean + 2 * H0stats$se
round(H0stats, 5)
H1stats = data.frame(
  mean = apply(epphsH1sims, MARGIN = 2, FUN = mean),
  sd = apply(epphsH1sims, MARGIN = 2, FUN = sd),
  se = apply(epphsH1sims, MARGIN = 2, FUN = sd) / numSims
)
H1stats$lower = H1stats$mean - 2 * H1stats$se
H1stats$upper = H1stats$mean + 2 * H1stats$se
round(H1stats, 5)

#
##
###
####
#####
######
#######
########
#########
##########
#########
########
#######
######
#####
####
###
##
#

################################################################################
# Scenario 2: True function is cubic
################################################################################
betaT = c(0, -0.75, 0, 1)
fT = function(x) betaT[1] + betaT[2] * x + betaT[3] * x^2 + betaT[4] * x^3

# seqmed settings
typeT = 4

# import the simulations
seqmed_file = paste0(output_home, "/seqmed/scenario2_seqmed_simulations", 
                     "_numSeq", numSeq, 
                     "_seqN", seqN,
                     "_numSims", numSims, 
                     ".rds")
bh_file = paste0(output_home, "/boxhill/scenario2_boxhill_simulations", 
                 "_N", N, 
                 "_MMEDinput", as.numeric(MMEDinputdata),
                 "_numSims", numSims, 
                 ".rds")
seqmeds = readRDS(seqmed_file)
bhs = readRDS(bh_file)

if(MMEDinputdata){
  # check if preliminary data are the same
  seqmed_x_input = c()
  seqmed_y_input = c()
  bh_x_input = c()
  bh_y_input = c()
  for(k in 1:numSims){
    seqmed.temp = seqmeds[[k]]
    bh_input.temp = bhs_inputs[[k]]
    seqmed_x_input = c(seqmed_x_input, seqmed.temp$D[1:seqN])
    seqmed_y_input = c(seqmed_y_input, seqmed.temp$y[1:seqN])
    bh_x_input = c(bh_x_input, bh_input.temp$x_input)
    bh_y_input = c(bh_y_input, bh_input.temp$y_input)
  }
  all.equal(seqmed_x_input, bh_x_input) # true
  all.equal(seqmed_y_input, bh_y_input) # true
}

################################################################################
# plot the designs
################################################################################
which.sim = 25

# plot a seqmed
seqmed1 = seqmeds[[which.sim]]

################################################################################
# plot the fits
################################################################################

# col_postmeans = c(1, 2, 3, 4)
seqmed1.postmean2 = matrix(NA, length(betaT), numSeq)
for(k in 1:numSeq){
  seqmed1.postmean2[ , k] = as.vector(postmean(seqmed1$y[1:(seqN * k)],
                                               seqmed1$D[1:(seqN * k)], (seqN * k),
                                               c(0, 0, 0, 0), diag(rep(sigmasq01, 4)), sigmasq, 4))
}

seqmed1.postmean2.1 = postmean(seqmed1$y, seqmed1$D, N,
                                   c(0, 0, 0, 0), diag(rep(sigmasq01, 4)), sigmasq, 4)
fTest = function(x){
  seqmed1.postmean2[1, numSeq] + 
    seqmed1.postmean2[2, numSeq] * x + 
    seqmed1.postmean2[3, numSeq] * x^2 + 
    seqmed1.postmean2[4, numSeq] * x^3
}
f2est = function(x){
  seqmed1$postmean1[1, numSeq] + 
    seqmed1$postmean1[2, numSeq] * x + 
    seqmed1$postmean1[3, numSeq] * x^2
}
f1est = function(x){
  seqmed1$postmean0[1, numSeq] + 
    seqmed1$postmean0[2, numSeq] * x
}

f1est_seq = sapply(x_seq, f1est)
f2est_seq = sapply(x_seq, f2est)
fTest_seq = sapply(x_seq, fTest)
fT_seq = sapply(x_seq, fT)

################################################################################
# plot the posterior probabilities of the hypotheses
################################################################################

# # first, calculate the posterior probabilities

models = list("H0" = list(mu0, V0, 2),
              "H1" = list(mu1, V1, 3),
              "H2" = list(betaT, diag(rep(sigmasq01, 4)), 4))
numModels = length(models)

# non-sequential methods
epphs_space = calcEPPH(
  space_filling, N, betaT, typeT, models, sigmasq, numSims, seed = 123, saveSims = TRUE)
epphs_dopt1 = calcEPPH(
  dopt_linear,  N, betaT, typeT, models, sigmasq, numSims, seed = 123, saveSims = TRUE)
epphs_dopt2 = calcEPPH(
  dopt_quadratic, N, betaT, typeT, models, sigmasq, numSims, seed = 123, saveSims = TRUE)
# save H0sims and H1sims
epphsH0sims = as.data.frame(
  cbind(
    SpaceFill = epphs_space$sims[1, ], 
    Dlinear = epphs_dopt1$sims[1, ], 
    Dquadratic = epphs_dopt2$sims[1, ]
  )
)
epphsH1sims = as.data.frame(
  cbind(
    SpaceFill = epphs_space$sims[2, ], 
    Dlinear = epphs_dopt1$sims[2, ], 
    Dquadratic = epphs_dopt2$sims[2, ]
  )
)
epphsHTsims = as.data.frame(
  cbind(
    SpaceFill = epphs_space$sims[3, ], 
    Dlinear = epphs_dopt1$sims[3, ], 
    Dquadratic = epphs_dopt2$sims[3, ]
  )
)
# seqmed
pphs_seqmed= array(NA, dim = c(numModels, numSeq, numSims))
for(k in 1:numSims){
  smmed_data_k = seqmeds[[k]]
  for(i in 1:numSeq){
    pphs_seqmed[ , i, k] = calcEPPHdata(smmed_data_k$y[1:(seqN * i)], 
                                        smmed_data_k$D[1:(seqN * i)], 
                                        N = seqN * i, models, sigmasq)
  }
}
# save H0sims and H1sims (at last stage)
epphsH0sims$SeqMED = pphs_seqmed[1, numSeq, ]
epphsH1sims$SeqMED = pphs_seqmed[2, numSeq, ]
epphsHTsims$SeqMED = pphs_seqmed[3, numSeq, ]
# get expected value (average)
epphs_seqmed = apply(pphs_seqmed, c(1,2), mean)

# box-hill
bhs.format = list()
for(k in 1:numSims){
  bh.temp = bhs[[k]]
  bhs.format[[k]] = list(
    D = c(bh.temp$x, bh.temp$x.new),
    y = c(bh.temp$y, bh.temp$y.new),
    post.probs = bh.temp$post.probs
  )
}
pphs_bh = array(NA, dim = c(numModels, N, numSims))
for(k in 1:numSims){
  bh_data_k = bhs.format[[k]]
  for(i in 1:N){
    pphs_bh[ , i, k] = calcEPPHdata(bh_data_k$y[1:i], 
                                    bh_data_k$D[1:i], 
                                    N = i, models, sigmasq)
  }
}
# save H0sims and H1sims (at last stage)
epphsH0sims$BoxHill = pphs_bh[1, numSeq, ]
epphsH1sims$BoxHill = pphs_bh[2, numSeq, ]
epphsHTsims$BoxHill = pphs_bh[3, numSeq, ]
# get expected value (average)
epphs_bh = apply(pphs_bh, c(1,2), mean)

# plot
ggdata00 = data.table(
  x = 1:numSeq, 
  Dlinear = rep(epphs_dopt1[[1]][1], numSeq), 
  Dquadratic = rep(epphs_dopt2[[1]][1], numSeq), 
  SpaceFill = rep(epphs_space[[1]][1], numSeq), 
  SeqMED = epphs_seqmed[ 1, ], 
  Hypothesis = rep("H0", numSeq), 
  key = "x"
)
ggdata01 = data.table(
  x = seq(1, numSeq, length.out = N), 
  BoxHill = epphs_bh[ 1, ], 
  Hypothesis = rep("H0", N), 
  key = "x"
)
ggdata10 = data.table(
  x = 1:numSeq, 
  Dlinear = rep(epphs_dopt1[[1]][2], numSeq), 
  Dquadratic = rep(epphs_dopt2[[1]][2], numSeq), 
  SpaceFill = rep(epphs_space[[1]][2], numSeq), 
  SeqMED = epphs_seqmed[ 2, ],
  Hypothesis = rep("H1", numSeq), 
  key = "x"
)
ggdata11 = data.table(
  x = seq(1, numSeq, length.out = N), 
  BoxHill = epphs_bh[ 2, ], 
  Hypothesis = rep("H1", N), 
  key = c("x")
)
ggdataT0 = data.table(
  x = 1:numSeq, 
  Dlinear = rep(epphs_dopt1[[1]][3], numSeq), 
  Dquadratic = rep(epphs_dopt2[[1]][3], numSeq), 
  SpaceFill = rep(epphs_space[[1]][3], numSeq), 
  SeqMED = epphs_seqmed[ 3, ],
  Hypothesis = rep("HT", numSeq), 
  key = "x"
)
ggdataT1 = data.table(
  x = seq(1, numSeq, length.out = N), 
  BoxHill = epphs_bh[ 3, ], 
  Hypothesis = rep("HT", N), 
  key = c("x")
)
ggdata.woBH = merge(ggdata00, ggdata10, all = TRUE, 
                    by = names(ggdata00))
ggdata.woBH = merge(ggdata.woBH, ggdataT0, all = TRUE, 
                    by = names(ggdata00))
ggdata.BH = merge(ggdata01, ggdata11, all = TRUE, 
                  by = names(ggdata01))
ggdata.BH = merge(ggdata.BH, ggdataT1, all = TRUE, 
                  by = names(ggdata01))
ggdata.tog = merge(ggdata.woBH, ggdata.BH, all = TRUE, 
                   by = c("x", "Hypothesis"))
ggdata.m2 = melt(ggdata.tog, id = c("x", "Hypothesis"), value.name = "epph", 
                 variable.name = "Design")
ggdata0.m2 = melt(ggdata.woBH, id = c("x", "Hypothesis"), value.name = "epph", 
                  variable.name = "Design")
ggdata1.m2 = melt(ggdata.BH, id = c("x", "Hypothesis"), value.name = "epph", 
                  variable.name = "Design")
epph.plt = ggplot(ggdata0.m2, aes(x = x, y = epph, color = Design, linetype = Design)) +
  facet_wrap(vars(Hypothesis)) + 
  geom_path() + 
  scale_linetype_manual(values=c("solid", "dashed", "dashed", "solid", "dashed")) +
  geom_point(data = ggdata0.m2[x == numSeq], aes(x = x, y = epph)) + 
  scale_x_continuous(breaks = 1:10) +
  theme_bw() + 
  theme(panel.grid.minor = element_blank(), axis.text.x = element_text(size = 5)) + 
  labs(y = "", x = "Stages") + 
  ylim(0, 1) + 
  geom_path(data = ggdata1.m2, mapping = aes(x = x, y = epph, color = Design, 
                                             linetype = Design), 
            inherit.aes = FALSE) + 
  geom_point(data = ggdata1.m2[x == numSeq], aes(x = x, y = epph))
epph.plt

# calculate standard errors
dim(epphsH0sims)
dim(epphsH1sims)
dim(epphsHTsims)
H0stats = data.frame(
  mean = apply(epphsH0sims, MARGIN = 2, FUN = mean),
  sd = apply(epphsH0sims, MARGIN = 2, FUN = sd),
  se = apply(epphsH0sims, MARGIN = 2, FUN = sd) / numSims
)
H0stats$lower = H0stats$mean - 2 * H0stats$se
H0stats$upper = H0stats$mean + 2 * H0stats$se
round(H0stats, 5)
H1stats = data.frame(
  mean = apply(epphsH1sims, MARGIN = 2, FUN = mean),
  sd = apply(epphsH1sims, MARGIN = 2, FUN = sd),
  se = apply(epphsH1sims, MARGIN = 2, FUN = sd) / numSims
)
H1stats$lower = H1stats$mean - 2 * H1stats$se
H1stats$upper = H1stats$mean + 2 * H1stats$se
round(H1stats, 5)
HTstats = data.frame(
  mean = apply(epphsHTsims, MARGIN = 2, FUN = mean),
  sd = apply(epphsHTsims, MARGIN = 2, FUN = sd),
  se = apply(epphsHTsims, MARGIN = 2, FUN = sd) / numSims
)
HTstats$lower = HTstats$mean - 2 * HTstats$se
HTstats$upper = HTstats$mean + 2 * HTstats$se
round(HTstats, 5)
