################################################################################
# last updated: 12/15/20
# purpose: to calculate and plot metrics for scenario 1 simulations

################################################################################
# Sources/Libraries
################################################################################
# tamu cluster
# home = "/scratch/user/kristynp/smed_ms/"
# output_home = paste(home,"run_designs/",sep="")
output_home = "run_designs/updated_simulations/lm"
functions_home = "functions"

# for seqmed design
source(paste(functions_home, "/SeqMED.R", sep = ""))
source(paste(functions_home, "/SeqMED_batch.R", sep = ""))
source(paste(functions_home, "/charge_function_q.R", sep = ""))
source(paste(functions_home, "/construct_design_matrix.R", sep = ""))
source(paste(functions_home, "/wasserstein_distance.R", sep = ""))
source(paste(functions_home, "/posterior_parameters.R", sep = ""))
source(paste(functions_home, "/simulate_y.R", sep = ""))

# for evaluating designs
source(paste(functions_home, "/simulate_y.R", sep = ""))
source(paste(functions_home, "/postprob_hypotheses.R", sep = ""))
source(paste(functions_home, "/posterior_mean_mse.R", sep = ""))
# source(paste(functions_home, "/plot_utils.R", sep = ""))
source(paste(functions_home, "/predictive_yhat_mse.R", sep = ""))

# for generating initial data
source(paste(functions_home, "/MMED.R", sep = ""))
source(paste(functions_home, "/variance_marginal_y.R", sep = ""))

# for box-hill deisign
source(paste(functions_home, "/boxhill.R", sep = ""))

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
numSims = 25

# simulation settings
numSeq = 10
seqN = 10
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
seqmeds = readRDS(paste0(output_home, "/scenario1_seqmed_simulations", 
                         "_numSeq", numSeq, 
                         "_seqN", seqN,
                         "_numSims", numSims,
                         ".rds", sep = ""))
bhs = readRDS(paste0(output_home, "/scenario1_boxhill_simulations", 
                     "_N", N, 
                     "_MMEDinput", as.numeric(MMEDinputdata),
                     "_numSims", numSims, 
                     ".rds", sep = ""))
bhs = bhs$bh_list

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
which.sim = 1

# plot a seqmed
seqmed1 = seqmeds[[which.sim]]
ggdata = data.frame(x = seqmed1$D, y = seqmed1$y)
design.seqmed = ggplot(ggdata) + 
  geom_histogram(binwidth = 0.12, closed = "right", aes(x =x, y = after_stat(density))) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
data.seqmed = ggplot(ggdata) + 
  geom_point(aes(x, y), color = gg_color_hue(1)) +
  stat_function(fun = fT) + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggarrange(design.seqmed, data.seqmed)
# ggsave("seqmed_d2.pdf",
#        plot = last_plot(),
#        device = "pdf",
#        path = image_path,
#        scale = 1,
#        width = 4.5,
#        height = 2,
#        units = c("in")
# )

################################################################################
# plot the wasserstein distance and model fits
################################################################################
numseq = 1e2
x_seq = seq(from = xmin, to = xmax, length.out = numseq)

w_seq = sapply(x_seq, function(x) WNlm(
  x, seqmed1$postmean0[, numSeq], seqmed1$postmean1[, numSeq], 
  diag(seqmed1$postvar0[, numSeq]), diag(seqmed1$postvar1[, numSeq]), 
  sigmasq, type01))
f1est = function(x) seqmed1$postmean1[1, seqN] + 
  seqmed1$postmean1[2, seqN] * x + seqmed1$postmean1[3, seqN] * x^2
f2est = function(x) seqmed1$postmean0[1, seqN] + 
  seqmed1$postmean0[2, seqN] * x
f1est_seq = sapply(x_seq, f1est)
f2est_seq = sapply(x_seq, f2est)
fT_seq = sapply(x_seq, fT)

ggdata = data.table::data.table(
  x = x_seq, 
  `Estimated Quadratic` = f1est_seq, 
  `Estimated Line` = f2est_seq, 
  `True Quadratic` = fT_seq, 
  `Wasserstein` = w_seq
)
ggdata = data.table::melt(ggdata, id = c("x"), value.name = "y", variable.name = "Function")

ggdata_ribbon = data.table::data.table(
  x = x_seq, 
  ymin = apply(cbind(f1est_seq, f2est_seq), 1, min), 
  ymax = apply(cbind(f1est_seq, f2est_seq), 1, max)
)
ggplot(ggdata, aes(x = x, y = y, color = Function, linetype = Function)) + 
  scale_linetype_manual(values = c(2, 2, 1, 1)) +
  scale_color_manual(values = c(gg_color_hue(3)[c(1, 3)], "black", gg_color_hue(3)[2])) + 
  geom_path() +
  geom_ribbon(data = ggdata_ribbon, mapping = aes(x = x, ymin = ymin, ymax = ymax), alpha = 0.2, 
              inherit.aes = FALSE) +
  theme_bw() + 
  theme(panel.grid.minor = element_blank())
# ggsave("w_d2.pdf",
#        plot = last_plot(),
#        device = "pdf",
#        path = image_path,
#        scale = 1,
#        width = 4.5,
#        height = 2,
#        units = c("in")
# )

################################################################################
# plot the posterior probabilities of the hypotheses
################################################################################

# first, calculate the posterior probabilities

# for the non-sequential methods;
epphs_space = calcExpPostProbH(
  space_filling, N, betaT, mu0, mu1, V0, V1, sigmasq, numSims = 100, 
  typeT, type01, seed = 123)
epphs_dopt1 = calcExpPostProbH(
  dopt_linear, N, betaT, mu0, mu1, V0, V1, sigmasq, numSims = 100, 
  typeT, type01, seed = 123)
epphs_dopt2 = calcExpPostProbH(
  dopt_quadratic, N, betaT, mu0, mu1, V0, V1, sigmasq, numSims = 100, 
  typeT, type01, seed = 123)

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
# then calculate, like seqmeds -- this is done in batches, even though it's
#   fully sequential, only for convenience!!!!!!!!!!!!!!! ######################
# checking whether the posterior probabilities match,
#   using calcExpPostProbH -- they do!
postprobs0seq = matrix(NA, numSeq, numSims)
postprobs1seq = matrix(NA, numSeq, numSims)
BF01seq = matrix(NA, numSeq, numSims)
for(k in 1:numSims){
  bh_data_k = bhs.format[[k]]
  for(i in 1:numSeq){
    changing_postprobs = calcExpPostProbH_data(
      bh_data_k$y[1:(seqN * i)], bh_data_k$D[1:(seqN * i)], 
      N = seqN * i, mu0, V0, mu1, V1, sigmasq, type01)
    postprobs0seq[i, k] = changing_postprobs[1]
    postprobs1seq[i, k] = changing_postprobs[2]
    BF01seq[i, k] = changing_postprobs[3]
  }
}
# get expected value (average)
epph0seq_bh = apply(postprobs0seq, 1, mean)
epph1seq_bh = apply(postprobs1seq, 1, mean)
BF01seq_bh = apply(BF01seq, 1, mean)

# plot 1
ggdata0 = data.table(
  x = 1:numSeq, 
  Dlinear = rep(epphs_dopt1[1], numSeq), 
  Dquadratic = rep(epphs_dopt2[1], numSeq), 
  SpaceFilling = rep(epphs_space[1], numSeq), 
  SeqMED = epph0seq_seqmed,
  BoxHill = epph0seq_bh, 
  Hypothesis = rep("H0", numSeq)
)
ggdata1 = data.table(
  x = 1:numSeq, 
  Dlinear = rep(epphs_dopt1[2], numSeq), 
  Dquadratic = rep(epphs_dopt2[2], numSeq), 
  SpaceFilling = rep(epphs_space[2], numSeq), 
  SeqMED = epph1seq_seqmed,
  BoxHill = epph1seq_bh, 
  Hypothesis = rep("H1", numSeq)
)
ggdata = rbind(ggdata0, ggdata1)
ggdata.melted = melt(ggdata, id = c("x", "Hypothesis"), value.name = "epph", 
                     variable.name = "Design")
epph.plt = ggplot(ggdata.melted, aes(x = x, y = epph, color = Design, linetype = Design)) +
  facet_wrap(vars(Hypothesis)) + 
  geom_path(size = 1) + 
  scale_linetype_manual(values=c(rep("dashed", 3), rep("solid", 3))) + 
  geom_point(data = ggdata.melted[x == numSeq], aes(x = x, y = epph), size = 2) + 
  theme_bw(base_size = 20) + 
  theme(panel.grid.minor = element_blank()) + 
  labs(y = "", x = "Stages") + 
  ylim(0, 1)
epph.plt

# plot 2
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
# get expected value (average)
epph0seq_bh2 = apply(postprobs0seq, 1, mean)
epph1seq_bh2 = apply(postprobs1seq, 1, mean)
BF01seq_bh2 = apply(BF01seq, 1, mean)
ggdata00 = data.table(
  x = 1:numSeq, 
  Dlinear = rep(epphs_dopt1[1], numSeq), 
  Dquadratic = rep(epphs_dopt2[1], numSeq), 
  SpaceFill = rep(epphs_space[1], numSeq), 
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
  Dlinear = rep(epphs_dopt1[2], numSeq), 
  Dquadratic = rep(epphs_dopt2[2], numSeq), 
  SpaceFill = rep(epphs_space[2], numSeq), 
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
  geom_point(data = ggdata0.m2[x == numSeq], aes(x = x, y = epph), size = 2) + 
  theme_bw() + 
  scale_x_continuous(breaks = 1:10) +
  theme(panel.grid.minor = element_blank(), axis.text.x = element_text(size = 5)) + 
  labs(y = "", x = "Stages") + 
  ylim(0, 1) + 
  geom_path(data = ggdata1.m2, mapping = aes(x = x, y = epph, color = Design, 
                                             linetype = Design), 
            inherit.aes = FALSE) + 
  geom_point(data = ggdata1.m2[x == numSeq], aes(x = x, y = epph), size = 2)
# ggsave("epph_d2.pdf",
#        plot = last_plot(),
#        device = "pdf",
#        path = image_path,
#        scale = 1,
#        width = 4.5,
#        height = 2,
#        units = c("in")
# )

################################################################################
# plot the MSE of beta-hat (posterior mean) of the hypotheses
################################################################################
sim.idx = 1

MSEbetahat_doptlin = getClosedMSE(
  dopt_linear, N, betaT, mu1, V1, sigmasq, typeT)$MSE_postmean
MSEbetahat_doptquad = getClosedMSE(
  dopt_quadratic, N, betaT, mu1, V1, sigmasq, typeT)$MSE_postmean
MSEbetahat_space = getClosedMSE(
  space_filling, N, betaT, mu1, V1, sigmasq, typeT)$MSE_postmean
MSEbetahat_seqmed = getClosedMSE(
  seqmeds[[sim.idx]]$D, N, betaT, mu1, V1, sigmasq, typeT)$MSE_postmean
MSEbetahat_bh = getClosedMSE(
  bhs.format[[sim.idx]]$D, N, betaT, mu1, V1, sigmasq, typeT)$MSE_postmean

b0 = c(MSEbetahat_doptlin[1], MSEbetahat_doptquad[1], MSEbetahat_space[1], 
       MSEbetahat_seqmed[1], MSEbetahat_bh[1])
b1 = c(MSEbetahat_doptlin[2], MSEbetahat_doptquad[2], MSEbetahat_space[2], 
       MSEbetahat_seqmed[2], MSEbetahat_bh[2])
b2 = c(MSEbetahat_doptlin[3], MSEbetahat_doptquad[3], MSEbetahat_space[3], 
       MSEbetahat_seqmed[3], MSEbetahat_bh[3])

ggdata = data.frame(
  Designs = rep(c("Dlinear", "Dquadratic", "SpaceFill", "SeqMED", "BoxHill"), 3), 
  MSE = c(b0, b1, b2), beta = rep(c("B0", "B1", "B2"), each = length(b0)))
mseb.plt = ggplot(ggdata, aes(x = Designs, y = MSE)) + 
  geom_bar(stat = "identity") +
  facet_wrap(vars(beta)) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.text.x = element_text(angle = 45, vjust = 0.5)) +
  labs(y = NULL)
mseb.plt
# ggsave("mseb_d2.pdf",
#        plot = last_plot(),
#        device = "pdf",
#        path = image_path,
#        scale = 1,
#        width = 4.5,
#        height = 2,
#        units = c("in")
# )

################################################################################
# plot the MSE of y-hat (posterior mean) of the hypotheses
################################################################################
sim.idx = 1
x_seq2 = seq(from = -1.25, to = 1.25, length.out = 1e4)

yhatmse_space = getClosedMSEyhat_seq(
  x_seq2, space_filling, N, betaT, typeT, mu1, V1, sigmasq, type01[2])
yhatmse_doptquad = getClosedMSEyhat_seq(
  x_seq2, dopt_quadratic, N, betaT, typeT, mu1, V1, sigmasq, type01[2])
yhatmse_doptlin = getClosedMSEyhat_seq(
  x_seq2, dopt_linear, N, betaT, typeT, mu1, V1, sigmasq, type01[2])
yhatmse_seqmed = getClosedMSEyhat_seq(
  x_seq2, seqmeds[[1]]$D, N, betaT, typeT, mu1, V1, sigmasq, type01[2])
yhatmse_bh = getClosedMSEyhat_seq(
  x_seq2, bhs.format[[1]]$D, N, betaT, typeT, mu1, V1, sigmasq, type01[2])

ylimarg = range(
  0, yhatmse_space$MSEyhat, yhatmse_doptquad$MSEyhat, yhatmse_seqmed$MSEyhat, 
  yhatmse_bh$MSEyhat)

ggdata = data.table(
  x = x_seq2, 
  Dlinear = yhatmse_doptlin$MSEyhat, 
  Dquadratic = yhatmse_doptquad$MSEyhat, 
  SpaceFill = yhatmse_space$MSEyhat, 
  SeqMED = yhatmse_seqmed$MSEyhat,
  BoxHill = yhatmse_bh$MSEyhat
)
ggdata = melt(ggdata, id = c("x"), value.name = "yhatmse", variable.name = "Design")
msey.plt = ggplot(ggdata, aes(x = x, y = yhatmse, color = Design)) +
  coord_cartesian(ylim = ylimarg, xlim = c(-1, 1)) +
  geom_path(size = 1) + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  labs(y = "", x = "x")
msey.plt
# ggsave("msey_d2.pdf",
#        plot = last_plot(),
#        device = "pdf",
#        path = image_path,
#        scale = 1,
#        width = 4.5,
#        height = 2,
#        units = c("in")
# )

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
seqmeds = readRDS(paste0(output_home, "/scenario2_seqmed_simulations", 
                         "_numSeq", numSeq, 
                         "_seqN", seqN,
                         "_numSims", numSims,
                         ".rds", sep = ""))
bhs = readRDS(paste0(output_home, "/scenario2_boxhill_simulations", 
                     "_N", N, 
                     "_MMEDinput", as.numeric(MMEDinputdata),
                     "_numSims", numSims, 
                     ".rds", sep = ""))
bhs = bhs$bh_list

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
which.sim = 2

# plot a seqmed
seqmed1 = seqmeds[[which.sim]]
ggdata = data.frame(x = seqmed1$D, y = seqmed1$y)
design.seqmed = ggplot(ggdata) + 
  geom_histogram(binwidth = 0.12, closed = "right", aes(x =x, y = after_stat(density))) + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
data.seqmed = ggplot(ggdata) + 
  geom_point(aes(x, y), color = gg_color_hue(1)) +
  stat_function(fun = fT) + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggarrange(design.seqmed, data.seqmed)
# ggsave("seqmed_d3.pdf",
#        plot = last_plot(),
#        device = "pdf",
#        path = image_path,
#        scale = 1,
#        width = 4.5,
#        height = 2,
#        units = c("in")
# )

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
fTest = function(x) seqmed1.postmean2[1, seqN] + 
  seqmed1.postmean2[2, seqN] * x + 
  seqmed1.postmean2[3, seqN] * x^2 + 
  seqmed1.postmean2[4, seqN] * x^3
f2est = function(x) seqmed1$postmean1[1, seqN] + 
  seqmed1$postmean1[2, seqN] * x + 
  seqmed1$postmean1[3, seqN] * x^2
f1est = function(x) seqmed1$postmean0[1, seqN] + 
  seqmed1$postmean0[2, seqN] * x

f1est_seq = sapply(x_seq, f1est)
f2est_seq = sapply(x_seq, f2est)
fTest_seq = sapply(x_seq, fTest)
fT_seq = sapply(x_seq, fT)

ggdata_est = data.table::data.table(
  x = x_seq, 
  `Est Line` = f1est_seq, 
  `Est Quadratic` = f2est_seq, 
  `Est Cubic` = fTest_seq
)
ggdata_est = data.table::melt(ggdata_est, id = c("x"), value.name = "y", variable.name = "Function")

ggdata_true = data.table::data.table(x = x_seq, y = fT_seq)

yrange = range(seqmed1$y)
ggplot(ggdata_est) + 
  facet_wrap(facets = vars(Function)) +
  geom_path(aes(x = x, y = y, color = Function)) + 
  geom_path(data = ggdata_true, aes(x, y, color = "True Cubic")) + 
  scale_color_manual(values = c(gg_color_hue(3), "black")) +
  ylim(range(seqmed1$y)) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.text.x = element_text(size = 5))
# ggsave("fits_d3.pdf",
#        plot = last_plot(),
#        device = "pdf",
#        path = image_path,
#        scale = 1,
#        width = 4.5,
#        height = 1.5,
#        units = c("in")
# )

################################################################################
# plot the posterior probabilities of the hypotheses
################################################################################

# # first, calculate the posterior probabilities

models = list("H0" = list(mu0, V0, 2),
              "H1" = list(mu1, V1, 3),
              "H2" = list(betaT, diag(rep(sigmasq01, 4)), 4))

# non-sequential methods
epphs_space = calcEPPH(
  space_filling, N, betaT, typeT, models, sigmasq, numSims = 100, seed = 123)
epphs_dopt1 = calcEPPH(
  dopt_linear,  N, betaT, typeT, models, sigmasq, numSims = 100, seed = 123)
epphs_dopt2 = calcEPPH(
  dopt_quadratic, N, betaT, typeT, models, sigmasq, numSims = 100, seed = 123)

numModels = length(models)

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
pphs_bh = array(NA, dim = c(numModels, numSeq, numSims))
for(k in 1:numSims){
  bh_data_k = bhs.format[[k]]
  for(i in 1:numSeq){
    pphs_bh[ , i, k] = calcEPPHdata(bh_data_k$y[1:(seqN * i)], 
                                    bh_data_k$D[1:(seqN * i)], 
                                    N = seqN * i, models, sigmasq)
  }
}
epphs_bh = apply(pphs_bh, c(1,2), mean)

# plot
pphs_bh = array(NA, dim = c(numModels, N, numSims))
for(k in 1:numSims){
  bh_data_k = bhs.format[[k]]
  for(i in 1:N){
    pphs_bh[ , i, k] = calcEPPHdata(bh_data_k$y[1:i], 
                                    bh_data_k$D[1:i], 
                                    N = i, models, sigmasq)
  }
}
epphs_bh2 = apply(pphs_bh, c(1,2), mean)
ggdata00 = data.table(
  x = 1:numSeq, 
  Dlinear = rep(epphs_dopt1[1], numSeq), 
  Dquadratic = rep(epphs_dopt2[1], numSeq), 
  SpaceFill = rep(epphs_space[1], numSeq), 
  SeqMED = epphs_seqmed[ 1, ], 
  Hypothesis = rep("H0", numSeq), 
  key = "x"
)
ggdata01 = data.table(
  x = seq(1, numSeq, length.out = N), 
  BoxHill = epphs_bh2[ 1, ], 
  Hypothesis = rep("H0", N), 
  key = "x"
)
ggdata10 = data.table(
  x = 1:numSeq, 
  Dlinear = rep(epphs_dopt1[2], numSeq), 
  Dquadratic = rep(epphs_dopt2[2], numSeq), 
  SpaceFill = rep(epphs_space[2], numSeq), 
  SeqMED = epphs_seqmed[ 2, ],
  Hypothesis = rep("H1", numSeq), 
  key = "x"
)
ggdata11 = data.table(
  x = seq(1, numSeq, length.out = N), 
  BoxHill = epphs_bh2[ 2, ], 
  Hypothesis = rep("H1", N), 
  key = c("x")
)
ggdataT0 = data.table(
  x = 1:numSeq, 
  Dlinear = rep(epphs_dopt1[2], numSeq), 
  Dquadratic = rep(epphs_dopt2[2], numSeq), 
  SpaceFill = rep(epphs_space[2], numSeq), 
  SeqMED = epphs_seqmed[ 3, ],
  Hypothesis = rep("HT", numSeq), 
  key = "x"
)
ggdataT1 = data.table(
  x = seq(1, numSeq, length.out = N), 
  BoxHill = epphs_bh2[ 3, ], 
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
  geom_point(data = ggdata0.m2[x == numSeq], aes(x = x, y = epph), size = 2) + 
  scale_x_continuous(breaks = 1:10) +
  theme_bw() + 
  theme(panel.grid.minor = element_blank(), axis.text.x = element_text(size = 5)) + 
  labs(y = "", x = "Stages") + 
  ylim(0, 1) + 
  geom_path(data = ggdata1.m2, mapping = aes(x = x, y = epph, color = Design, 
                                             linetype = Design), 
            inherit.aes = FALSE) + 
  geom_point(data = ggdata1.m2[x == numSeq], aes(x = x, y = epph), size = 2)
epph.plt
# ggsave("epph_d3.pdf",
#        plot = last_plot(),
#        device = "pdf",
#        path = image_path,
#        scale = 1,
#        width = 4.5,
#        height = 2,
#        units = c("in")
# )

################################################################################
# plot the MSE of beta-hat (posterior mean) of the hypotheses
################################################################################
# given H2
sim.idx = 1

# define new priors
mu2 = rep(0, 4)
V2 = diag(rep(sigmasq01, length(mu2)))

MSEbetahat_doptlin = getClosedMSE(
  dopt_linear, N, betaT, mu2, V2, sigmasq, typeT)$MSE_postmean
MSEbetahat_doptquad = getClosedMSE(
  dopt_quadratic, N, betaT, mu2, V2, sigmasq, typeT)$MSE_postmean
MSEbetahat_space = getClosedMSE(
  space_filling, N, betaT, mu2, V2, sigmasq, typeT)$MSE_postmean
MSEbetahat_seqmed = getClosedMSE(
  seqmeds[[sim.idx]]$D, N, betaT, mu2, V2, sigmasq, typeT)$MSE_postmean
MSEbetahat_bh = getClosedMSE(
  bhs.format[[sim.idx]]$D, N, betaT, mu2, V2, sigmasq, typeT)$MSE_postmean

b0 = c(MSEbetahat_doptlin[1], MSEbetahat_doptquad[1], MSEbetahat_space[1], 
       MSEbetahat_seqmed[1], MSEbetahat_bh[1])
b1 = c(MSEbetahat_doptlin[2], MSEbetahat_doptquad[2], MSEbetahat_space[2], 
       MSEbetahat_seqmed[2], MSEbetahat_bh[2])
b2 = c(MSEbetahat_doptlin[3], MSEbetahat_doptquad[3], MSEbetahat_space[3], 
       MSEbetahat_seqmed[3], MSEbetahat_bh[3])
b3 = c(MSEbetahat_doptlin[4], MSEbetahat_doptquad[4], MSEbetahat_space[4], 
       MSEbetahat_seqmed[4], MSEbetahat_bh[4])

ggdata = data.frame(
  Designs = rep(c("Dlinear", "Dquadratic", "SpaceFill", "SeqMED", "BoxHill"), 4), 
  MSE = c(b0, b1, b2, b3), beta = rep(c("B0", "B1", "B2", "B3"), each = length(b0)))
mseb.plt = ggplot(ggdata, aes(x = Designs, y = MSE)) + 
  geom_bar(stat = "identity") +
  facet_wrap(vars(beta)) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.text.x = element_text(angle = 45, vjust = 0.5)) +
  labs(y = NULL)
mseb.plt
# ggsave("mseb_d3.pdf",
#        plot = last_plot(),
#        device = "pdf",
#        path = image_path,
#        scale = 1,
#        width = 4.5,
#        height = 3,
#        units = c("in")
# )

################################################################################
# plot the MSE of y-hat (posterior mean) of the hypotheses
################################################################################

# given H2
sim.idx = 1
x_seq2 = seq(from = -1.25, to = 1.25, length.out = 1e4)

# define new priors
mu2 = rep(0, 4)
V2 = diag(rep(sigmasq01, length(mu2)))

yhatmse_space = getClosedMSEyhat_seq(
  x_seq2, space_filling, N, betaT, typeT, mu2, V2, sigmasq, typeT)
yhatmse_doptquad = getClosedMSEyhat_seq(
  x_seq2, dopt_quadratic, N, betaT, typeT, mu2, V2, sigmasq, typeT)
yhatmse_doptlin = getClosedMSEyhat_seq(
  x_seq2, dopt_linear, N, betaT, typeT, mu2, V2, sigmasq, typeT)
yhatmse_seqmed = getClosedMSEyhat_seq(
  x_seq2, seqmeds[[sim.idx]]$D, N, betaT, typeT, mu2, V2, sigmasq, typeT)
yhatmse_bh = getClosedMSEyhat_seq(
  x_seq2, bhs.format[[sim.idx]]$D, N, betaT, typeT, mu2, V2, sigmasq, typeT)

ylimarg = range(
  0, yhatmse_space$MSEyhat, yhatmse_doptquad$MSEyhat, yhatmse_seqmed$MSEyhat, 
  yhatmse_bh$MSEyhat)
ylimarg = c(0, 0.12)

ggdata = data.table(
  x = x_seq2, 
  Dlinear = yhatmse_doptlin$MSEyhat, 
  Dquadratic = yhatmse_doptquad$MSEyhat, 
  SpaceFill = yhatmse_space$MSEyhat, 
  SeqMED = yhatmse_seqmed$MSEyhat, 
  BoxHill = yhatmse_bh$MSEyhat
)
ggdata = melt(ggdata, id = c("x"), value.name = "yhatmse", variable.name = "Design")
msey.plt = ggplot(ggdata, aes(x = x, y = yhatmse, color = Design)) +
  coord_cartesian(ylim = ylimarg, xlim = c(-1, 1)) +
  geom_path() + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  labs(y = "", x = "x")
msey.plt
# ggsave("msey_d3.pdf",
#        plot = last_plot(),
#        device = "pdf",
#        path = image_path,
#        scale = 1,
#        width = 4.5,
#        height = 2,
#        units = c("in")
# )




