################################################################################
# last updated: 03/17/21
# purpose: to calculate and plot metrics for gp simulations



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
# directories
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
source(paste(functions_home, "/predictive_yhat_mse.R", sep = ""))

library(Matrix)
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
numSims = 100
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
f1est = function(x) seqmed1$postmean1[1, numSeq] +
  seqmed1$postmean1[2, numSeq] * x + seqmed1$postmean1[3, numSeq] * x^2
f2est = function(x) seqmed1$postmean0[1, numSeq] +
  seqmed1$postmean0[2, numSeq] * x
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
  space_filling, N, betaT, mu0, mu1, V0, V1, sigmasq, numSims, 
  typeT, type01, seed = 123)
epphs_dopt1 = calcExpPostProbH(
  dopt_linear, N, betaT, mu0, mu1, V0, V1, sigmasq, numSims, 
  typeT, type01, seed = 123)
epphs_dopt2 = calcExpPostProbH(
  dopt_quadratic, N, betaT, mu0, mu1, V0, V1, sigmasq, numSims, 
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
  geom_point(data = ggdata0.m2[x == numSeq], aes(x = x, y = epph)) + 
  theme_bw() + 
  scale_x_continuous(breaks = c(0, 25, 50, 75, 100)) +
  theme(panel.grid.minor = element_blank(), axis.text.x = element_text(size = 5)) + 
  labs(y = "", x = "Stages") + 
  ylim(0, 1) + 
  geom_path(data = ggdata1.m2, mapping = aes(x = x, y = epph, color = Design, 
                                             linetype = Design), 
            inherit.aes = FALSE) + 
  geom_point(data = ggdata1.m2[x == numSeq], aes(x = x, y = epph))
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

MSEbetahat_doptlin = getMSEBeta(
  dopt_linear, N, betaT, mu1, V1, sigmasq, typeT)$MSE_postmean
MSEbetahat_doptquad = getMSEBeta(
  dopt_quadratic, N, betaT, mu1, V1, sigmasq, typeT)$MSE_postmean
MSEbetahat_space = getMSEBeta(
  space_filling, N, betaT, mu1, V1, sigmasq, typeT)$MSE_postmean
MSEbetahat_seqmed = getMSEBeta(
  seqmeds[[which.sim]]$D, N, betaT, mu1, V1, sigmasq, typeT)$MSE_postmean
MSEbetahat_bh = getMSEBeta(
  bhs.format[[which.sim]]$D, N, betaT, mu1, V1, sigmasq, typeT)$MSE_postmean

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
x_seq2 = seq(from = -1.25, to = 1.25, length.out = 1e4)

getMSEYhat.d2 = function(x, design){
  getMSEYhat(x, design, N, betaT, typeT, mu1, V1, sigmasq, type01[2])$MSEyhat
}

yhatmse_space = sapply(x_seq2, getMSEYhat.d2, design = space_filling)
yhatmse_doptquad = sapply(x_seq2, getMSEYhat.d2, design = dopt_quadratic)
yhatmse_doptlin = sapply(x_seq2, getMSEYhat.d2, design = dopt_linear)
yhatmse_seqmed = sapply(x_seq2, getMSEYhat.d2, design = seqmeds[[which.sim]]$D)
yhatmse_bh = sapply(x_seq2, getMSEYhat.d2, design = bhs.format[[which.sim]]$D)

ylimarg = range(0, yhatmse_space, yhatmse_doptquad, yhatmse_seqmed, yhatmse_bh)

ggdata = data.table(
  x = x_seq2, 
  Dlinear = yhatmse_doptlin, 
  Dquadratic = yhatmse_doptquad, 
  SpaceFill = yhatmse_space, 
  SeqMED = yhatmse_seqmed,
  BoxHill = yhatmse_bh
)
ggdata = melt(ggdata, id = c("x"), value.name = "yhatmse", variable.name = "Design")
msey.plt = ggplot(ggdata, aes(x = x, y = yhatmse, color = Design)) +
  coord_cartesian(ylim = ylimarg, xlim = c(-1, 1)) +
  geom_path() + 
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
  space_filling, N, betaT, typeT, models, sigmasq, numSims, seed = 123)
epphs_dopt1 = calcEPPH(
  dopt_linear,  N, betaT, typeT, models, sigmasq, numSims, seed = 123)
epphs_dopt2 = calcEPPH(
  dopt_quadratic, N, betaT, typeT, models, sigmasq, numSims, seed = 123)

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
pphs_bh = array(NA, dim = c(numModels, N, numSims))
for(k in 1:numSims){
  bh_data_k = bhs.format[[k]]
  for(i in 1:N){
    pphs_bh[ , i, k] = calcEPPHdata(bh_data_k$y[1:i], 
                                    bh_data_k$D[1:i], 
                                    N = i, models, sigmasq)
  }
}
epphs_bh = apply(pphs_bh, c(1,2), mean)

# plot
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
  BoxHill = epphs_bh[ 1, ], 
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
  BoxHill = epphs_bh[ 2, ], 
  Hypothesis = rep("H1", N), 
  key = c("x")
)
ggdataT0 = data.table(
  x = 1:numSeq, 
  Dlinear = rep(epphs_dopt1[3], numSeq), 
  Dquadratic = rep(epphs_dopt2[3], numSeq), 
  SpaceFill = rep(epphs_space[3], numSeq), 
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
  scale_x_continuous(breaks = c(0, 25, 50, 75, 100)) +
  theme_bw() + 
  theme(panel.grid.minor = element_blank(), axis.text.x = element_text(size = 5)) + 
  labs(y = "", x = "Stages") + 
  ylim(0, 1) + 
  geom_path(data = ggdata1.m2, mapping = aes(x = x, y = epph, color = Design, 
                                             linetype = Design), 
            inherit.aes = FALSE) + 
  geom_point(data = ggdata1.m2[x == numSeq], aes(x = x, y = epph))
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

# define new priors
mu2 = rep(0, 4)
V2 = diag(rep(sigmasq01, length(mu2)))

MSEbetahat_doptlin = getMSEBeta(
  dopt_linear, N, betaT, mu2, V2, sigmasq, typeT)$MSE_postmean
MSEbetahat_doptquad = getMSEBeta(
  dopt_quadratic, N, betaT, mu2, V2, sigmasq, typeT)$MSE_postmean
MSEbetahat_space = getMSEBeta(
  space_filling, N, betaT, mu2, V2, sigmasq, typeT)$MSE_postmean
MSEbetahat_seqmed = getMSEBeta(
  seqmeds[[which.sim]]$D, N, betaT, mu2, V2, sigmasq, typeT)$MSE_postmean
MSEbetahat_bh = getMSEBeta(
  bhs.format[[which.sim]]$D, N, betaT, mu2, V2, sigmasq, typeT)$MSE_postmean

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
x_seq2 = seq(from = -1.25, to = 1.25, length.out = 1e4)

# define new priors
mu2 = rep(0, 4)
V2 = diag(rep(sigmasq01, length(mu2)))

getMSEYhat.d3 = function(x, design){
  getMSEYhat(x, design, N, betaT, typeT, mu2, V2, sigmasq, typeT)$MSEyhat
}

yhatmse_space = sapply(x_seq2, getMSEYhat.d3, design = space_filling)
yhatmse_doptquad = sapply(x_seq2, getMSEYhat.d3, design = dopt_quadratic)
yhatmse_doptlin = sapply(x_seq2, getMSEYhat.d3, design = dopt_linear)
yhatmse_seqmed = sapply(x_seq2, getMSEYhat.d3, design = seqmeds[[which.sim]]$D)
yhatmse_bh = sapply(x_seq2, getMSEYhat.d3, design = bhs.format[[which.sim]]$D)

# yhatmse_space = getMSEYhat_seq(
#   x_seq2, space_filling, N, betaT, typeT, mu2, V2, sigmasq, typeT)
# yhatmse_doptquad = getMSEYhat_seq(
#   x_seq2, dopt_quadratic, N, betaT, typeT, mu2, V2, sigmasq, typeT)
# yhatmse_doptlin = getMSEYhat_seq(
#   x_seq2, dopt_linear, N, betaT, typeT, mu2, V2, sigmasq, typeT)
# yhatmse_seqmed = getMSEYhat_seq(
#   x_seq2, seqmeds[[which.sim]]$D, N, betaT, typeT, mu2, V2, sigmasq, typeT)
# yhatmse_bh = getMSEYhat_seq(
#   x_seq2, bhs.format[[which.sim]]$D, N, betaT, typeT, mu2, V2, sigmasq, typeT)

ylimarg = range(0, yhatmse_space, yhatmse_doptquad, yhatmse_seqmed, yhatmse_bh)
ylimarg = c(0, 0.12)

ggdata = data.table(
  x = x_seq2, 
  Dlinear = yhatmse_doptlin, 
  Dquadratic = yhatmse_doptquad, 
  SpaceFill = yhatmse_space, 
  SeqMED = yhatmse_seqmed, 
  BoxHill = yhatmse_bh
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
###########
############
#############
##############
###############
###############
##############
#############
############
###########
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
################################################################################
################################################################################
################################################################################
################################################################################
# GP Plots
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################

################################################################################
# Sources/Libraries
################################################################################

# directories
output_home = "run_designs/updated_simulations/gp"
functions_home = "functions"

# for seqmed design
source(paste(functions_home, "/SeqMEDgp.R", sep = ""))
source(paste(functions_home, "/SeqMEDgp_batch.R", sep = ""))
source(paste(functions_home, "/charge_function_q.R", sep = ""))
source(paste(functions_home, "/wasserstein_distance.R", sep = ""))
source(paste(functions_home, "/charge_function_q.R", sep = ""))
source(paste(functions_home, "/covariance_functions.R", sep = ""))
source(paste(functions_home, "/gp_predictive.R", sep = ""))

# for box-hill design
source(paste(functions_home, "/boxhill.R", sep = ""))
source(paste(functions_home, "/boxhill_gp.R", sep = ""))
source(paste(functions_home, "/kl_divergence.R", sep = ""))

library(mvtnorm)

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

# helper functions

getRSS01 = function(
  x, 
  x.idx,
  x.new, 
  x.new.idx,
  y.new = NULL,
  function.values,
  type, 
  l, 
  nugget
){
  y = function.values[x.idx]
  if(is.null(y.new)) y.new = function.values[x.new.idx]
  pred0.tmp = getGPPredictive(x.new, x, y, type[1], l[1], nugget)
  pred1.tmp = getGPPredictive(x.new, x, y, type[2], l[2], nugget)
  RSS0.tmp = sum((pred0.tmp$pred_mean - y.new)^2)  
  RSS1.tmp = sum((pred1.tmp$pred_mean - y.new)^2)
  RSS01 = RSS0.tmp / RSS1.tmp
  RSS01
}

################################################################################
# simulation settings, shared for both scenarios
################################################################################

# simulations settings
numSims = 100
seed = 12
Nin = 6
numSeq = 15
seqN = 1
Nnew = numSeq * seqN
Nttl = Nin + Nnew
xmin = 0
xmax = 1
numx = 10^3 + 1
x_seq = seq(from = xmin, to = xmax, length.out = numx)

# SeqMED settings
nuggetSM = 1e-10

# boxhill settings
prior_probs = rep(1 / 2, 2)
nuggetBH = 1e-10

################################################################################
# input data
################################################################################

# 1. make space-filling design
# space_filling = seq(from = xmin, to = xmax, length.out = Nttl)
space_filling_idx = c(1, 1 + ((numx - 1)/(Nttl - 1)) * 1:((numx - 1) / ((numx - 1)/(Nttl - 1))))
space_filling = x_seq[space_filling_idx]

# input set 1 (extrapolation)
x_in1_idx = space_filling_idx[1:Nin]
x_in1 = x_seq[x_in1_idx]
x_spacefill1_idx = space_filling_idx[-c(1:Nin)]
x_spacefill1 = x_seq[x_spacefill1_idx]
# all.equal(space_filling, c(x_in1, x_spacefill1))

# input set 2 (increasing spread)
x_in2_idx = space_filling_idx[c(1, 2, 4, 7, 12, 21)]
x_in2 = x_seq[x_in2_idx]
x_spacefill2_idx = space_filling_idx[-c(1, 2, 4, 7, 12, 21)]
x_spacefill2 = x_seq[x_spacefill2_idx]
# all.equal(space_filling, sort(c(x_in2, x_spacefill2)))

# input set 3 (space-filling / even coverage)
x_in3_idx = c(1, 1 + ((numx - 1)/(Nin - 1)) * 
                1:((numx - 1) / ((numx - 1)/(Nin - 1))))
x_in3 = x_seq[x_in3_idx]
x_spacefill3_idx = space_filling_idx[!(space_filling_idx %in% x_in3_idx)]
x_spacefill3 = x_seq[x_spacefill3_idx]
# all.equal(space_filling, sort(c(x_in3, x_spacefill3)))

# input set 4 (uniform / random)
x_spacefill4_idx = floor(c(1, 1 + ((numx - 1)/(Nnew - 1)) * 
                             1:((numx - 1) / ((numx - 1)/(Nnew - 1)))))
x_spacefill4 = x_seq[x_spacefill4_idx]

# make sets
input.set = list(
  input1 = list(x.new = x_in1, x.new.idx = x_in1_idx), 
  input2 = list(x.new = x_in2, x.new.idx = x_in2_idx), 
  input3 = list(x.new = x_in3, x.new.idx = x_in3_idx), 
  input4 = list(x.new = NULL, x.new.idx = NULL)
)
spacefill.set = list(
  input1 = list(x.new = x_spacefill1, x.new.idx = x_spacefill1_idx), 
  input2 = list(x.new = x_spacefill2, x.new.idx = x_spacefill2_idx), 
  input3 = list(x.new = x_spacefill3, x.new.idx = x_spacefill3_idx), 
  input4 = list(x.new = x_spacefill4, x.new.idx = x_spacefill4_idx)
)

################################################################################
# Scenario 1: Squared exponential vs. matern, true = matern
################################################################################
type01 = c("squaredexponential", "matern")
l01= c(0.01, 0.01)
# generate matern functions
set.seed(seed)
null_cov = getCov(x_seq, x_seq, type01[2], l01[2])
null_mean = rep(0, numx)

# bh settings
model0 = list(type = type01[1], l = l01[1])
model1 = list(type = type01[2], l = l01[2])

################################################################################
# Plot initial data
################################################################################
# different input points

# plot
ggdata = data.table(
  `Extrapolation` = x_in1, 
  `Inc Spread` = x_in2, 
  `Even Coverage` = x_in3, 
  `Random` = runif(Nin)
)
ggdata = melt(ggdata, measure.vars = 1:4)
ggdata$variable = factor(ggdata$variable)
ggdata$Type = "Input"
ggdata.domain = data.frame(
  variable = rep(unique(ggdata$variable), each = 2), 
  value = rep(c(0, 1), 2)
)
ggdata.domain$variable = factor(ggdata.domain$variable)
plt0 = ggplot() + 
  geom_path(data = ggdata.domain, aes(x = value, y = variable), 
            linetype = 2, color = "gray", inherit.aes = FALSE) + 
  geom_point(data = ggdata, aes(x = value, y = variable, col = Type)) +
  scale_color_manual(values = gg_color_hue(2)[1]) + 
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  labs(y = "", x = "x")
plt0
# ggsave("initial.pdf",
#        plot = last_plot(),
#        device = "pdf",
#        path = image_path,
#        scale = 1,
#        width = 4.5,
#        height = 2,
#        units = c("in")
# )

# what do they look like with the spacefilling points?

ggdata2 = data.table(
  `Extrapolation` = x_spacefill1, 
  `Inc Spread` = x_spacefill2, 
  `Even Coverage` = x_spacefill3, 
  `Random` = x_spacefill4
)
ggdata2 = melt(ggdata2, measure.vars = 1:4)
ggdata2$variable = factor(ggdata2$variable)
ggdata$Type = "Input       "
ggdata2$Type = "SpaceFill"
ggdata.ttl = rbind(ggdata, ggdata2)
plt0.1 = ggplot() + 
  geom_path(data = ggdata.domain, aes(x = value, y = variable), 
            linetype = 2, color = "gray", inherit.aes = FALSE) + 
  geom_point(data = ggdata, aes(x = value, y = variable, col = Type)) +
  scale_color_manual(values = gg_color_hue(2)[1]) + 
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  labs(y = "", x = "x")
plt0.2 = ggplot() + 
  geom_path(data = ggdata.domain, aes(x = value, y = variable), 
            linetype = 2, color = "gray", inherit.aes = FALSE) + 
  geom_point(data = ggdata.ttl, aes(x = value, y = variable, col = Type)) +
  scale_color_manual(values = gg_color_hue(4)[c(1,4)]) + 
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  labs(y = "", x = "x")
ggarrange(plt0.1, plt0.2, nrow = 2)
# ggsave("initial_and_spacefill2.pdf",
#        plot = last_plot(),
#        device = "pdf",
#        path = image_path,
#        scale = 1,
#        width = 4.5,
#        height = 4,
#        units = c("in")
# )
################################################################################
# Plot the demo design
################################################################################

# demo settings
l01.demo = c(0.1, 0.1)
nugget.demo = NULL

# read in sims
gaussianvsmatern_train4sims = readRDS(
  paste0(
    "run_designs/gp",
    "/gp_demo/gaussianvsmatern_train4sims.rds"))

# get sim info
sim_ind = 1
x_in = gaussianvsmatern_train4sims$x_train[ , sim_ind]
x_in_idx = gaussianvsmatern_train4sims$x_train_ind[ , sim_ind]
mmed_gp = gaussianvsmatern_train4sims$mmed_gp_list[[sim_ind]]
y_seq = gaussianvsmatern_train4sims$sim_fns[ , sim_ind]
y_in = y_seq[x_in_idx]

newpts = mmed_gp$addD
truey = y_seq[mmed_gp$indices]

H0_predfn = getGPPredictive(x_seq, x_in, y_in, type01[1], l01.demo[1],
                            nugget = nugget.demo)
H1_predfn = getGPPredictive(x_seq, x_in, y_in, type01[2], l01.demo[2],
                            nugget = nugget.demo)

# get w_seq
Kinv0 = solve(getCov(x_in, x_in, type01[1], l01.demo[1]))
Kinv1 = solve(getCov(x_in, x_in, type01[2], l01.demo[2]))
w_seq = sapply(x_seq, FUN = function(x1) 
  WNgp(x1, Kinv0, Kinv1, x_in, y_in, var_e = 1, type01, l01.demo))
# plot
err0 = 2 * sqrt(diag(H0_predfn$pred_var))
err1 = 2 * sqrt(diag(H1_predfn$pred_var))
ggdata = data.table(
  x = x_seq, 
  `True Function` = y_seq, 
  Wasserstein = w_seq, 
  `H0 Predictive` = H0_predfn$pred_mean, 
  `H1 Predictive` = H1_predfn$pred_mean,
  lower0 = H0_predfn$pred_mean - err0, 
  lower1 = H1_predfn$pred_mean - err1, 
  upper0 = H0_predfn$pred_mean + err0,
  upper1 = H1_predfn$pred_mean + err1
)
yrange = range(ggdata$lower0, ggdata$lower1, 
               ggdata$upper0, ggdata$upper1)
yrange[1] = yrange[1] - 1
ggdata$Wasserstein = ggdata$Wasserstein - abs(yrange[1])
ggdata$zero1 = NA
ggdata$zero2 = NA
ggdata.melted = melt(ggdata, id.vars = c("x"), 
                     measure.vars = c("True Function", "Wasserstein", "H0 Predictive", "H1 Predictive"))
ggdata.lower = melt(ggdata, id.vars = c("x"), 
                    measure.vars = c("zero1", "zero2", "lower0", "lower1"))
ggdata.upper = melt(ggdata, id.vars = c("x"), 
                    measure.vars = c("zero1", "zero2", "upper0", "upper1"))
ggdata.melted = cbind(ggdata.melted, 
                      lower = ggdata.lower$value, 
                      upper = ggdata.upper$value)
ggdata_pts = data.table(
  x = c(x_in, newpts), 
  y = c(y_in, truey), 
  color = c(rep(gg_color_hue(2)[2], length(x_in)), 
            rep(gg_color_hue(2)[1], length(newpts))), 
  shape = c(rep(8, length(x_in)), 
            rep(16, length(newpts)))
)
ggplot(data = ggdata.melted, aes(x = x, y =value, color = variable), 
       linetype = 1) + 
  geom_path() + 
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = variable), 
              alpha = 0.1, linetype = 0) +
  scale_linetype_manual(values = c(1, 1, 2, 2)) + 
  scale_fill_manual(values = c(NA, NA, "#00BFC4", "#C77CFF")) + 
  scale_color_manual(values = c(1, "gray", "#00BFC4", "#C77CFF")) + 
  geom_point(data = ggdata_pts, mapping = aes(x = x, y = y), 
             inherit.aes = FALSE, color = ggdata_pts$color, 
             shape = ggdata_pts$shape,
             size = 2) +
  geom_point(data = ggdata_pts, mapping = aes(x = x, y = yrange[1]), 
             inherit.aes = FALSE, color = ggdata_pts$color, 
             shape = ggdata_pts$shape, 
             size = 2) +
  scale_y_continuous(limits = yrange) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  labs(y = "y", x = "x", fill = "Function", color = "Function")
# ggsave("gvm.pdf",
#        plot = last_plot(),
#        device = "pdf",
#        path = image_path,
#        scale = 1,
#        width = 4.5,
#        height = 2.5,
#        units = c("in")
# )

################################################################################
# Metrics
################################################################################

seqmeds = list()
for(i in 1:4){
  seqmeds[[i]] = readRDS(
    paste0(
      output_home, 
      "/seqmed/scenario1_seqmed_simulations", 
      "_input", i, 
      "_Nin6_Nnew15_numSims100.rds"))
}

boxhills = list()
for(i in 1:4){
  boxhills[[i]] = readRDS(
    paste0(
      output_home, 
      "/boxhill/scenario1_boxhill_simulations", 
      "_input", i, 
      "_Nin6_Nnew15_numSims100.rds"))
}

# RSS Ratio (0/1)

# see helper function getRSS01

# # check that seqmeds and boxhills use the same values # # # # # # # # # # # #
# seqmeds.function.values = matrix(NA, nrow = numx * 4, ncol = numSims)
# boxhills.function.values = matrix(NA, nrow = numx * 4, ncol = numSims)
# for(i in 1:numSims){
#   seqmeds.function.values[, i] = c(
#     seqmeds[[1]]$function.values.list[ , i],
#     seqmeds[[2]]$function.values.list[ , i],
#     seqmeds[[3]]$function.values.list[ , i],
#     seqmeds[[4]]$function.values.list[ , i]
#   )
#   boxhills.function.values[, i] = c(
#     boxhills[[1]]$function.values.list[ , i],
#     boxhills[[2]]$function.values.list[ , i],
#     boxhills[[3]]$function.values.list[ , i],
#     boxhills[[4]]$function.values.list[ , i]
#   )
# }
# all.equal(seqmeds.function.values, boxhills.function.values) # TRUE!
# 
# # check that each of seqmeds function.values.lists are equal # # # # # # # # #
# seqmed.function.values1 = matrix(NA, nrow = numx, ncol = numSims)
# seqmed.function.values2 = matrix(NA, nrow = numx, ncol = numSims)
# seqmed.function.values3 = matrix(NA, nrow = numx, ncol = numSims)
# seqmed.function.values4 = matrix(NA, nrow = numx, ncol = numSims)
# for(i in 1:numSims){
#   seqmed.function.values1[, i] = seqmeds[[1]]$function.values.list[ , i]
#   seqmed.function.values2[, i] = seqmeds[[2]]$function.values.list[ , i]
#   seqmed.function.values3[, i] = seqmeds[[3]]$function.values.list[ , i]
#   seqmed.function.values4[, i] = seqmeds[[4]]$function.values.list[ , i]
# }
# all.equal(seqmed.function.values1, seqmed.function.values2) # TRUE!
# all.equal(seqmed.function.values1, seqmed.function.values3) # TRUE!
# all.equal(seqmed.function.values1, seqmed.function.values4) # TRUE!
f.vals.mat = seqmeds[[1]]$function.values.list

# # check that random designs are the same in the seqmed and boxhill sims # # # 
# seqmeds.random.mat = matrix(NA, nrow = Nin, ncol = numSims)
# boxhills.random.mat = matrix(NA, nrow = Nin, ncol = numSims)
# seqmeds.random.idx.mat = matrix(NA, nrow = Nin, ncol = numSims)
# boxhills.random.idx.mat = matrix(NA, nrow = Nin, ncol = numSims)
# for(i in 1:numSims){
#   seqmeds.random.mat[, i] = seqmeds[[4]]$design.list[[i]]$x
#   seqmeds.random.idx.mat[, i] = seqmeds[[4]]$design.list[[i]]$x.idx
#   boxhills.random.mat[, i] = seqmeds[[4]]$design.list[[i]]$x
#   boxhills.random.idx.mat[, i] = seqmeds[[4]]$design.list[[i]]$x.idx
# }
# all.equal(seqmeds.random.mat, boxhills.random.mat) # TRUE!
# all.equal(seqmeds.random.idx.mat, boxhills.random.idx.mat) # TRUE!
random.mat = matrix(NA, nrow = Nin, ncol = numSims)
random.idx.mat = matrix(NA, nrow = Nin, ncol = numSims)
for(i in 1:numSims){
  random.mat[, i] = seqmeds[[4]]$design.list[[i]]$x
  random.idx.mat[, i] = seqmeds[[4]]$design.list[[i]]$x.idx
}

# RSS01 for SeqMED
RSS01_in1_seqmed = rep(NA, numSims)
RSS01_in2_seqmed = rep(NA, numSims)
RSS01_in3_seqmed = rep(NA, numSims)
RSS01_in4_seqmed = rep(NA, numSims)

for(i in 1:numSims){
  RSS01_in1_seqmed[i] = getRSS01(
    x = x_in1, x.idx = x_in1_idx, 
    x.new = seqmeds[[1]]$design.list[[i]]$x.new, 
    x.new.idx = seqmeds[[1]]$design.list[[i]]$x.new.idx, 
    y.new = seqmeds[[1]]$design.list[[i]]$y.new, 
    function.values = f.vals.mat[ , i], type = type01, l = l01, nugget = nuggetSM)
  # plot(x_seq, seqmeds[[1]]$function.values.list[ , i], type = "l")
  # points(seqmeds[[1]]$x0, seqmeds[[1]]$design.list[[i]]$y[1:Nin], col = 3)
  # points(seqmeds[[1]]$design.list[[i]]$D[-c(1:Nin)], 
  #        seqmeds[[1]]$design.list[[i]]$y[-c(1:Nin)], col = 2)
  RSS01_in2_seqmed[i] = getRSS01(
    x = x_in2, x.idx = x_in2_idx, 
    x.new = seqmeds[[2]]$design.list[[i]]$x.new, 
    x.new.idx = seqmeds[[2]]$design.list[[i]]$x.new.idx, 
    y.new = seqmeds[[2]]$design.list[[i]]$y.new, 
    function.values = f.vals.mat[ , i], type = type01, l = l01, nugget = nuggetSM)
  RSS01_in3_seqmed[i] = getRSS01(
    x = x_in3, x.idx = x_in3_idx, 
    x.new = seqmeds[[3]]$design.list[[i]]$x.new, 
    x.new.idx = seqmeds[[3]]$design.list[[i]]$x.new.idx, 
    y.new = seqmeds[[3]]$design.list[[i]]$y.new, 
    function.values = f.vals.mat[ , i], type = type01, l = l01, nugget = nuggetSM)
  RSS01_in4_seqmed[i] = getRSS01(
    x = random.mat[, i], x.idx = random.idx.mat[, i], 
    x.new = seqmeds[[4]]$design.list[[i]]$x.new, 
    x.new.idx = seqmeds[[4]]$design.list[[i]]$x.new.idx, 
    y.new = seqmeds[[4]]$design.list[[i]]$y.new, 
    function.values = f.vals.mat[ , i], type = type01, l = l01, nugget = nuggetSM)
}

# RSS01 for Box-Hill
RSS01_in1_boxhill = rep(NA, numSims)
RSS01_in2_boxhill = rep(NA, numSims)
RSS01_in3_boxhill = rep(NA, numSims)
RSS01_in4_boxhill = rep(NA, numSims)
for(i in 1:numSims){
  RSS01_in1_boxhill[i] = getRSS01(
    x = x_in1, x.idx = x_in1_idx, 
    x.new = as.vector(na.omit(boxhills[[1]]$design.list[[i]]$x.new)), 
    x.new.idx = as.vector(na.omit(boxhills[[1]]$design.list[[i]]$x.new.idx)), 
    y.new = as.vector(na.omit(boxhills[[1]]$design.list[[i]]$y.new)), 
    function.values = f.vals.mat[ , i], type = type01, l = l01, nugget = nuggetSM)
  RSS01_in2_boxhill[i] = getRSS01(
    x = x_in2, x.idx = x_in2_idx, 
    x.new = as.vector(na.omit(boxhills[[2]]$design.list[[i]]$x.new)), 
    x.new.idx = as.vector(na.omit(boxhills[[2]]$design.list[[i]]$x.new.idx)), 
    y.new = as.vector(na.omit(boxhills[[2]]$design.list[[i]]$y.new)), 
    function.values = f.vals.mat[ , i], type = type01, l = l01, nugget = nuggetSM)
  RSS01_in3_boxhill[i] = getRSS01(
    x = x_in3, x.idx = x_in3_idx, 
    x.new = as.vector(na.omit(boxhills[[3]]$design.list[[i]]$x.new)), 
    x.new.idx = as.vector(na.omit(boxhills[[3]]$design.list[[i]]$x.new.idx)), 
    y.new = as.vector(na.omit(boxhills[[3]]$design.list[[i]]$y.new)), 
    function.values = f.vals.mat[ , i], type = type01, l = l01, nugget = nuggetSM)
  RSS01_in4_boxhill[i] = getRSS01(
    x = random.mat[, i], x.idx = random.idx.mat[, i], 
    x.new = as.vector(na.omit(boxhills[[4]]$design.list[[i]]$x.new)), 
    x.new.idx = as.vector(na.omit(boxhills[[4]]$design.list[[i]]$x.new.idx)), 
    y.new = as.vector(na.omit(boxhills[[4]]$design.list[[i]]$y.new)), 
    function.values = f.vals.mat[ , i], type = type01, l = l01, nugget = nuggetSM)
}

# RSS01 for Space-filling
RSS01_in1_spacefilling = rep(NA, numSims)
RSS01_in2_spacefilling = rep(NA, numSims)
RSS01_in3_spacefilling = rep(NA, numSims)
RSS01_in4_spacefilling = rep(NA, numSims)
for(i in 1:numSims){
  RSS01_in1_spacefilling[i] = getRSS01(
    x = x_in1, x.idx = x_in1_idx, 
    x.new = x_spacefill1, x.new.idx = x_spacefill1_idx, 
    y.new = NULL, 
    function.values = f.vals.mat[ , i], type = type01, l = l01, nugget = nuggetSM)
  RSS01_in2_spacefilling[i] = getRSS01(
    x = x_in2, x.idx = x_in2_idx, 
    x.new = x_spacefill2, x.new.idx = x_spacefill2_idx,  
    y.new = NULL, 
    function.values = f.vals.mat[ , i], type = type01, l = l01, nugget = nuggetSM)
  RSS01_in3_spacefilling[i] = getRSS01(
    x = x_in3, x.idx = x_in3_idx,  
    x.new = x_spacefill3, x.new.idx = x_spacefill3_idx, 
    y.new = NULL, 
    function.values = f.vals.mat[ , i], type = type01, l = l01, nugget = nuggetSM)
  RSS01_in4_spacefilling[i] = getRSS01(
    x = random.mat[, i], x.idx = random.idx.mat[, i], 
    x.new = x_spacefill4, x.new.idx = x_spacefill4_idx, 
    y.new = NULL, 
    function.values = f.vals.mat[ , i], type = type01, l = l01, nugget = nuggetSM)
}

# RSS01 for Uniformly-distributed Random
RSS01_in1_random = rep(NA, numSims)
RSS01_in2_random = rep(NA, numSims)
RSS01_in3_random = rep(NA, numSims)
RSS01_in4_random = rep(NA, numSims)
for(i in 1:numSims){
  # first, get uniformly sampled random points
  x_uniform_idx = sample(1:numx, Nnew)
  x_uniform = x_seq[x_uniform_idx]
  # now calculate RSS01
  RSS01_in1_random[i] = getRSS01(
    x = x_in1, x.idx = x_in1_idx, 
    x.new = x_uniform, x.new.idx = x_uniform_idx, 
    y.new = NULL, 
    function.values = f.vals.mat[ , i], type = type01, l = l01, nugget = nuggetSM)
  RSS01_in2_random[i] = getRSS01(
    x = x_in2, x.idx = x_in2_idx, 
    x.new = x_uniform, x.new.idx = x_uniform_idx, 
    y.new = NULL, 
    function.values = f.vals.mat[ , i], type = type01, l = l01, nugget = nuggetSM)
  RSS01_in3_random[i] = getRSS01(
    x = x_in3, x.idx = x_in3_idx, 
    x.new = x_uniform, x.new.idx = x_uniform_idx, 
    y.new = NULL, 
    function.values = f.vals.mat[ , i], type = type01, l = l01, nugget = nuggetSM)
  RSS01_in4_random[i] = getRSS01(
    x = random.mat[, i], x.idx = random.idx.mat[, i], 
    x.new = x_uniform, x.new.idx = x_uniform_idx, 
    y.new = NULL, 
    function.values = f.vals.mat[ , i], type = type01, l = l01, nugget = nuggetSM)
}

# get the medians
RSS01_in1 = c(
  SeqMED = median(RSS01_in1_seqmed), 
  BoxHill = median(RSS01_in1_boxhill), 
  SpaceFill = median(RSS01_in1_spacefilling), 
  Random = median(RSS01_in1_random)
  )
RSS01_in2 = c(
  SeqMED = median(RSS01_in2_seqmed), 
  BoxHill = median(RSS01_in2_boxhill), 
  SpaceFill = median(RSS01_in2_spacefilling), 
  Random = median(RSS01_in2_random)
)
RSS01_in3 = c(
  SeqMED = median(RSS01_in3_seqmed), 
  BoxHill = median(RSS01_in3_boxhill), 
  SpaceFill = median(RSS01_in3_spacefilling), 
  Random = median(RSS01_in3_random)
)
RSS01_in4 = c(
  SeqMED = median(RSS01_in4_seqmed), 
  BoxHill = median(RSS01_in4_boxhill), 
  SpaceFill = median(RSS01_in4_spacefilling), 
  Random = median(RSS01_in4_random)
)

# plot
ggdata = data.table(
  `Extrapolation` = RSS01_in1,
  `Inc Spread` = RSS01_in2,
  `Even Coverage` = RSS01_in3,
  `Random` = RSS01_in4,
  Design = c("SeqMED", "Boxhill", "SpaceFill", "Random")
)
ggdata = melt(ggdata, id.vars = c("Design"))
RSS01.plt = ggplot(ggdata, aes(x = variable, y = value, group = Design, 
                   color = Design)) + 
  geom_point() + 
  geom_path(linetype = 2) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.text.x = element_text(size = 5)) +
  labs(y = "RSS0/RSS1", x = "Initial Data")
RSS01.plt
# ggsave("gvm_medianlogrss01.pdf",
#        plot = last_plot(),
#        device = "pdf",
#        path = image_path,
#        scale = 1,
#        width = 4.5,
#        height = 2,
#        units = c("in")
# )

# Posterior Probability of H1

# see boxhill_gp.R function Evidence_gp()
# also see boxhill.R function getHypothesesPosteriors()

# PPH1 for SeqMED
PPH1_in1_seqmed = rep(NA, numSims)
PPH1_in2_seqmed = rep(NA, numSims)
PPH1_in3_seqmed = rep(NA, numSims)
PPH1_in4_seqmed = rep(NA, numSims)
for(i in 1:numSims){
  # for input 1
  y.tmp = c(seqmeds[[1]]$design.list[[i]]$y, seqmeds[[1]]$design.list[[i]]$y.new)
  x.tmp = c(seqmeds[[1]]$design.list[[i]]$x, seqmeds[[1]]$design.list[[i]]$x.new)
  PPH1_in1_seqmed[i] = getHypothesesPosteriors(
    prior.probs = prior_probs, 
    evidences = c(
      Evidence_gp(y.tmp, x.tmp, model0, nuggetSM),
      Evidence_gp(y.tmp, x.tmp, model1, nuggetSM)
    )
  )[2]
  # for input 2
  y.tmp = c(seqmeds[[2]]$design.list[[i]]$y, seqmeds[[2]]$design.list[[i]]$y.new)
  x.tmp = c(seqmeds[[2]]$design.list[[i]]$x, seqmeds[[2]]$design.list[[i]]$x.new)
  PPH1_in2_seqmed[i] = getHypothesesPosteriors(
    prior.probs = prior_probs, 
    evidences = c(
      Evidence_gp(y.tmp, x.tmp, model0, nuggetSM),
      Evidence_gp(y.tmp, x.tmp, model1, nuggetSM)
    )
  )[2]
  # for input 3
  y.tmp = c(seqmeds[[3]]$design.list[[i]]$y, seqmeds[[3]]$design.list[[i]]$y.new)
  x.tmp = c(seqmeds[[3]]$design.list[[i]]$x, seqmeds[[3]]$design.list[[i]]$x.new)
  PPH1_in3_seqmed[i] = getHypothesesPosteriors(
    prior.probs = prior_probs, 
    evidences = c(
      Evidence_gp(y.tmp, x.tmp, model0, nuggetSM),
      Evidence_gp(y.tmp, x.tmp, model1, nuggetSM)
    )
  )[2]
  # for input 4
  y.tmp = c(seqmeds[[4]]$design.list[[i]]$y, seqmeds[[4]]$design.list[[i]]$y.new)
  x.tmp = c(seqmeds[[4]]$design.list[[i]]$x, seqmeds[[4]]$design.list[[i]]$x.new)
  PPH1_in4_seqmed[i] = getHypothesesPosteriors(
    prior.probs = prior_probs, 
    evidences = c(
      Evidence_gp(y.tmp, x.tmp, model0, nuggetSM),
      Evidence_gp(y.tmp, x.tmp, model1, nuggetSM)
    )
  )[2]
}

# PPH1 for Box-Hill
PPH1_in1_boxhill = rep(NA, numSims)
PPH1_in2_boxhill = rep(NA, numSims)
PPH1_in3_boxhill = rep(NA, numSims)
PPH1_in4_boxhill = rep(NA, numSims)
for(i in 1:numSims){
  # input 1
  y.tmp = c(boxhills[[1]]$design.list[[i]]$y, 
            as.vector(na.omit(boxhills[[1]]$design.list[[i]]$y.new)))
  x.tmp = c(boxhills[[1]]$design.list[[i]]$x, 
            as.vector(na.omit(boxhills[[1]]$design.list[[i]]$x.new)))
  PPH1_in1_boxhill[i] = getHypothesesPosteriors(
    prior.probs = prior_probs, 
    evidences = c(
      Evidence_gp(y.tmp, x.tmp, model0, nuggetBH),
      Evidence_gp(y.tmp, x.tmp, model1, nuggetBH)
    )
  )[2]
  # input 2
  y.tmp = c(boxhills[[2]]$design.list[[i]]$y, 
            as.vector(na.omit(boxhills[[2]]$design.list[[i]]$y.new)))
  x.tmp = c(boxhills[[2]]$design.list[[i]]$x, 
            as.vector(na.omit(boxhills[[2]]$design.list[[i]]$x.new)))
  PPH1_in2_boxhill[i] = getHypothesesPosteriors(
    prior.probs = prior_probs, 
    evidences = c(
      Evidence_gp(y.tmp, x.tmp, model0, nuggetBH),
      Evidence_gp(y.tmp, x.tmp, model1, nuggetBH)
    )
  )[2]
  # for input 3
  y.tmp = c(boxhills[[3]]$design.list[[i]]$y, 
            as.vector(na.omit(boxhills[[3]]$design.list[[i]]$y.new)))
  x.tmp = c(boxhills[[3]]$design.list[[i]]$x, 
            as.vector(na.omit(boxhills[[3]]$design.list[[i]]$x.new)))
  PPH1_in3_boxhill[i] = getHypothesesPosteriors(
    prior.probs = prior_probs, 
    evidences = c(
      Evidence_gp(y.tmp, x.tmp, model0, nuggetBH),
      Evidence_gp(y.tmp, x.tmp, model1, nuggetBH)
    )
  )[2]
  # input 4
  y.tmp = c(boxhills[[4]]$design.list[[i]]$y, 
            as.vector(na.omit(boxhills[[4]]$design.list[[i]]$y.new)))
  x.tmp = c(boxhills[[4]]$design.list[[i]]$x, 
            as.vector(na.omit(boxhills[[4]]$design.list[[i]]$x.new)))
  PPH1_in4_boxhill[i] = getHypothesesPosteriors(
    prior.probs = prior_probs, 
    evidences = c(
      Evidence_gp(y.tmp, x.tmp, model0, nuggetBH),
      Evidence_gp(y.tmp, x.tmp, model1, nuggetBH)
    )
  )[2]
}

# PPH1 for Space-filling
PPH1_in1_spacefilling = rep(NA, numSims)
PPH1_in2_spacefilling = rep(NA, numSims)
PPH1_in3_spacefilling = rep(NA, numSims)
PPH1_in4_spacefilling = rep(NA, numSims)
nugget = nuggetSM
for(i in 1:numSims){
  # input 1
  fn.vals.tmp = seqmeds[[1]]$function.values.list[ , i]
  y.tmp = c(fn.vals.tmp[seqmeds[[1]]$x0.idx], fn.vals.tmp[x_spacefill1_idx])
  x.tmp = c(seqmeds[[1]]$x0, x_spacefill1)
  PPH1_in1_spacefilling[i] = getHypothesesPosteriors(
    prior.probs = prior_probs, 
    evidences = c(
      Evidence_gp(y.tmp, x.tmp, model0, nugget),
      Evidence_gp(y.tmp, x.tmp, model1, nugget)
    )
  )[2]
  # input 2
  fn.vals.tmp = seqmeds[[2]]$function.values.list[ , i]
  y.tmp = c(fn.vals.tmp[seqmeds[[2]]$x0.idx], fn.vals.tmp[x_spacefill2_idx])
  x.tmp = c(seqmeds[[2]]$x0, x_spacefill2)
  PPH1_in2_spacefilling[i] = getHypothesesPosteriors(
    prior.probs = prior_probs, 
    evidences = c(
      Evidence_gp(y.tmp, x.tmp, model0, nugget),
      Evidence_gp(y.tmp, x.tmp, model1, nugget)
    )
  )[2]
  # input 3
  fn.vals.tmp = seqmeds[[3]]$function.values.list[ , i]
  y.tmp = c(fn.vals.tmp[seqmeds[[3]]$x0.idx], fn.vals.tmp[x_spacefill3_idx])
  x.tmp = c(seqmeds[[3]]$x0, x_spacefill3)
  PPH1_in3_spacefilling[i] = getHypothesesPosteriors(
    prior.probs = prior_probs, 
    evidences = c(
      Evidence_gp(y.tmp, x.tmp, model0, nugget),
      Evidence_gp(y.tmp, x.tmp, model1, nugget)
    )
  )[2]
  # input 4
  fn.vals.tmp = seqmeds[[4]]$function.values.list[ , i]
  y.tmp = c(fn.vals.tmp[seqmeds[[4]]$design.list[[i]]$x.idx], 
            fn.vals.tmp[x_spacefill4_idx])
  x.tmp = c(seqmeds[[4]]$design.list[[i]]$x, x_spacefill4)
  PPH1_in4_spacefilling[i] = getHypothesesPosteriors(
    prior.probs = prior_probs, 
    evidences = c(
      Evidence_gp(y.tmp, x.tmp, model0, nugget),
      Evidence_gp(y.tmp, x.tmp, model1, nugget)
    )
  )[2]
}

# PPH1 for Uniformly-distributed Random
PPH1_in1_random = rep(NA, numSims)
PPH1_in2_random = rep(NA, numSims)
PPH1_in3_random = rep(NA, numSims)
PPH1_in4_random = rep(NA, numSims)
for(i in 1:numSims){
  # first, get uniformly sampled random points
  x_uniform_idx = sample(1:numx, Nnew)
  x_uniform = x_seq[x_uniform_idx]
  # now calculate PPH1
  # input 1
  fn.vals.tmp = seqmeds[[1]]$function.values.list[ , i]
  y.tmp = c(fn.vals.tmp[seqmeds[[1]]$x0.idx], fn.vals.tmp[x_uniform_idx])
  x.tmp = c(seqmeds[[1]]$x0, x_uniform)
  PPH1_in1_random[i] = getHypothesesPosteriors(
    prior.probs = prior_probs, 
    evidences = c(
      Evidence_gp(y.tmp, x.tmp, model0, nugget),
      Evidence_gp(y.tmp, x.tmp, model1, nugget)
    )
  )[2]
  # input 2
  fn.vals.tmp = seqmeds[[2]]$function.values.list[ , i]
  y.tmp = c(fn.vals.tmp[seqmeds[[2]]$x0.idx], fn.vals.tmp[x_uniform_idx])
  x.tmp = c(seqmeds[[2]]$x0, x_uniform)
  PPH1_in2_random[i] = getHypothesesPosteriors(
    prior.probs = prior_probs, 
    evidences = c(
      Evidence_gp(y.tmp, x.tmp, model0, nugget),
      Evidence_gp(y.tmp, x.tmp, model1, nugget)
    )
  )[2]
  # input 3
  fn.vals.tmp = seqmeds[[3]]$function.values.list[ , i]
  y.tmp = c(fn.vals.tmp[seqmeds[[3]]$x0.idx], fn.vals.tmp[x_uniform_idx])
  x.tmp = c(seqmeds[[3]]$x0, x_uniform)
  PPH1_in3_random[i] = getHypothesesPosteriors(
    prior.probs = prior_probs, 
    evidences = c(
      Evidence_gp(y.tmp, x.tmp, model0, nugget),
      Evidence_gp(y.tmp, x.tmp, model1, nugget)
    )
  )[2]
  # input 4
  fn.vals.tmp = seqmeds[[4]]$function.values.list[ , i]
  y.tmp = c(fn.vals.tmp[seqmeds[[4]]$design.list[[i]]$x.idx], 
            fn.vals.tmp[x_uniform_idx])
  x.tmp = c(seqmeds[[4]]$design.list[[i]]$x, x_uniform)
  PPH1_in4_random[i] = getHypothesesPosteriors(
    prior.probs = prior_probs, 
    evidences = c(
      Evidence_gp(y.tmp, x.tmp, model0, nugget),
      Evidence_gp(y.tmp, x.tmp, model1, nugget)
    )
  )[2]
}

# get the medians
PPH1_in1 = c(
  SeqMED = median(PPH1_in1_seqmed), 
  BoxHill = median(PPH1_in1_boxhill), 
  SpaceFill = median(PPH1_in1_spacefilling), 
  Random = median(PPH1_in1_random)
)
PPH1_in2 = c(
  SeqMED = median(PPH1_in2_seqmed), 
  BoxHill = median(PPH1_in2_boxhill), 
  SpaceFill = median(PPH1_in2_spacefilling), 
  Random = median(PPH1_in2_random)
)
PPH1_in3 = c(
  SeqMED = median(PPH1_in3_seqmed), 
  BoxHill = median(PPH1_in3_boxhill), 
  SpaceFill = median(PPH1_in3_spacefilling), 
  Random = median(PPH1_in3_random)
)
PPH1_in4 = c(
  SeqMED = median(PPH1_in4_seqmed), 
  BoxHill = median(PPH1_in4_boxhill), 
  SpaceFill = median(PPH1_in4_spacefilling), 
  Random = median(PPH1_in4_random)
)

ggdata = data.table(
  Extrapolation = PPH1_in1,
  `Inc Spread` = PPH1_in2,
  `Even Coverage` = PPH1_in3,
  `Random` = PPH1_in4,
  Design = c("SeqMED", "BoxHill", "SpaceFill", "Random")
)
ggdata = melt(ggdata, id.vars = c("Design"))
PPH1.plt = ggplot(ggdata, aes(x = variable, y = value, group = Design, 
                   color = Design)) + 
  geom_point() + 
  geom_path(linetype = 2) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.text.x = element_text(size = 5)) +
  labs(y = "P(H1|X, Y)", x = "Initial Data")
PPH1.plt
# ggsave("gvm_medianpph1.pdf",
#        plot = last_plot(),
#        device = "pdf",
#        path = image_path,
#        scale = 1,
#        width = 4.5,
#        height = 2,
#        units = c("in")
# )

# in the same plot
ggdata1 = data.table(
  `Extrapolation` = RSS01_in1,
  `Inc Spread` = RSS01_in2,
  `Even Coverage` = RSS01_in3,
  `Random` = RSS01_in4,
  Design = c("SeqMED", "BoxHill", "SpaceFill", "Random")
)
ggdata1 = melt(ggdata1, id.vars = c("Design"))
ggdata1$Metric = "RSS0/RSS1"
ggdata2 = data.table(
  `Extrapolation` = PPH1_in1,
  `Inc Spread` = PPH1_in2,
  `Even Coverage` = PPH1_in3,
  `Random` = PPH1_in4,
  Design = c("SeqMED", "BoxHill", "SpaceFill", "Random")
)


ggdata2 = melt(ggdata2, id.vars = c("Design"))
ggdata2$Metric = "P(H1|X,Y)"
ggdata3 = rbind(ggdata1, ggdata2)
ggdata3$Metric = factor(ggdata3$Metric)
plt_gvm3 = ggplot(ggdata3, aes(x = variable, y = value, group = Design, 
                               color = Design)) + 
  facet_wrap(vars(Metric), scales = "free") + 
  geom_point() + 
  geom_path(linetype = 2) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.title = element_blank(), 
        axis.text.x = element_text(angle = 45, vjust = 0.5))
plt_gvm3
# ggsave("gvm_medianlogrss01pph1.pdf",
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
# Scenario 2: matern vs. periodic, true = periodic
################################################################################
type01 = c("matern", "periodic")
l01= c(0.01, 0.1)
# generate matern functions
set.seed(seed)
null_cov = getCov(x_seq, x_seq, type01[2], l01[2])
null_mean = rep(0, numx)

# bh settings
model0 = list(type = type01[1], l = l01[1])
model1 = list(type = type01[2], l = l01[2])

################################################################################
# Plot the demo design
################################################################################

# demo settings
l01.demo = c(0.1, 0.5)
nugget.demo = NULL

# read in sims
maternvsperiodic_train4sims = readRDS(
  paste0(
    "run_designs/gp",
    "/gp_demo/maternvsperiodic_train4sims.rds"))

# get sim info
sim_ind = 1
x_in = maternvsperiodic_train4sims$x_train[ , sim_ind]
x_in_idx = maternvsperiodic_train4sims$x_train_ind[ , sim_ind]
mmed_gp = maternvsperiodic_train4sims$mmed_gp_list[[sim_ind]]
y_seq = maternvsperiodic_train4sims$sim_fns[ , sim_ind]
y_in = y_seq[x_in_idx]

newpts = mmed_gp$addD
truey = y_seq[mmed_gp$indices]

H0_predfn = getGPPredictive(x_seq, x_in, y_in, type01[1], l01.demo[1],
                            nugget = nugget.demo)
H1_predfn = getGPPredictive(x_seq, x_in, y_in, type01[2], l01.demo[2],
                            nugget = nugget.demo)

# get w_seq
Kinv0 = solve(getCov(x_in, x_in, type01[1], l01.demo[1]))
Kinv1 = solve(getCov(x_in, x_in, type01[2], l01.demo[2]))
w_seq = sapply(x_seq, FUN = function(x1) 
  WNgp(x1, Kinv0, Kinv1, x_in, y_in, var_e = 1, type01, l01.demo))

# plot
err0 = 2 * sqrt(diag(H0_predfn$pred_var))
err1 = 2 * sqrt(diag(H1_predfn$pred_var))
ggdata = data.table(
  x = x_seq, 
  `True Function` = y_seq, 
  Wasserstein = w_seq, 
  `H0 Predictive` = H0_predfn$pred_mean, 
  `H1 Predictive` = H1_predfn$pred_mean,
  lower0 = H0_predfn$pred_mean - err0, 
  lower1 = H1_predfn$pred_mean - err1, 
  upper0 = H0_predfn$pred_mean + err0,
  upper1 = H1_predfn$pred_mean + err1
)
yrange = range(ggdata$lower0, ggdata$lower1, 
               ggdata$upper0, ggdata$upper1)
yrange[1] = yrange[1] - 1
ggdata$Wasserstein = ggdata$Wasserstein * 0.25 - abs(yrange[1])
ggdata$zero1 = NA
ggdata$zero2 = NA
ggdata.melted = melt(ggdata, id.vars = c("x"), 
                     measure.vars = c("True Function", "Wasserstein", "H0 Predictive", "H1 Predictive"))
ggdata.lower = melt(ggdata, id.vars = c("x"), 
                    measure.vars = c("zero1", "zero2", "lower0", "lower1"))
ggdata.upper = melt(ggdata, id.vars = c("x"), 
                    measure.vars = c("zero1", "zero2", "upper0", "upper1"))
ggdata.melted = cbind(ggdata.melted, 
                      lower = ggdata.lower$value, 
                      upper = ggdata.upper$value)
ggdata_pts = data.table(
  x = c(x_in, newpts), 
  y = c(y_in, truey), 
  color = c(rep(gg_color_hue(2)[2], length(x_in)), 
            rep(gg_color_hue(2)[1], length(newpts))), 
  shape = c(rep(8, length(x_in)), 
            rep(16, length(newpts)))
)
ggplot(data = ggdata.melted, aes(x = x, y =value, color = variable), 
       linetype = 1) + 
  geom_path() + 
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = variable), 
              alpha = 0.1, linetype = 0) +
  scale_linetype_manual(values = c(1, 1, 2, 2)) + 
  scale_fill_manual(values = c(NA, NA, "#00BFC4", "#C77CFF")) + 
  scale_color_manual(values = c(1, "gray", "#00BFC4", "#C77CFF")) + 
  geom_point(data = ggdata_pts, mapping = aes(x = x, y = y), 
             inherit.aes = FALSE, color = ggdata_pts$color, 
             shape = ggdata_pts$shape,
             size = 2) +
  geom_point(data = ggdata_pts, mapping = aes(x = x, y = yrange[1]), 
             inherit.aes = FALSE, color = ggdata_pts$color, 
             shape = ggdata_pts$shape, 
             size = 2) +
  scale_y_continuous(limits = yrange) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  labs(y = "y", x = "x", fill = "Function", color = "Function")
# ggsave("mvp.pdf",
#        plot = last_plot(),
#        device = "pdf",
#        path = image_path,
#        scale = 1,
#        width = 4.5,
#        height = 2.5,
#        units = c("in")
# )

################################################################################
# Metrics
################################################################################

seqmeds = list()
for(i in 1:4){
  seqmeds[[i]] = readRDS(
    paste0(
      output_home, 
      "/seqmed/scenario2_seqmed_simulations", 
      "_input", i, 
      "_Nin6_Nnew15_numSims100.rds"))
}

boxhills = list()
for(i in 1:4){
  boxhills[[i]] = readRDS(
    paste0(
      output_home, 
      "/boxhill/scenario2_boxhill_simulations", 
      "_input", i, 
      "_Nin6_Nnew15_numSims100.rds"))
}

# RSS Ratio (0/1)

# see helper function getRSS01

# # check that seqmeds and boxhills use the same values # # # # # # # # # # # #
# seqmeds.function.values = matrix(NA, nrow = numx * 4, ncol = numSims)
# boxhills.function.values = matrix(NA, nrow = numx * 4, ncol = numSims)
# for(i in 1:numSims){
#   seqmeds.function.values[, i] = c(
#     seqmeds[[1]]$function.values.list[ , i],
#     seqmeds[[2]]$function.values.list[ , i],
#     seqmeds[[3]]$function.values.list[ , i],
#     seqmeds[[4]]$function.values.list[ , i]
#   )
#   boxhills.function.values[, i] = c(
#     boxhills[[1]]$function.values.list[ , i],
#     boxhills[[2]]$function.values.list[ , i],
#     boxhills[[3]]$function.values.list[ , i],
#     boxhills[[4]]$function.values.list[ , i]
#   )
# }
# all.equal(seqmeds.function.values, boxhills.function.values) # TRUE!
# 
# # check that each of seqmeds function.values.lists are equal # # # # # # # # #
# seqmed.function.values1 = matrix(NA, nrow = numx, ncol = numSims)
# seqmed.function.values2 = matrix(NA, nrow = numx, ncol = numSims)
# seqmed.function.values3 = matrix(NA, nrow = numx, ncol = numSims)
# seqmed.function.values4 = matrix(NA, nrow = numx, ncol = numSims)
# for(i in 1:numSims){
#   seqmed.function.values1[, i] = seqmeds[[1]]$function.values.list[ , i]
#   seqmed.function.values2[, i] = seqmeds[[2]]$function.values.list[ , i]
#   seqmed.function.values3[, i] = seqmeds[[3]]$function.values.list[ , i]
#   seqmed.function.values4[, i] = seqmeds[[4]]$function.values.list[ , i]
# }
# all.equal(seqmed.function.values1, seqmed.function.values2) # TRUE!
# all.equal(seqmed.function.values1, seqmed.function.values3) # TRUE!
# all.equal(seqmed.function.values1, seqmed.function.values4) # TRUE!
f.vals.mat = seqmeds[[1]]$function.values.list

# # check that random designs are the same in the seqmed and boxhill sims # # # 
# seqmeds.random.mat = matrix(NA, nrow = Nin, ncol = numSims)
# boxhills.random.mat = matrix(NA, nrow = Nin, ncol = numSims)
# seqmeds.random.idx.mat = matrix(NA, nrow = Nin, ncol = numSims)
# boxhills.random.idx.mat = matrix(NA, nrow = Nin, ncol = numSims)
# for(i in 1:numSims){
#   seqmeds.random.mat[, i] = seqmeds[[4]]$design.list[[i]]$x
#   seqmeds.random.idx.mat[, i] = seqmeds[[4]]$design.list[[i]]$x.idx
#   boxhills.random.mat[, i] = seqmeds[[4]]$design.list[[i]]$x
#   boxhills.random.idx.mat[, i] = seqmeds[[4]]$design.list[[i]]$x.idx
# }
# all.equal(seqmeds.random.mat, boxhills.random.mat) # TRUE!
# all.equal(seqmeds.random.idx.mat, boxhills.random.idx.mat) # TRUE!
random.mat = matrix(NA, nrow = Nin, ncol = numSims)
random.idx.mat = matrix(NA, nrow = Nin, ncol = numSims)
for(i in 1:numSims){
  random.mat[, i] = seqmeds[[4]]$design.list[[i]]$x
  random.idx.mat[, i] = seqmeds[[4]]$design.list[[i]]$x.idx
}

# RSS01 for SeqMED
RSS01_in1_seqmed = rep(NA, numSims)
RSS01_in2_seqmed = rep(NA, numSims)
RSS01_in3_seqmed = rep(NA, numSims)
RSS01_in4_seqmed = rep(NA, numSims)

for(i in 1:numSims){
  RSS01_in1_seqmed[i] = getRSS01(
    x = x_in1, x.idx = x_in1_idx, 
    x.new = seqmeds[[1]]$design.list[[i]]$x.new, 
    x.new.idx = seqmeds[[1]]$design.list[[i]]$x.new.idx, 
    y.new = seqmeds[[1]]$design.list[[i]]$y.new, 
    function.values = f.vals.mat[ , i], type = type01, l = l01, nugget = nuggetSM)
  # plot(x_seq, seqmeds[[1]]$function.values.list[ , i], type = "l")
  # points(seqmeds[[1]]$x0, seqmeds[[1]]$design.list[[i]]$y[1:Nin], col = 3)
  # points(seqmeds[[1]]$design.list[[i]]$D[-c(1:Nin)], 
  #        seqmeds[[1]]$design.list[[i]]$y[-c(1:Nin)], col = 2)
  RSS01_in2_seqmed[i] = getRSS01(
    x = x_in2, x.idx = x_in2_idx, 
    x.new = seqmeds[[2]]$design.list[[i]]$x.new, 
    x.new.idx = seqmeds[[2]]$design.list[[i]]$x.new.idx, 
    y.new = seqmeds[[2]]$design.list[[i]]$y.new, 
    function.values = f.vals.mat[ , i], type = type01, l = l01, nugget = nuggetSM)
  RSS01_in3_seqmed[i] = getRSS01(
    x = x_in3, x.idx = x_in3_idx, 
    x.new = seqmeds[[3]]$design.list[[i]]$x.new, 
    x.new.idx = seqmeds[[3]]$design.list[[i]]$x.new.idx, 
    y.new = seqmeds[[3]]$design.list[[i]]$y.new, 
    function.values = f.vals.mat[ , i], type = type01, l = l01, nugget = nuggetSM)
  RSS01_in4_seqmed[i] = getRSS01(
    x = random.mat[, i], x.idx = random.idx.mat[, i], 
    x.new = seqmeds[[4]]$design.list[[i]]$x.new, 
    x.new.idx = seqmeds[[4]]$design.list[[i]]$x.new.idx, 
    y.new = seqmeds[[4]]$design.list[[i]]$y.new, 
    function.values = f.vals.mat[ , i], type = type01, l = l01, nugget = nuggetSM)
}

# RSS01 for Box-Hill
RSS01_in1_boxhill = rep(NA, numSims)
RSS01_in2_boxhill = rep(NA, numSims)
RSS01_in3_boxhill = rep(NA, numSims)
RSS01_in4_boxhill = rep(NA, numSims)
for(i in 1:numSims){
  RSS01_in1_boxhill[i] = getRSS01(
    x = x_in1, x.idx = x_in1_idx, 
    x.new = as.vector(na.omit(boxhills[[1]]$design.list[[i]]$x.new)), 
    x.new.idx = as.vector(na.omit(boxhills[[1]]$design.list[[i]]$x.new.idx)), 
    y.new = as.vector(na.omit(boxhills[[1]]$design.list[[i]]$y.new)), 
    function.values = f.vals.mat[ , i], type = type01, l = l01, nugget = nuggetSM)
  RSS01_in2_boxhill[i] = getRSS01(
    x = x_in2, x.idx = x_in2_idx, 
    x.new = as.vector(na.omit(boxhills[[2]]$design.list[[i]]$x.new)), 
    x.new.idx = as.vector(na.omit(boxhills[[2]]$design.list[[i]]$x.new.idx)), 
    y.new = as.vector(na.omit(boxhills[[2]]$design.list[[i]]$y.new)), 
    function.values = f.vals.mat[ , i], type = type01, l = l01, nugget = nuggetSM)
  RSS01_in3_boxhill[i] = getRSS01(
    x = x_in3, x.idx = x_in3_idx, 
    x.new = as.vector(na.omit(boxhills[[3]]$design.list[[i]]$x.new)), 
    x.new.idx = as.vector(na.omit(boxhills[[3]]$design.list[[i]]$x.new.idx)), 
    y.new = as.vector(na.omit(boxhills[[3]]$design.list[[i]]$y.new)), 
    function.values = f.vals.mat[ , i], type = type01, l = l01, nugget = nuggetSM)
  RSS01_in4_boxhill[i] = getRSS01(
    x = random.mat[, i], x.idx = random.idx.mat[, i], 
    x.new = as.vector(na.omit(boxhills[[4]]$design.list[[i]]$x.new)), 
    x.new.idx = as.vector(na.omit(boxhills[[4]]$design.list[[i]]$x.new.idx)), 
    y.new = as.vector(na.omit(boxhills[[4]]$design.list[[i]]$y.new)), 
    function.values = f.vals.mat[ , i], type = type01, l = l01, nugget = nuggetSM)
}

# RSS01 for Space-filling
RSS01_in1_spacefilling = rep(NA, numSims)
RSS01_in2_spacefilling = rep(NA, numSims)
RSS01_in3_spacefilling = rep(NA, numSims)
RSS01_in4_spacefilling = rep(NA, numSims)
for(i in 1:numSims){
  RSS01_in1_spacefilling[i] = getRSS01(
    x = x_in1, x.idx = x_in1_idx, 
    x.new = x_spacefill1, x.new.idx = x_spacefill1_idx, 
    y.new = NULL, 
    function.values = f.vals.mat[ , i], type = type01, l = l01, nugget = nuggetSM)
  RSS01_in2_spacefilling[i] = getRSS01(
    x = x_in2, x.idx = x_in2_idx, 
    x.new = x_spacefill2, x.new.idx = x_spacefill2_idx,  
    y.new = NULL, 
    function.values = f.vals.mat[ , i], type = type01, l = l01, nugget = nuggetSM)
  RSS01_in3_spacefilling[i] = getRSS01(
    x = x_in3, x.idx = x_in3_idx,  
    x.new = x_spacefill3, x.new.idx = x_spacefill3_idx, 
    y.new = NULL, 
    function.values = f.vals.mat[ , i], type = type01, l = l01, nugget = nuggetSM)
  RSS01_in4_spacefilling[i] = getRSS01(
    x = random.mat[, i], x.idx = random.idx.mat[, i], 
    x.new = x_spacefill4, x.new.idx = x_spacefill4_idx, 
    y.new = NULL, 
    function.values = f.vals.mat[ , i], type = type01, l = l01, nugget = nuggetSM)
}

# RSS01 for Uniformly-distributed Random
RSS01_in1_random = rep(NA, numSims)
RSS01_in2_random = rep(NA, numSims)
RSS01_in3_random = rep(NA, numSims)
RSS01_in4_random = rep(NA, numSims)
for(i in 1:numSims){
  # first, get uniformly sampled random points
  x_uniform_idx = sample(1:numx, Nnew)
  x_uniform = x_seq[x_uniform_idx]
  # now calculate RSS01
  RSS01_in1_random[i] = getRSS01(
    x = x_in1, x.idx = x_in1_idx, 
    x.new = x_uniform, x.new.idx = x_uniform_idx, 
    y.new = NULL, 
    function.values = f.vals.mat[ , i], type = type01, l = l01, nugget = nuggetSM)
  RSS01_in2_random[i] = getRSS01(
    x = x_in2, x.idx = x_in2_idx, 
    x.new = x_uniform, x.new.idx = x_uniform_idx, 
    y.new = NULL, 
    function.values = f.vals.mat[ , i], type = type01, l = l01, nugget = nuggetSM)
  RSS01_in3_random[i] = getRSS01(
    x = x_in3, x.idx = x_in3_idx, 
    x.new = x_uniform, x.new.idx = x_uniform_idx, 
    y.new = NULL, 
    function.values = f.vals.mat[ , i], type = type01, l = l01, nugget = nuggetSM)
  RSS01_in4_random[i] = getRSS01(
    x = random.mat[, i], x.idx = random.idx.mat[, i], 
    x.new = x_uniform, x.new.idx = x_uniform_idx, 
    y.new = NULL, 
    function.values = f.vals.mat[ , i], type = type01, l = l01, nugget = nuggetSM)
}

# get the medians
RSS01_in1 = c(
  SeqMED = median(RSS01_in1_seqmed), 
  BoxHill = median(RSS01_in1_boxhill), 
  SpaceFill = median(RSS01_in1_spacefilling), 
  Random = median(RSS01_in1_random)
)
RSS01_in2 = c(
  SeqMED = median(RSS01_in2_seqmed), 
  BoxHill = median(RSS01_in2_boxhill), 
  SpaceFill = median(RSS01_in2_spacefilling), 
  Random = median(RSS01_in2_random)
)
RSS01_in3 = c(
  SeqMED = median(RSS01_in3_seqmed), 
  BoxHill = median(RSS01_in3_boxhill), 
  SpaceFill = median(RSS01_in3_spacefilling), 
  Random = median(RSS01_in3_random)
)
RSS01_in4 = c(
  SeqMED = median(RSS01_in4_seqmed), 
  BoxHill = median(RSS01_in4_boxhill), 
  SpaceFill = median(RSS01_in4_spacefilling), 
  Random = median(RSS01_in4_random)
)

# plot
ggdata = data.table(
  `Extrapolation` = RSS01_in1,
  `Inc Spread` = RSS01_in2,
  `Even Coverage` = RSS01_in3,
  `Random` = RSS01_in4,
  Design = c("SeqMED", "Boxhill", "SpaceFill", "Random")
)
ggdata = melt(ggdata, id.vars = c("Design"))
RSS01.plt = ggplot(ggdata, aes(x = variable, y = value, group = Design, 
                               color = Design)) + 
  geom_point() + 
  geom_path(linetype = 2) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.text.x = element_text(size = 5)) +
  labs(y = "RSS0/RSS1", x = "Initial Data")
RSS01.plt
# ggsave("mvp_medianlogrss01.pdf",
#        plot = last_plot(),
#        device = "pdf",
#        path = image_path,
#        scale = 1,
#        width = 4.5,
#        height = 2,
#        units = c("in")
# )

# Posterior Probability of H1

# see boxhill_gp.R function Evidence_gp()
# also see boxhill.R function getHypothesesPosteriors()

# PPH1 for SeqMED
PPH1_in1_seqmed = rep(NA, numSims)
PPH1_in2_seqmed = rep(NA, numSims)
PPH1_in3_seqmed = rep(NA, numSims)
PPH1_in4_seqmed = rep(NA, numSims)
for(i in 1:numSims){
  # for input 1
  y.tmp = c(seqmeds[[1]]$design.list[[i]]$y, seqmeds[[1]]$design.list[[i]]$y.new)
  x.tmp = c(seqmeds[[1]]$design.list[[i]]$x, seqmeds[[1]]$design.list[[i]]$x.new)
  PPH1_in1_seqmed[i] = getHypothesesPosteriors(
    prior.probs = prior_probs, 
    evidences = c(
      Evidence_gp(y.tmp, x.tmp, model0, nuggetSM),
      Evidence_gp(y.tmp, x.tmp, model1, nuggetSM)
    )
  )[2]
  # for input 2
  y.tmp = c(seqmeds[[2]]$design.list[[i]]$y, seqmeds[[2]]$design.list[[i]]$y.new)
  x.tmp = c(seqmeds[[2]]$design.list[[i]]$x, seqmeds[[2]]$design.list[[i]]$x.new)
  PPH1_in2_seqmed[i] = getHypothesesPosteriors(
    prior.probs = prior_probs, 
    evidences = c(
      Evidence_gp(y.tmp, x.tmp, model0, nuggetSM),
      Evidence_gp(y.tmp, x.tmp, model1, nuggetSM)
    )
  )[2]
  # for input 3
  y.tmp = c(seqmeds[[3]]$design.list[[i]]$y, seqmeds[[3]]$design.list[[i]]$y.new)
  x.tmp = c(seqmeds[[3]]$design.list[[i]]$x, seqmeds[[3]]$design.list[[i]]$x.new)
  PPH1_in3_seqmed[i] = getHypothesesPosteriors(
    prior.probs = prior_probs, 
    evidences = c(
      Evidence_gp(y.tmp, x.tmp, model0, nuggetSM),
      Evidence_gp(y.tmp, x.tmp, model1, nuggetSM)
    )
  )[2]
  # for input 4
  y.tmp = c(seqmeds[[4]]$design.list[[i]]$y, seqmeds[[4]]$design.list[[i]]$y.new)
  x.tmp = c(seqmeds[[4]]$design.list[[i]]$x, seqmeds[[4]]$design.list[[i]]$x.new)
  PPH1_in4_seqmed[i] = getHypothesesPosteriors(
    prior.probs = prior_probs, 
    evidences = c(
      Evidence_gp(y.tmp, x.tmp, model0, nuggetSM),
      Evidence_gp(y.tmp, x.tmp, model1, nuggetSM)
    )
  )[2]
}

# PPH1 for Box-Hill
PPH1_in1_boxhill = rep(NA, numSims)
PPH1_in2_boxhill = rep(NA, numSims)
PPH1_in3_boxhill = rep(NA, numSims)
PPH1_in4_boxhill = rep(NA, numSims)
for(i in 1:numSims){
  # input 1
  y.tmp = c(boxhills[[1]]$design.list[[i]]$y, 
            as.vector(na.omit(boxhills[[1]]$design.list[[i]]$y.new)))
  x.tmp = c(boxhills[[1]]$design.list[[i]]$x, 
            as.vector(na.omit(boxhills[[1]]$design.list[[i]]$x.new)))
  PPH1_in1_boxhill[i] = getHypothesesPosteriors(
    prior.probs = prior_probs, 
    evidences = c(
      Evidence_gp(y.tmp, x.tmp, model0, nuggetBH),
      Evidence_gp(y.tmp, x.tmp, model1, nuggetBH)
    )
  )[2]
  # input 2
  y.tmp = c(boxhills[[2]]$design.list[[i]]$y, 
            as.vector(na.omit(boxhills[[2]]$design.list[[i]]$y.new)))
  x.tmp = c(boxhills[[2]]$design.list[[i]]$x, 
            as.vector(na.omit(boxhills[[2]]$design.list[[i]]$x.new)))
  PPH1_in2_boxhill[i] = getHypothesesPosteriors(
    prior.probs = prior_probs, 
    evidences = c(
      Evidence_gp(y.tmp, x.tmp, model0, nuggetBH),
      Evidence_gp(y.tmp, x.tmp, model1, nuggetBH)
    )
  )[2]
  # for input 3
  y.tmp = c(boxhills[[3]]$design.list[[i]]$y, 
            as.vector(na.omit(boxhills[[3]]$design.list[[i]]$y.new)))
  x.tmp = c(boxhills[[3]]$design.list[[i]]$x, 
            as.vector(na.omit(boxhills[[3]]$design.list[[i]]$x.new)))
  PPH1_in3_boxhill[i] = getHypothesesPosteriors(
    prior.probs = prior_probs, 
    evidences = c(
      Evidence_gp(y.tmp, x.tmp, model0, nuggetBH),
      Evidence_gp(y.tmp, x.tmp, model1, nuggetBH)
    )
  )[2]
  # input 4
  y.tmp = c(boxhills[[4]]$design.list[[i]]$y, 
            as.vector(na.omit(boxhills[[4]]$design.list[[i]]$y.new)))
  x.tmp = c(boxhills[[4]]$design.list[[i]]$x, 
            as.vector(na.omit(boxhills[[4]]$design.list[[i]]$x.new)))
  PPH1_in4_boxhill[i] = getHypothesesPosteriors(
    prior.probs = prior_probs, 
    evidences = c(
      Evidence_gp(y.tmp, x.tmp, model0, nuggetBH),
      Evidence_gp(y.tmp, x.tmp, model1, nuggetBH)
    )
  )[2]
}

# PPH1 for Space-filling
PPH1_in1_spacefilling = rep(NA, numSims)
PPH1_in2_spacefilling = rep(NA, numSims)
PPH1_in3_spacefilling = rep(NA, numSims)
PPH1_in4_spacefilling = rep(NA, numSims)
nugget = nuggetSM
for(i in 1:numSims){
  # input 1
  fn.vals.tmp = seqmeds[[1]]$function.values.list[ , i]
  y.tmp = c(fn.vals.tmp[seqmeds[[1]]$x0.idx], fn.vals.tmp[x_spacefill1_idx])
  x.tmp = c(seqmeds[[1]]$x0, x_spacefill1)
  PPH1_in1_spacefilling[i] = getHypothesesPosteriors(
    prior.probs = prior_probs, 
    evidences = c(
      Evidence_gp(y.tmp, x.tmp, model0, nugget),
      Evidence_gp(y.tmp, x.tmp, model1, nugget)
    )
  )[2]
  # input 2
  fn.vals.tmp = seqmeds[[2]]$function.values.list[ , i]
  y.tmp = c(fn.vals.tmp[seqmeds[[2]]$x0.idx], fn.vals.tmp[x_spacefill2_idx])
  x.tmp = c(seqmeds[[2]]$x0, x_spacefill2)
  PPH1_in2_spacefilling[i] = getHypothesesPosteriors(
    prior.probs = prior_probs, 
    evidences = c(
      Evidence_gp(y.tmp, x.tmp, model0, nugget),
      Evidence_gp(y.tmp, x.tmp, model1, nugget)
    )
  )[2]
  # input 3
  fn.vals.tmp = seqmeds[[3]]$function.values.list[ , i]
  y.tmp = c(fn.vals.tmp[seqmeds[[3]]$x0.idx], fn.vals.tmp[x_spacefill3_idx])
  x.tmp = c(seqmeds[[3]]$x0, x_spacefill3)
  PPH1_in3_spacefilling[i] = getHypothesesPosteriors(
    prior.probs = prior_probs, 
    evidences = c(
      Evidence_gp(y.tmp, x.tmp, model0, nugget),
      Evidence_gp(y.tmp, x.tmp, model1, nugget)
    )
  )[2]
  # input 4
  fn.vals.tmp = seqmeds[[4]]$function.values.list[ , i]
  y.tmp = c(fn.vals.tmp[seqmeds[[4]]$design.list[[i]]$x.idx], 
            fn.vals.tmp[x_spacefill4_idx])
  x.tmp = c(seqmeds[[4]]$design.list[[i]]$x, x_spacefill4)
  PPH1_in4_spacefilling[i] = getHypothesesPosteriors(
    prior.probs = prior_probs, 
    evidences = c(
      Evidence_gp(y.tmp, x.tmp, model0, nugget),
      Evidence_gp(y.tmp, x.tmp, model1, nugget)
    )
  )[2]
}

# PPH1 for Uniformly-distributed Random
PPH1_in1_random = rep(NA, numSims)
PPH1_in2_random = rep(NA, numSims)
PPH1_in3_random = rep(NA, numSims)
PPH1_in4_random = rep(NA, numSims)
for(i in 1:numSims){
  # first, get uniformly sampled random points
  x_uniform_idx = sample(1:numx, Nnew)
  x_uniform = x_seq[x_uniform_idx]
  # now calculate PPH1
  # input 1
  fn.vals.tmp = seqmeds[[1]]$function.values.list[ , i]
  y.tmp = c(fn.vals.tmp[seqmeds[[1]]$x0.idx], fn.vals.tmp[x_uniform_idx])
  x.tmp = c(seqmeds[[1]]$x0, x_uniform)
  PPH1_in1_random[i] = getHypothesesPosteriors(
    prior.probs = prior_probs, 
    evidences = c(
      Evidence_gp(y.tmp, x.tmp, model0, nugget),
      Evidence_gp(y.tmp, x.tmp, model1, nugget)
    )
  )[2]
  # input 2
  fn.vals.tmp = seqmeds[[2]]$function.values.list[ , i]
  y.tmp = c(fn.vals.tmp[seqmeds[[2]]$x0.idx], fn.vals.tmp[x_uniform_idx])
  x.tmp = c(seqmeds[[2]]$x0, x_uniform)
  PPH1_in2_random[i] = getHypothesesPosteriors(
    prior.probs = prior_probs, 
    evidences = c(
      Evidence_gp(y.tmp, x.tmp, model0, nugget),
      Evidence_gp(y.tmp, x.tmp, model1, nugget)
    )
  )[2]
  # input 3
  fn.vals.tmp = seqmeds[[3]]$function.values.list[ , i]
  y.tmp = c(fn.vals.tmp[seqmeds[[3]]$x0.idx], fn.vals.tmp[x_uniform_idx])
  x.tmp = c(seqmeds[[3]]$x0, x_uniform)
  PPH1_in3_random[i] = getHypothesesPosteriors(
    prior.probs = prior_probs, 
    evidences = c(
      Evidence_gp(y.tmp, x.tmp, model0, nugget),
      Evidence_gp(y.tmp, x.tmp, model1, nugget)
    )
  )[2]
  # input 4
  fn.vals.tmp = seqmeds[[4]]$function.values.list[ , i]
  y.tmp = c(fn.vals.tmp[seqmeds[[4]]$design.list[[i]]$x.idx], 
            fn.vals.tmp[x_uniform_idx])
  x.tmp = c(seqmeds[[4]]$design.list[[i]]$x, x_uniform)
  PPH1_in4_random[i] = getHypothesesPosteriors(
    prior.probs = prior_probs, 
    evidences = c(
      Evidence_gp(y.tmp, x.tmp, model0, nugget),
      Evidence_gp(y.tmp, x.tmp, model1, nugget)
    )
  )[2]
}

# get the medians
PPH1_in1 = c(
  SeqMED = median(PPH1_in1_seqmed), 
  BoxHill = median(PPH1_in1_boxhill), 
  SpaceFill = median(PPH1_in1_spacefilling), 
  Random = median(PPH1_in1_random)
)
PPH1_in2 = c(
  SeqMED = median(PPH1_in2_seqmed), 
  BoxHill = median(PPH1_in2_boxhill), 
  SpaceFill = median(PPH1_in2_spacefilling), 
  Random = median(PPH1_in2_random)
)
PPH1_in3 = c(
  SeqMED = median(PPH1_in3_seqmed), 
  BoxHill = median(PPH1_in3_boxhill), 
  SpaceFill = median(PPH1_in3_spacefilling), 
  Random = median(PPH1_in3_random)
)
PPH1_in4 = c(
  SeqMED = median(PPH1_in4_seqmed), 
  BoxHill = median(PPH1_in4_boxhill), 
  SpaceFill = median(PPH1_in4_spacefilling), 
  Random = median(PPH1_in4_random)
)

ggdata = data.table(
  Extrapolation = PPH1_in1,
  `Inc Spread` = PPH1_in2,
  `Even Coverage` = PPH1_in3,
  `Random` = PPH1_in4,
  Design = c("SeqMED", "BoxHill", "SpaceFill", "Random")
)
ggdata = melt(ggdata, id.vars = c("Design"))
PPH1.plt = ggplot(ggdata, aes(x = variable, y = value, group = Design, 
                              color = Design)) + 
  geom_point() + 
  geom_path(linetype = 2) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.text.x = element_text(size = 5)) +
  labs(y = "P(H1|X, Y)", x = "Initial Data")
PPH1.plt
# ggsave("mvp_medianpph1.pdf",
#        plot = last_plot(),
#        device = "pdf",
#        path = image_path,
#        scale = 1,
#        width = 4.5,
#        height = 2,
#        units = c("in")
# )

# in the same plot
ggdata1 = data.table(
  `Extrapolation` = RSS01_in1,
  `Inc Spread` = RSS01_in2,
  `Even Coverage` = RSS01_in3,
  `Random` = RSS01_in4,
  Design = c("SeqMED", "BoxHill", "SpaceFill", "Random")
)
ggdata1 = melt(ggdata1, id.vars = c("Design"))
ggdata1$Metric = "RSS0/RSS1"
ggdata2 = data.table(
  `Extrapolation` = PPH1_in1,
  `Inc Spread` = PPH1_in2,
  `Even Coverage` = PPH1_in3,
  `Random` = PPH1_in4,
  Design = c("SeqMED", "BoxHill", "SpaceFill", "Random")
)


ggdata2 = melt(ggdata2, id.vars = c("Design"))
ggdata2$Metric = "P(H1|X,Y)"
ggdata3 = rbind(ggdata1, ggdata2)
ggdata3$Metric = factor(ggdata3$Metric)
plt_gvm3 = ggplot(ggdata3, aes(x = variable, y = value, group = Design, 
                               color = Design)) + 
  facet_wrap(vars(Metric), scales = "free") + 
  geom_point() + 
  geom_path(linetype = 2) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.title = element_blank(), 
        axis.text.x = element_text(angle = 45, vjust = 0.5))
plt_gvm3
# ggsave("mvp_medianlogrss01pph1.pdf",
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
###########
############
#############
##############
###############
###############
##############
#############
############
###########
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
################################################################################
################################################################################
################################################################################
################################################################################
# LM VS Plots
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################

################################################################################
# Sources/Libraries
################################################################################
# directories
output_home = "run_designs/smmed/lmvs_sims"
functions_home = "functions"

# for seqmed design
source(paste(functions_home, "/SeqMEDvs.R", sep = ""))
source(paste(functions_home, "/SeqMEDvs_batch.R", sep = ""))
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

library(Matrix)
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
# simulation settings, shared for both scenarios (2 vs. 3 factors)
################################################################################

mu_full = c(0.5, 0.5, 0.5) #

# settings
xmin = -1
xmax = 1
sigmasq01 = 0.5
sigmasq = 0.5
numCandidates = 5000 #
xmin = -1
xmax = 1
p = 3
k = 4 * p
initN = 5
pfull = length(mu_full)

# sequential settings
numSims = 50
numSeq = 9
N_seq = 3
alpha_seq = 1

# hypotheses
indices0 = c(1, 2) #
indices1 = 1:length(mu_full)
mu0 = rep(0, length(indices0))
mu1 = rep(0, length(indices1))
V0 = diag(rep(sigmasq01,length(mu0)))
V1 = diag(rep(sigmasq01,length(mu1)))
f0 = function(x) mu0 %*% x[indices0]
f1 = function(x) mu1 %*% x[indices1]

################################################################################
# Scenario 1: 2 factors are significant
################################################################################

# settings
betaT = mu_full[indices0]
indicesT = indices0
fT = function(x) betaT %*% x[indicesT]
set.seed(12)
initD = matrix(runif(n = pfull * initN, min = xmin, max = xmax), nrow = initN, ncol = pfull)
inity = as.vector(simulateYvs(initD[ , indicesT], initN, betaT, sigmasq, 1, seed = 123))

# import H0 sims
numSimsH0 = numSims
simsH0 = list()
for(i in 1:numSimsH0){
  simsH0[[i]] = readRDS(
    paste0(
      output_home, 
      "/lmvs2v3H0_sim", 
      i, 
      ".rds"))
}
smmedvsH0 = simsH0[[1]]

Nnew = numSeq * N_seq
Ntot = Nnew + initN

################################################################################
# non-sequential designs
################################################################################

# random design
set.seed(12345)
x_random = matrix(runif(n = pfull * Nnew, min = xmin, max = xmax), 
                  nrow = Nnew, ncol = pfull)
y_random = simulateYvs(x_random[ , indicesT], Nnew, betaT, sigmasq, 1, seed = 12)
x_random = rbind(initD, x_random)
y_random = c(inity, y_random)

# doptimal design
dopt_pts = as.matrix(expand.grid(c(-1, 1), c(-1, 1), c(-1, 1)))
dopt_pts = dopt_pts[sample(1:dim(dopt_pts)[1], replace = FALSE), ]
x_doptimal = matrix(NA, nrow = Nnew, ncol = pfull)
for(i in 0:(Nnew - 1)){
  x_doptimal[ i + 1, ] = as.matrix(dopt_pts[ 1 + (i %% dim(dopt_pts)[1]), ])
}
y_doptimal = simulateYvs(x_doptimal[ , indicesT], Nnew, betaT, sigmasq, 1, seed = 12)
x_doptimal = rbind(initD, x_doptimal)
y_doptimal = c(inity, y_doptimal)

# 3 level factorial design
factorial_pts = as.matrix(expand.grid(c(-1, 0, 1), c(-1, 0, 1), c(-1, 0, 1)))
factorial_pts = factorial_pts[sample(1:dim(factorial_pts)[1], replace = FALSE), ]
x_factorial = factorial_pts
y_factorial = simulateYvs(x_factorial[ , indicesT], Nnew, betaT, sigmasq, 1, seed = 12)
x_factorial = rbind(initD, x_factorial)
y_factorial = c(inity, y_factorial)

################################################################################
# Metrics
################################################################################

# Posterior Probability of H0 (sequential)

numSims_pph = 50
models = list("H0" = list(mu0, V0, NULL, indices0),
              "H1" = list(mu1, V1, NULL, indices1))
#EPPHs for other designs
EPPHs_rand = calcEPPH(x_random, Ntot, betaT, NULL, models, sigmasq, numSims_pph, indicesT, seed = 123)
EPPHs_dopt = calcEPPH(x_doptimal, Ntot, betaT, NULL, models, sigmasq, numSims_pph, indicesT, seed = 123)
EPPHs_fact = calcEPPH(x_factorial, Ntot, betaT, NULL, models, sigmasq, numSims_pph, indicesT, seed = 123)
#EPPHs for smmed
EPPHs_smmed = calcEPPHseqdata(smmedvsH0$y, smmedvsH0$D, models, sigmasq, initN, numSeq, N_seq)
#EPPHs for smmed sims
idxlast = numSeq + 1
EPPH0_smmedsims = matrix(NA, idxlast, numSimsH0)
EPPH1_smmedsims = matrix(NA, idxlast, numSimsH0)
for(i in 1:numSimsH0){
  simH0.temp = simsH0[[i]]
  EPPH_temp = calcEPPHseqdata(simH0.temp$y, simH0.temp$D, models, sigmasq, initN, numSeq, N_seq)
  EPPH0_smmedsims[ , i] = EPPH_temp[1 , ]
  EPPH1_smmedsims[ , i] = EPPH_temp[2 , ]
}
EPPH0_smmedsims_mean = apply(EPPH0_smmedsims, 1, mean)
EPPH1_smmedsims_mean = apply(EPPH1_smmedsims, 1, mean)
EPPH0_smmedsims_med = apply(EPPH0_smmedsims, 1, median)
EPPH1_smmedsims_med = apply(EPPH1_smmedsims, 1, median)

ggdata0 = data.table(
  x = 1:idxlast, 
  Random = rep(EPPHs_rand[1], idxlast), 
  DOptimal = rep(EPPHs_dopt[1], idxlast), 
  Factorial3 = rep(EPPHs_fact[1], idxlast), 
  SeqMED = EPPH0_smmedsims_mean, 
  Hypothesis = rep("H0", idxlast)
)
ggdata1 = data.table(
  x = 1:idxlast, 
  Random = rep(EPPHs_rand[2], idxlast), 
  DOptimal = rep(EPPHs_dopt[2], idxlast), 
  Factorial3 = rep(EPPHs_fact[2], idxlast), 
  SeqMED = EPPH1_smmedsims_mean, 
  Hypothesis = rep("H1", idxlast)
)
ggdata = rbind(ggdata0, ggdata1)
ggdata.melted = melt(ggdata, id = c("x", "Hypothesis"), value.name = "epph", 
                     variable.name = "Design")
plt = ggplot(ggdata.melted, aes(x = x, y = epph, color = Design, linetype = Design)) +
  facet_wrap(vars(Hypothesis)) + 
  geom_path() + 
  scale_linetype_manual(values=c(rep("dashed", 3), "solid")) + 
  geom_point(data = ggdata.melted[x == 10], aes(x = x, y = epph)) + 
  theme_bw() + 
  theme(panel.grid.minor = element_blank()) + 
  labs(y = "", x = "Stages")
plt
# ggsave("h0_epph.pdf",
#        plot = last_plot(),
#        device = "pdf",
#        path = image_path,
#        scale = 1,
#        width = 4.5,
#        height = 2,
#        units = c("in")
# )

# Posterior Probability of H0 (sequential)

numSims_seqPPH = numSimsH0

rand_seqPPH0 = matrix(NA, nrow = numSeq + 1, ncol = numSims_seqPPH)
dopt_seqPPH0 = matrix(NA, nrow = numSeq + 1, ncol = numSims_seqPPH)
fact_seqPPH0 = matrix(NA, nrow = numSeq + 1, ncol = numSims_seqPPH)
smmed_seqPPH0 = matrix(NA, nrow = numSeq + 1, ncol = numSims_seqPPH)

rand_seqPPH1 = matrix(NA, nrow = numSeq + 1, ncol = numSims_seqPPH)
dopt_seqPPH1 = matrix(NA, nrow = numSeq + 1, ncol = numSims_seqPPH)
fact_seqPPH1 = matrix(NA, nrow = numSeq + 1, ncol = numSims_seqPPH)
smmed_seqPPH1 = matrix(NA, nrow = numSeq + 1, ncol = numSims_seqPPH)

for(j in 1:numSims_seqPPH){
  set.seed(j + 50)
  smmedsim = simsH0[[j]]
  initD = smmedsim$initD
  inity = smmedsim$y[1:initN]
  
  # random design
  x_random = matrix(runif(n = pfull * Nnew, min = xmin, max = xmax), 
                    nrow = Nnew, ncol = pfull)
  y_random = simulateYvs(x_random[ , indicesT], Nnew, betaT, sigmasq, 1, seed = 12)
  #
  x_random = rbind(initD, x_random)
  y_random = c(inity, y_random)
  
  # doptimal design
  dopt_pts = as.matrix(expand.grid(c(-1, 1), c(-1, 1), c(-1, 1)))
  dopt_pts = dopt_pts[sample(1:dim(dopt_pts)[1], replace = FALSE), ]
  x_doptimal = matrix(NA, nrow = Nnew, ncol = pfull)
  for(i in 0:(Nnew - 1)){
    x_doptimal[ i + 1, ] = as.matrix(dopt_pts[ 1 + (i %% dim(dopt_pts)[1]), ])
  }
  y_doptimal = simulateYvs(x_doptimal[ , indicesT], Nnew, betaT, sigmasq, 1, seed = 12)
  #
  x_doptimal = rbind(initD, x_doptimal)
  y_doptimal = c(inity, y_doptimal)
  
  # 3 level factorial design
  factorial_pts = as.matrix(expand.grid(c(-1, 0, 1), c(-1, 0, 1), c(-1, 0, 1)))
  factorial_pts = factorial_pts[sample(1:dim(factorial_pts)[1], replace = FALSE), ]
  ###
  x_factorial = factorial_pts
  ###
  y_factorial = simulateYvs(x_factorial[ , indicesT], Nnew, betaT, sigmasq, 1, seed = 12)
  #
  x_factorial = rbind(initD, x_factorial)
  y_factorial = c(inity, y_factorial)
  
  rand_seqPPH.temp = calcEPPHseqdata(y_random, x_random, models, sigmasq, initN, numSeq, N_seq)
  rand_seqPPH0[ , j] = rand_seqPPH.temp[1, ]
  rand_seqPPH1[ , j] = rand_seqPPH.temp[2, ]
  dopt_seqPPH.temp = calcEPPHseqdata(y_doptimal, x_doptimal, models, sigmasq, initN, numSeq, N_seq)
  dopt_seqPPH0[ , j] = dopt_seqPPH.temp[1, ]
  dopt_seqPPH1[ , j] = dopt_seqPPH.temp[2, ]
  fact_seqPPH.temp = calcEPPHseqdata(y_factorial, x_factorial, models, sigmasq, initN, numSeq, N_seq)
  fact_seqPPH0[ , j] = fact_seqPPH.temp[1, ]
  fact_seqPPH1[ , j] = fact_seqPPH.temp[2, ]
  smmed_seqPPH.temp = calcEPPHseqdata(smmedsim$y, smmedsim$D, models, sigmasq, initN, numSeq, N_seq)
  smmed_seqPPH0[ , j] = smmed_seqPPH.temp[1, ]
  smmed_seqPPH1[ , j] = smmed_seqPPH.temp[2, ]
}

# mean
rand_seqPPH0_mean = apply(rand_seqPPH0, 1, mean)
dopt_seqPPH0_mean = apply(dopt_seqPPH0, 1, mean)
fact_seqPPH0_mean = apply(fact_seqPPH0, 1, mean)
smmed_seqPPH0_mean = apply(smmed_seqPPH0, 1, mean)
rand_seqPPH1_mean = apply(rand_seqPPH1, 1, mean)
dopt_seqPPH1_mean = apply(dopt_seqPPH1, 1, mean)
fact_seqPPH1_mean = apply(fact_seqPPH1, 1, mean)
smmed_seqPPH1_mean = apply(smmed_seqPPH1, 1, mean)

# plot
ggdata0 = data.table(
  x = 1:idxlast, 
  Random = rand_seqPPH0_mean, 
  DOptimal = dopt_seqPPH0_mean, 
  Factorial3 = fact_seqPPH0_mean, 
  SeqMED = smmed_seqPPH0_mean, 
  Hypothesis = rep("H0", idxlast)
)
ggdata1 = data.table(
  x = 1:idxlast, 
  Random = rand_seqPPH1_mean, 
  DOptimal = dopt_seqPPH1_mean, 
  Factorial3 = fact_seqPPH1_mean, 
  SeqMED = smmed_seqPPH1_mean, 
  Hypothesis = rep("H1", idxlast)
)
ggdata = rbind(ggdata0, ggdata1)
ggdata.melted = melt(ggdata, id = c("x", "Hypothesis"), value.name = "epph", 
                     variable.name = "Design")
plt1 = ggplot(ggdata.melted, aes(x = x, y = epph, color = Design, linetype = Design)) +
  facet_wrap(vars(Hypothesis)) + 
  geom_path() + 
  scale_linetype_manual(values=c(rep("dashed", 3), "solid")) + 
  geom_point(data = ggdata.melted[x == 10], aes(x = x, y = epph)) + 
  theme_bw() + 
  theme(panel.grid.minor = element_blank()) + 
  labs(y = "", x = "Stages")
plt1
# ggsave("h0_epph_seq.pdf",
#        plot = last_plot(),
#        device = "pdf",
#        path = image_path,
#        scale = 1,
#        width = 4.5,
#        height = 2,
#        units = c("in")
# )

# MSE beta
smmedvs = smmedvsH0
hyp_mu = mu0
hyp_V = V0
hyp_ind = indices0
mseBn_smmed = getMSEBeta(smmedvs$D, Ntot, betaT, hyp_mu, hyp_V, sigmasq, NULL, hyp_ind)$MSE_postmean
mseBn_rand = getMSEBeta(x_random, Ntot, betaT, hyp_mu, hyp_V, sigmasq, NULL, hyp_ind)$MSE_postmean
mseBn_dopt = getMSEBeta(x_doptimal, Ntot, betaT, hyp_mu, hyp_V, sigmasq, NULL, hyp_ind)$MSE_postmean
mseBn_fact = getMSEBeta(x_factorial, Ntot, betaT, hyp_mu, hyp_V, sigmasq, NULL, hyp_ind)$MSE_postmean

# plot
b1 = c(mseBn_smmed[1], mseBn_rand[1], mseBn_dopt[1], mseBn_fact[1])
b2 = c(mseBn_smmed[2], mseBn_rand[2], mseBn_dopt[2], mseBn_fact[2])
ggdesigns = c("SeqMED", "Random", "Doptimal", "Factorial3")
ggdata = data.frame(Designs = factor(rep(ggdesigns, 2), 
                                     levels = ggdesigns[c(2, 1, 4, 3)]), 
                    MSE = c(b1, b2), beta = rep(c("B1", "B2"), each = length(b1)))
ggplot(ggdata, aes(x = Designs, y = MSE)) + 
  geom_bar(stat = "identity") +
  facet_wrap(vars(beta)) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.text.x = element_text(angle = 45, vjust = 0.5)) +
  labs(y = NULL)
# ggsave("h0_msebeta.pdf",
#        plot = last_plot(),
#        device = "pdf",
#        path = image_path,
#        scale = 1,
#        width = 4.5,
#        height = 2,
#        units = c("in")
# )

# the design

numBreaks = 8
smmed = smmedvsH0
Ntot = dim(smmed$D)[1]

maxcounts = rep(NA, length(mu_full))
for(i in 1:length(mu_full)){
  marginal = i
  h = hist(smmed$D[ (initN + 1):Ntot, marginal], breaks = numBreaks, plot = FALSE)
  maxcounts[i] = max(h$counts)
}

marginals = matrix(NA, nrow = length((initN + 1):Ntot), ncol = 3)
for(i in 1:(dim(marginals)[2])) {
  marginals[, i] = smmed$D[ (initN + 1):Ntot, i]
}
colnames(marginals) = paste("Variable", 1:3, sep = " ")
marginals = as.data.table(marginals)
marginals.tall = melt(marginals, measure.vars = 1:3)
ggplot(marginals.tall, aes(x = value)) + 
  facet_wrap(vars(variable)) +
  geom_histogram(binwidth = 0.12, closed = "right", aes(y = after_stat(density))) + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) + 
  labs(x = "x")
# ggsave("h0_marginals2.pdf",
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
# Scenario 2: all 3 factors are significant
################################################################################

betaT = mu_full[indices1]
indicesT = indices1
fT = function(x) betaT %*% x[indicesT]
set.seed(12)
initD = matrix(runif(n = pfull * initN, min = xmin, max = xmax), nrow = initN, ncol = pfull)
inity = as.vector(simulateYvs(initD[ , indicesT], initN, betaT, sigmasq, 1, seed = 123))

# import H1 sims
numSimsH1 = numSims
simsH1 = list()
for(i in 1:numSimsH0){
  simsH1[[i]] = readRDS(
    paste0(
      output_home, 
      "/lmvs2v3H1_sim", 
      i, 
      ".rds"))
}
smmedvsH1 = simsH1[[1]]

################################################################################
# non-sequential designs
################################################################################

Nnew = numSeq * N_seq
Ntot = Nnew + initN

# random design
set.seed(123)
x_random = matrix(runif(n = pfull * Nnew, min = xmin, max = xmax), 
                  nrow = Nnew, ncol = pfull)
y_random = simulateYvs(x_random[ , indicesT], Nnew, betaT, sigmasq, 1, seed = 12)
#
x_random = rbind(initD, x_random)
y_random = c(inity, y_random)

# postmean1 of random design
postmean0_random = matrix(NA, nrow = length(mu0), numSeq)
postmean1_random = matrix(NA, nrow = length(mu1), numSeq)
for(i in 1:numSeq){
  ind_seq = 1:(initN + i * N_seq)
  postmean1_random[ , i] = postmean(y_random[ind_seq], x_random[ind_seq , indices1], length(ind_seq), mu1, V1, sigmasq, NULL)
  postmean0_random[ , i] = postmean(y_random[ind_seq], x_random[ind_seq , indices0], length(ind_seq), mu0, V0, sigmasq, NULL)
  
}

# doptimal design
dopt_pts = as.matrix(expand.grid(c(-1, 1), c(-1, 1), c(-1, 1)))
dopt_pts = dopt_pts[sample(1:dim(dopt_pts)[1], replace = FALSE), ]
x_doptimal = matrix(NA, nrow = Nnew, ncol = pfull)
for(i in 0:(Nnew - 1)){
  x_doptimal[ i + 1, ] = as.matrix(dopt_pts[ 1 + (i %% dim(dopt_pts)[1]), ])
}
y_doptimal = simulateYvs(x_doptimal[ , indicesT], Nnew, betaT, sigmasq, 1, seed = 12)
#
x_doptimal = rbind(initD, x_doptimal)
y_doptimal = c(inity, y_doptimal)

# postmean1 of doptimal design
postmean0_doptimal = matrix(NA, nrow = length(mu0), numSeq)
postmean1_doptimal = matrix(NA, nrow = length(mu1), numSeq)
for(i in 1:numSeq){
  ind_seq = 1:(initN + i * N_seq)
  postmean1_doptimal[ , i] = postmean(y_doptimal[ind_seq], x_doptimal[ind_seq , indices1], length(ind_seq), mu1, V1, sigmasq, NULL)
  postmean0_doptimal[ , i] = postmean(y_doptimal[ind_seq], x_doptimal[ind_seq , indices0], length(ind_seq), mu0, V0, sigmasq, NULL)
}

# 3 level factorial design
factorial_pts = as.matrix(expand.grid(c(-1, 0, 1), c(-1, 0, 1), c(-1, 0, 1)))
factorial_pts = factorial_pts[sample(1:dim(factorial_pts)[1], replace = FALSE), ]
###
x_factorial = factorial_pts
###
y_factorial = simulateYvs(x_factorial[ , indicesT], Nnew, betaT, sigmasq, 1, seed = 12)
#
x_factorial = rbind(initD, x_factorial)
y_factorial = c(inity, y_factorial)

# postmean1 of doptimal design
postmean0_factorial = matrix(NA, nrow = length(mu0), numSeq)
postmean1_factorial = matrix(NA, nrow = length(mu1), numSeq)
for(i in 1:numSeq){
  ind_seq = 1:(initN + i * N_seq)
  postmean1_factorial[ , i] = postmean(y_factorial[ind_seq], x_factorial[ind_seq , indices1], length(ind_seq), mu1, V1, sigmasq, NULL)
  postmean0_factorial[ , i] = postmean(y_factorial[ind_seq], x_factorial[ind_seq , indices0], length(ind_seq), mu0, V0, sigmasq, NULL)
}



# -- EPPH non-sequential --- #

numSims_pph = 50
models = list("H0" = list(mu0, V0, NULL, indices0),
              "H1" = list(mu1, V1, NULL, indices1))
#EPPHs for other designs
EPPHs_rand = calcEPPH(x_random, Ntot, betaT, NULL, models, sigmasq, numSims_pph, indicesT, seed = 123)
EPPHs_dopt = calcEPPH(x_doptimal, Ntot, betaT, NULL, models, sigmasq, numSims_pph, indicesT, seed = 123)
EPPHs_fact = calcEPPH(x_factorial, Ntot, betaT, NULL, models, sigmasq, numSims_pph, indicesT, seed = 123)
#EPPHs for smmed
EPPHs_smmed = calcEPPHseqdata(smmedvsH1$y, smmedvsH1$D, models, sigmasq, initN, numSeq, N_seq)

#EPPHs for smmed sims
idxlast = numSeq + 1
EPPH0_smmedsims = matrix(NA, idxlast, numSimsH1)
EPPH1_smmedsims = matrix(NA, idxlast, numSimsH1)
for(i in 1:numSimsH1){
  simH1.temp = simsH1[[i]]
  EPPH_temp = calcEPPHseqdata(simH1.temp$y, simH1.temp$D, models, sigmasq, initN, numSeq, N_seq)
  EPPH0_smmedsims[ , i] = EPPH_temp[1 , ]
  EPPH1_smmedsims[ , i] = EPPH_temp[2 , ]
}
EPPH0_smmedsims_mean = apply(EPPH0_smmedsims, 1, mean)
EPPH1_smmedsims_mean = apply(EPPH1_smmedsims, 1, mean)

EPPH0_smmedsims_med = apply(EPPH0_smmedsims, 1, median)
EPPH1_smmedsims_med = apply(EPPH1_smmedsims, 1, median)

ggdata0 = data.table(
  x = 1:idxlast, 
  Random = rep(EPPHs_rand[1], idxlast), 
  DOptimal = rep(EPPHs_dopt[1], idxlast), 
  Factorial3 = rep(EPPHs_fact[1], idxlast), 
  SeqMED = EPPH0_smmedsims_mean, 
  Hypothesis = rep("H0", idxlast)
)
ggdata1 = data.table(
  x = 1:idxlast, 
  Random = rep(EPPHs_rand[2], idxlast), 
  DOptimal = rep(EPPHs_dopt[2], idxlast), 
  Factorial3 = rep(EPPHs_fact[2], idxlast), 
  SeqMED = EPPH1_smmedsims_mean, 
  Hypothesis = rep("H1", idxlast)
)
ggdata = rbind(ggdata0, ggdata1)
ggdata.melted = melt(ggdata, id = c("x", "Hypothesis"), value.name = "epph", 
                     variable.name = "Design")
plt = ggplot(ggdata.melted, aes(x = x, y = epph, color = Design, linetype = Design)) +
  facet_wrap(vars(Hypothesis)) + 
  geom_path() + 
  scale_linetype_manual(values=c(rep("dashed", 3), "solid")) + 
  geom_point(data = ggdata.melted[x == 10], aes(x = x, y = epph)) + 
  theme_bw() + 
  theme(panel.grid.minor = element_blank()) + 
  labs(y = "", x = "Stages")
plt
# ggsave("h1_epph.pdf",
#        plot = last_plot(),
#        device = "pdf",
#        path = image_path,
#        scale = 1,
#        width = 4.5,
#        height = 2,
#        units = c("in")
# )


# --- EPPH sequential --- #

numSims_seqPPH = numSimsH1

idxlast = numSeq + 1
rand_seqPPH0 = matrix(NA, nrow = idxlast, ncol = numSims_seqPPH)
dopt_seqPPH0 = matrix(NA, nrow = idxlast, ncol = numSims_seqPPH)
fact_seqPPH0 = matrix(NA, nrow = idxlast, ncol = numSims_seqPPH)
smmed_seqPPH0 = matrix(NA, nrow = idxlast, ncol = numSims_seqPPH)

rand_seqPPH1 = matrix(NA, nrow = idxlast, ncol = numSims_seqPPH)
dopt_seqPPH1 = matrix(NA, nrow = idxlast, ncol = numSims_seqPPH)
fact_seqPPH1 = matrix(NA, nrow = idxlast, ncol = numSims_seqPPH)
smmed_seqPPH1 = matrix(NA, nrow = idxlast, ncol = numSims_seqPPH)

for(j in 1:numSims_seqPPH){
  set.seed(j + 50)
  smmedsim = simsH1[[j]]
  initD = smmedsim$initD
  inity = smmedsim$y[1:initN]
  
  # random design
  x_random = matrix(runif(n = pfull * Nnew, min = xmin, max = xmax), 
                    nrow = Nnew, ncol = pfull)
  y_random = simulateYvs(x_random[ , indicesT], Nnew, betaT, sigmasq, 1, seed = 12)
  #
  x_random = rbind(initD, x_random)
  y_random = c(inity, y_random)
  
  # doptimal design
  dopt_pts = as.matrix(expand.grid(c(-1, 1), c(-1, 1), c(-1, 1)))
  dopt_pts = dopt_pts[sample(1:dim(dopt_pts)[1], replace = FALSE), ]
  x_doptimal = matrix(NA, nrow = Nnew, ncol = pfull)
  for(i in 0:(Nnew - 1)){
    x_doptimal[ i + 1, ] = as.matrix(dopt_pts[ 1 + (i %% dim(dopt_pts)[1]), ])
  }
  y_doptimal = simulateYvs(x_doptimal[ , indicesT], Nnew, betaT, sigmasq, 1, seed = 12)
  #
  x_doptimal = rbind(initD, x_doptimal)
  y_doptimal = c(inity, y_doptimal)
  
  # 3 level factorial design
  factorial_pts = as.matrix(expand.grid(c(-1, 0, 1), c(-1, 0, 1), c(-1, 0, 1)))
  factorial_pts = factorial_pts[sample(1:dim(factorial_pts)[1], replace = FALSE), ]
  ###
  x_factorial = factorial_pts
  ###
  y_factorial = simulateYvs(x_factorial[ , indicesT], Nnew, betaT, sigmasq, 1, seed = 12)
  #
  x_factorial = rbind(initD, x_factorial)
  y_factorial = c(inity, y_factorial)
  
  rand_seqPPH.temp = calcEPPHseqdata(y_random, x_random, models, sigmasq, initN, numSeq, N_seq)
  rand_seqPPH0[ , j] = rand_seqPPH.temp[1, ]
  rand_seqPPH1[ , j] = rand_seqPPH.temp[2, ]
  dopt_seqPPH.temp = calcEPPHseqdata(y_doptimal, x_doptimal, models, sigmasq, initN, numSeq, N_seq)
  dopt_seqPPH0[ , j] = dopt_seqPPH.temp[1, ]
  dopt_seqPPH1[ , j] = dopt_seqPPH.temp[2, ]
  fact_seqPPH.temp = calcEPPHseqdata(y_factorial, x_factorial, models, sigmasq, initN, numSeq, N_seq)
  fact_seqPPH0[ , j] = fact_seqPPH.temp[1, ]
  fact_seqPPH1[ , j] = fact_seqPPH.temp[2, ]
  smmed_seqPPH.temp = calcEPPHseqdata(smmedsim$y, smmedsim$D, models, sigmasq, initN, numSeq, N_seq)
  smmed_seqPPH0[ , j] = smmed_seqPPH.temp[1, ]
  smmed_seqPPH1[ , j] = smmed_seqPPH.temp[2, ]
}

# mean
rand_seqPPH0_mean = apply(rand_seqPPH0, 1, mean)
dopt_seqPPH0_mean = apply(dopt_seqPPH0, 1, mean)
fact_seqPPH0_mean = apply(fact_seqPPH0, 1, mean)
smmed_seqPPH0_mean = apply(smmed_seqPPH0, 1, mean)
rand_seqPPH1_mean = apply(rand_seqPPH1, 1, mean)
dopt_seqPPH1_mean = apply(dopt_seqPPH1, 1, mean)
fact_seqPPH1_mean = apply(fact_seqPPH1, 1, mean)
smmed_seqPPH1_mean = apply(smmed_seqPPH1, 1, mean)

#plot
ggdata0 = data.table(
  x = 1:idxlast, 
  Random = rand_seqPPH0_mean, 
  DOptimal = dopt_seqPPH0_mean, 
  Factorial3 = fact_seqPPH0_mean, 
  SeqMED = smmed_seqPPH0_mean, 
  Hypothesis = rep("H0", idxlast)
)
ggdata1 = data.table(
  x = 1:idxlast, 
  Random = rand_seqPPH1_mean, 
  DOptimal = dopt_seqPPH1_mean, 
  Factorial3 = fact_seqPPH1_mean, 
  SeqMED = smmed_seqPPH1_mean, 
  Hypothesis = rep("H1", idxlast)
)
ggdata = rbind(ggdata0, ggdata1)
ggdata.melted = melt(ggdata, id = c("x", "Hypothesis"), value.name = "epph", 
                     variable.name = "Design")
plt1 = ggplot(ggdata.melted, aes(x = x, y = epph, color = Design, linetype = Design)) +
  facet_wrap(vars(Hypothesis)) + 
  geom_path() + 
  scale_linetype_manual(values=c(rep("dashed", 3), "solid")) + 
  geom_point(data = ggdata.melted[x == 10], aes(x = x, y = epph)) + 
  theme_bw() + 
  theme(panel.grid.minor = element_blank()) + 
  labs(y = "", x = "Stages")
plt1
# ggsave("h1_epph_seq.pdf",
#        plot = last_plot(),
#        device = "pdf",
#        path = image_path,
#        scale = 1,
#        width = 4.5,
#        height = 2,
#        units = c("in")
# )

# MSE Beta

smmedvs = smmedvsH1
hyp_mu = mu1
hyp_V = V1
hyp_ind = indices1

mseBn_smmed = getMSEBeta(smmedvs$D, Ntot, betaT, hyp_mu, hyp_V, sigmasq, NULL, hyp_ind)$MSE_postmean
mseBn_rand = getMSEBeta(x_random, Ntot, betaT, hyp_mu, hyp_V, sigmasq, NULL, hyp_ind)$MSE_postmean
mseBn_dopt = getMSEBeta(x_doptimal, Ntot, betaT, hyp_mu, hyp_V, sigmasq, NULL, hyp_ind)$MSE_postmean
mseBn_fact = getMSEBeta(x_factorial, Ntot, betaT, hyp_mu, hyp_V, sigmasq, NULL, hyp_ind)$MSE_postmean

b1 = c(mseBn_smmed[1], mseBn_rand[1], mseBn_dopt[1], mseBn_fact[1])
b2 = c(mseBn_smmed[2], mseBn_rand[2], mseBn_dopt[2], mseBn_fact[2])
b3 = c(mseBn_smmed[3], mseBn_rand[3], mseBn_dopt[3], mseBn_fact[3])

ggdesigns = c("SeqMED", "Random", "Doptimal", "Factorial3")
ggdata = data.frame(Designs = factor(rep(ggdesigns, 3), 
                                     levels = ggdesigns[c(2, 1, 4, 3)]), 
                    MSE = c(b1, b2, b3), beta = rep(c("B1", "B2", "B3"), each = length(b1)))
ggplot(ggdata, aes(x = Designs, y = MSE)) + 
  geom_bar(stat = "identity") +
  facet_wrap(vars(beta)) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.text.x = element_text(angle = 45, vjust = 0.5)) +
  labs(y = NULL)
# ggsave("h1_msebeta.pdf",
#        plot = last_plot(),
#        device = "pdf",
#        path = image_path,
#        scale = 1,
#        width = 4.5,
#        height = 2,
#        units = c("in")
# )


# --- the design --- #

# smmed = smmedvsH1
# Ntot = dim(smmed$D)[1]
# 
# marginals = matrix(NA, nrow = length((initN + 1):Ntot), ncol = 3)
# for(i in 1:(dim(marginals)[2])) {
#   marginals[, i] = smmed$D[ (initN + 1):Ntot, i]
# }
# colnames(marginals) = paste("Marginal", 1:3, sep = " ")
# marginals = as.data.table(marginals)
# marginals.tall = melt(marginals, measure.vars = 1:3)
# ggplot(marginals.tall, aes(x = value)) + 
#   facet_wrap(vars(variable)) +
#   geom_histogram(binwidth = 0.12, closed = "right", aes(y = after_stat(density))) + 
#   theme_bw() + 
#   theme(panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank()) + 
#   labs(x = "x")
# ggsave("h1_marginals.pdf",
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
###########
############
#############
##############
###############
###############
##############
#############
############
###########
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
################################################################################
################################################################################
################################################################################
################################################################################
# GP VS Plots
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################

################################################################################
# Sources/Libraries
################################################################################
# tamu cluster
output_home = "run_designs/gp/gpvs_sims"

# --- Sources/Libraries --- #
source(paste(functions_home, "/SeqMEDvs.R", sep = ""))
source(paste(functions_home, "/SeqMEDvs_batch.R", sep = ""))

# --- Helper Functions for Evaluation Metrics --- #
add_errorbands = function(xs, ys, MoE, color){
  y_lower = ys - MoE
  y_upper = ys + MoE
  polygon(c(xs,rev(xs)),c(y_lower,rev(y_upper)),col=color, border = NA)
}

getAICBIC = function(model = NULL, p, y, X, beta, sigmasq){
  # if(sum(class(model) == "glmnet") > 0){ # if glm
  #   tLL = model$nulldev - deviance(model)
  #   k = model$df
  #   n = model$nobs
  #   aic = -tLL + 2* k + 2* k * (k + 1) / (n - k - 1)
  #   bic = -tLL + log(n) * k
  # }
  # if(class(model) == "lm"){
  if(length(y) != dim(X)[1]) stop("length of y is not equal to number of rows in X")
  n = length(y)
  eval_criteria = function(crit_const){
    - 2 * dmvnorm(y, mean = X %*% beta, sigma = diag(rep(sigmasq, length(y)))) + crit_const * p
  }
  aic = eval_criteria(2)
  bic = eval_criteria(log(n))
  # }
  return(c("aic" = aic, "bic" = bic))
}

logjointlik = function(x, y, type_arg, l_arg, nugget = NULL){
  joint_var = getCov(x, x, type_arg, l_arg)
  if(!is.null(nugget)) joint_var = joint_var + diag(rep(nugget, length(y)))
  return(dmvnorm(y, mean = rep(0, length(y)), sigma = joint_var, log = TRUE))
}

logjointlik_star = function(x_star, y_star, x_input, y_input, type_arg, l_arg, nugget = NULL){
  y = c(y_input, y_star)
  x = rbind(x_input, x_star)
  logjointlik(x, y, type_arg, l_arg, nugget)
}

calcPPH = function(x, y, indices0, indices1, type_args, l_args, nugget){
  initD0.temp = x[ , indices0, drop = FALSE]
  initD1.temp = x[ , indices1, drop = FALSE]
  loglikH0 = logjointlik(initD0.temp, y, type_args[1], l_args[1], nugget)
  loglikH1 = logjointlik(initD1.temp, y, type_args[2], l_args[2], nugget)
  PPH0 = exp(loglikH0) / (exp(loglikH0) + exp(loglikH1))
  PPH1 = exp(loglikH1) / (exp(loglikH0) + exp(loglikH1))
  return(c("PPH0" = PPH0, "PPH1" = PPH1))
}

calcPPH_star = function(candidates, true_y, design_indices, inputx, inputy, 
                        indices0, indices1, hyp_type, hyp_l, nugget_term){
  loglikH0 = logjointlik_star(x_star = candidates[ design_indices, indices0, drop = FALSE], 
                              y_star = true_y[ design_indices], 
                              inputx[ , indices0, drop = FALSE], inputy, 
                              hyp_type[1], hyp_l[1], nugget = nugget_term)
  loglikH1 = logjointlik_star(x_star = candidates[ design_indices, indices1, drop = FALSE], 
                              y_star = true_y[ design_indices], 
                              inputx[ , indices1, drop = FALSE], inputy, 
                              hyp_type[2], hyp_l[2], nugget = nugget_term)
  PPH0 = exp(loglikH0) / (exp(loglikH0) + exp(loglikH1))
  PPH1 = exp(loglikH1) / (exp(loglikH0) + exp(loglikH1))
  return(c("PPH0" = PPH0, "PPH1" = PPH1))
}

################################################################################
# simulation settings, shared for both scenarios (2 vs. 3 factors)
################################################################################

numSims = 25

# settings
xmin = 0; xmax = 1
l01= c(0.1, 0.1)
type01 = c(1, 1)

N = 3
p = 2
k = 4 * p
alpha = 1
indices0 = c(1)
indices1 = c(1, 2)
N2 = 9
nugget = 1e-1

# read in grid info
grid12 = readRDS(paste(output_home, "/grid_gpvs12.rds", sep = ""))
# the grid of inputs
numx = grid12[[1]]
x_seq = grid12[[2]]
x_grid = as.matrix(grid12[[3]])
null_cov1d = grid12[[4]]
null_mean1d = grid12[[5]]
null_cov2d = grid12[[6]]
null_mean2d = grid12[[7]]

# input points (randomly selected)
set.seed(2)
# set.seed(6)
x_input_ind = sample(1:dim(x_grid)[1], N)
x_input = x_grid[x_input_ind, ]

order_x = order(x_grid[ , 1])

# read in simulations
gpvs_sims = readRDS(paste0(output_home, "/gpvs_sims/gpvs_sims.rds"))
mmedgpf1d_sims = gpvs_sims$mmedgp_f1dsims
smmedgpf1d_sims = gpvs_sims$smmedgp_f1dsims
mmedgpf2d_sims = gpvs_sims$mmedgp_f2dsims
smmedgpf2d_sims = gpvs_sims$smmedgp_f2dsims

x_rand_ind_mat = matrix(NA, nrow = N2, ncol = numSims)
set.seed(1)
for(j in 1:numSims) x_rand_ind_mat[ , j] = sample(1:dim(x_grid)[1], N2)

################################################################################
# Other designs
################################################################################

# random points
set.seed(1990)
# set.seed(1997)
x_rand_ind = sample(1:dim(x_grid)[1], N2)
x_rand = x_grid[x_rand_ind, ]

# points where x2 = 1
gridlen = length(x_seq)
x_at1_ind = (1:N2) * floor(gridlen / N2) + gridlen * (gridlen - 1)
x_at1 = x_grid[x_at1_ind, ]

# points on diagonal
x_diag_ind = (1:N2) * floor(gridlen / N2) + gridlen * ((1:N2) * floor(gridlen / N2) - 1)
x_diag = x_grid[x_diag_ind, ]

# grid, only works when have 9 pts
beg_mid_end = c(1, ceiling(gridlen/(sqrt(9) - 1)) * 1:(sqrt(9) - 2), gridlen)
x_sf_ind = c(beg_mid_end, beg_mid_end + floor((gridlen * (gridlen - 1))) / 2, beg_mid_end + gridlen * (gridlen - 1))
x_sf = x_grid[x_sf_ind, ]

################################################################################
# Scenario 1: 1 dimension
################################################################################

smmed_PPH0_mat = matrix(NA, nrow = N2 + 1, ncol = numSims)
mmed_PPH0_mat = matrix(NA, nrow = N2 + 1, ncol = numSims)
rand_PPH0_mat = matrix(NA, nrow = N2 + 1, ncol = numSims)
x_at1_PPH0_mat = matrix(NA, nrow = N2 + 1, ncol = numSims)
x_diag_PPH0_mat = matrix(NA, nrow = N2 + 1, ncol = numSims)
x_sf_PPH0_mat = matrix(NA, nrow = N2 + 1, ncol = numSims)

smmed_PPH1_mat = matrix(NA, nrow = N2 + 1, ncol = numSims)
mmed_PPH1_mat = matrix(NA, nrow = N2 + 1, ncol = numSims)
rand_PPH1_mat = matrix(NA, nrow = N2 + 1, ncol = numSims)
x_at1_PPH1_mat = matrix(NA, nrow = N2 + 1, ncol = numSims)
x_diag_PPH1_mat = matrix(NA, nrow = N2 + 1, ncol = numSims)
x_sf_PPH1_mat = matrix(NA, nrow = N2 + 1, ncol = numSims)

# for each new point in each design
set.seed(1)
for(j in 1:numSims){
  mmed_gp_vs = mmedgpf1d_sims[[j]]
  smmed_gp_vs = smmedgpf1d_sims[[j]]
  x_rand_ind_mat[ , j] = sample(1:dim(x_grid)[1], N2)
  x_rand_ind = x_rand_ind_mat[ , j]
  
  x_input_ind = gpvs_sims$x_input_ind[ , j]
  x_input = x_grid[x_input_ind, ]
  y_seq = gpvs_sims$y_seq_mat1d[ , j]
  y_input = y_seq[x_input_ind]
  
  # from initial data
  initPPH = calcPPH(x_input, y_input, indices0, indices1, type01, l01, nugget)
  
  smmed_PPH0_mat[1, j] = initPPH[1]
  mmed_PPH0_mat[1, j] = initPPH[1]
  rand_PPH0_mat[1, j] = initPPH[1]
  x_at1_PPH0_mat[1, j] = initPPH[1]
  x_diag_PPH0_mat[1, j] = initPPH[1]
  x_sf_PPH0_mat[1, j] = initPPH[1]
  
  smmed_PPH1_mat[1, j] = initPPH[2]
  mmed_PPH1_mat[1, j] = initPPH[2]
  rand_PPH1_mat[1, j] = initPPH[2]
  x_at1_PPH1_mat[1, j] = initPPH[2]
  x_diag_PPH1_mat[1, j] = initPPH[2]
  x_sf_PPH1_mat[1, j] = initPPH[2]
  
  for(i in 1:N2){
    which_ind = 1:i
    
    # smmed postprobs
    x_ind = smmed_gp_vs$indices[which_ind]
    smmedgpPPH.temp = calcPPH_star(x_grid, y_seq, x_ind, x_input, y_input,
                                   indices0, indices1, type01, l01, nugget)
    smmed_PPH0_mat[i + 1, j] = smmedgpPPH.temp[1]
    smmed_PPH1_mat[i + 1, j] = smmedgpPPH.temp[2]
    # mmed postprobs
    x_ind = mmed_gp_vs$indices[which_ind]
    mmedgpPPH.temp = calcPPH_star(x_grid, y_seq, x_ind, x_input, y_input,
                                  indices0, indices1, type01, l01, nugget)
    mmed_PPH0_mat[i + 1, j] = mmedgpPPH.temp[1]
    mmed_PPH1_mat[i + 1, j] = mmedgpPPH.temp[2]
    # random postprobs
    x_ind = x_rand_ind[which_ind]
    randPPH.temp = calcPPH_star(x_grid, y_seq, x_ind, x_input, y_input,
                                indices0, indices1, type01, l01, nugget)
    rand_PPH0_mat[i + 1, j] = randPPH.temp[1]
    rand_PPH1_mat[i + 1, j] = randPPH.temp[2]
    
    # x_at1 postprobs
    x_ind = x_at1_ind[which_ind]
    at1PPH.temp = calcPPH_star(x_grid, y_seq, x_ind, x_input, y_input,
                               indices0, indices1, type01, l01, nugget)
    x_at1_PPH0_mat[i + 1, j] = at1PPH.temp[1]
    x_at1_PPH1_mat[i + 1, j] = at1PPH.temp[2]
    
    # x_diag postprobs
    x_ind = x_diag_ind[which_ind]
    diagPPH.temp = calcPPH_star(x_grid, y_seq, x_ind, x_input, y_input,
                                indices0, indices1, type01, l01, nugget)
    x_diag_PPH0_mat[i + 1, j] = diagPPH.temp[1]
    x_diag_PPH1_mat[i + 1, j] = diagPPH.temp[2]
    
    # x_sf postprobs
    x_ind = x_sf_ind[which_ind]
    sfPPH.temp = calcPPH_star(x_grid, y_seq, x_ind, x_input, y_input,
                              indices0, indices1, type01, l01, nugget)
    x_sf_PPH0_mat[i + 1, j] = sfPPH.temp[1]
    x_sf_PPH1_mat[i + 1, j] = sfPPH.temp[2]
  }
}

par(mfrow = c(1, 2))

# by mean #
smmed_PPH0_avg = apply(smmed_PPH0_mat, 1, mean, na.rm = TRUE)
mmed_PPH0_avg = apply(mmed_PPH0_mat, 1, mean, na.rm = TRUE)
rand_PPH0_avg = apply(rand_PPH0_mat, 1, mean, na.rm = TRUE)
x_at1_PPH0_avg = apply(x_at1_PPH0_mat, 1, mean, na.rm = TRUE)
x_diag_PPH0_avg = apply(x_diag_PPH0_mat, 1, mean, na.rm = TRUE)
x_sf_PPH0_avg = apply(x_sf_PPH0_mat, 1, mean, na.rm = TRUE)

smmed_PPH1_avg = apply(smmed_PPH1_mat, 1, mean, na.rm = TRUE)
mmed_PPH1_avg = apply(mmed_PPH1_mat, 1, mean, na.rm = TRUE)
rand_PPH1_avg = apply(rand_PPH1_mat, 1, mean, na.rm = TRUE)
x_at1_PPH1_avg = apply(x_at1_PPH1_mat, 1, mean, na.rm = TRUE)
x_diag_PPH1_avg = apply(x_diag_PPH1_mat, 1, mean, na.rm = TRUE)
x_sf_PPH1_avg = apply(x_sf_PPH1_mat, 1, mean, na.rm = TRUE)

# by median #
smmed_PPH0_median = apply(smmed_PPH0_mat, 1, median, na.rm = TRUE)
mmed_PPH0_median = apply(mmed_PPH0_mat, 1, median, na.rm = TRUE)
rand_PPH0_median = apply(rand_PPH0_mat, 1, median, na.rm = TRUE)
x_at1_PPH0_median = apply(x_at1_PPH0_mat, 1, median, na.rm = TRUE)
x_diag_PPH0_median = apply(x_diag_PPH0_mat, 1, median, na.rm = TRUE)
x_sf_PPH0_median = apply(x_sf_PPH0_mat, 1, median, na.rm = TRUE)

smmed_PPH1_median = apply(smmed_PPH1_mat, 1, median, na.rm = TRUE)
mmed_PPH1_median = apply(mmed_PPH1_mat, 1, median, na.rm = TRUE)
rand_PPH1_median = apply(rand_PPH1_mat, 1, median, na.rm = TRUE)
x_at1_PPH1_median = apply(x_at1_PPH1_mat, 1, median, na.rm = TRUE)
x_diag_PPH1_median = apply(x_diag_PPH1_mat, 1, median, na.rm = TRUE)
x_sf_PPH1_median = apply(x_sf_PPH1_mat, 1, median, na.rm = TRUE)

idxlast = N2 + 1
ggdata0 = data.table(
  x = 1:(N2 + 1), 
  Random = rand_PPH0_avg, 
  `X2=1` = x_at1_PPH0_avg, 
  Diagonal = x_diag_PPH0_avg, 
  SpaceFilling = x_sf_PPH0_avg,
  SeqMED = smmed_PPH0_avg, 
  Hypothesis = rep("H0", idxlast)
)
ggdata1 = data.table(
  x = 1:(N2 + 1), 
  Random = rand_PPH1_avg, 
  `X2=1` = x_at1_PPH1_avg, 
  Diagonal = x_diag_PPH1_avg, 
  SpaceFilling = x_sf_PPH1_avg,
  SeqMED = smmed_PPH1_avg, 
  Hypothesis = rep("H1", idxlast)
)
ggdata = rbind(ggdata0, ggdata1)
ggdata.melted = melt(ggdata, id = c("x", "Hypothesis"), value.name = "epph", 
                     variable.name = "Design")
ggdata.gpvs.H0true = ggdata.melted
plt = ggplot(ggdata.melted, aes(x = x, y = epph, color = Design, linetype = Design)) +
  facet_wrap(vars(Hypothesis)) + 
  geom_path() + 
  scale_linetype_manual(values=c(rep("dashed", 4), "solid")) + 
  geom_point(data = ggdata.melted[x == 10], aes(x = x, y = epph)) + 
  theme_bw() + 
  scale_x_continuous(breaks = 1:10) +
  theme(panel.grid.minor = element_blank(), axis.text.x = element_text(size = 5)) + 
  labs(y = "", x = "Stages") + 
  ylim(0, 1)
plt
# ggsave("gpvsh0_epph.pdf",
#        plot = last_plot(),
#        device = "pdf",
#        path = image_path,
#        scale = 1,
#        width = 4.5,
#        height = 2,
#        units = c("in")
# )

################################################################################
# Scenario 2: 2 dimensions
################################################################################

smmed_PPH0_mat = matrix(NA, nrow = N2 + 1, ncol = numSims)
mmed_PPH0_mat = matrix(NA, nrow = N2 + 1, ncol = numSims)
rand_PPH0_mat = matrix(NA, nrow = N2 + 1, ncol = numSims)
x_at1_PPH0_mat = matrix(NA, nrow = N2 + 1, ncol = numSims)
x_diag_PPH0_mat = matrix(NA, nrow = N2 + 1, ncol = numSims)
x_sf_PPH0_mat = matrix(NA, nrow = N2 + 1, ncol = numSims)

smmed_PPH1_mat = matrix(NA, nrow = N2 + 1, ncol = numSims)
mmed_PPH1_mat = matrix(NA, nrow = N2 + 1, ncol = numSims)
rand_PPH1_mat = matrix(NA, nrow = N2 + 1, ncol = numSims)
x_at1_PPH1_mat = matrix(NA, nrow = N2 + 1, ncol = numSims)
x_diag_PPH1_mat = matrix(NA, nrow = N2 + 1, ncol = numSims)
x_sf_PPH1_mat = matrix(NA, nrow = N2 + 1, ncol = numSims)

# for each new point in each design
set.seed(1)
for(j in 1:numSims){
  mmed_gp_vs = mmedgpf2d_sims[[j]]
  smmed_gp_vs = smmedgpf2d_sims[[j]]
  x_rand_ind = x_rand_ind_mat[ , j] # created in the beginning
  
  x_input_ind = gpvs_sims$x_input_ind[ , j]
  x_input = x_grid[x_input_ind, ]
  y_seq = gpvs_sims$y_seq_mat2d[ , j]
  y_input = y_seq[x_input_ind]
  
  # from initial data
  initPPH = calcPPH(x_input, y_input, indices0, indices1, type01, l01, nugget)
  
  smmed_PPH0_mat[1, j] = initPPH[1]
  mmed_PPH0_mat[1, j] = initPPH[1]
  rand_PPH0_mat[1, j] = initPPH[1]
  x_at1_PPH0_mat[1, j] = initPPH[1]
  x_diag_PPH0_mat[1, j] = initPPH[1]
  x_sf_PPH0_mat[1, j] = initPPH[1]
  
  smmed_PPH1_mat[1, j] = initPPH[2]
  mmed_PPH1_mat[1, j] = initPPH[2]
  rand_PPH1_mat[1, j] = initPPH[2]
  x_at1_PPH1_mat[1, j] = initPPH[2]
  x_diag_PPH1_mat[1, j] = initPPH[2]
  x_sf_PPH1_mat[1, j] = initPPH[2]
  
  for(i in 1:N2){
    which_ind = 1:i
    
    # smmed postprobs
    x_ind = smmed_gp_vs$indices[which_ind]
    smmedgpPPH.temp = calcPPH_star(x_grid, y_seq, x_ind, x_input, y_input,
                                   indices0, indices1, type01, l01, nugget)
    smmed_PPH0_mat[i + 1, j] = smmedgpPPH.temp[1]
    smmed_PPH1_mat[i + 1, j] = smmedgpPPH.temp[2]
    
    # mmed postprobs
    x_ind = mmed_gp_vs$indices[which_ind]
    mmedgpPPH.temp = calcPPH_star(x_grid, y_seq, x_ind, x_input, y_input,
                                  indices0, indices1, type01, l01, nugget)
    mmed_PPH0_mat[i + 1, j] = mmedgpPPH.temp[1]
    mmed_PPH1_mat[i + 1, j] = mmedgpPPH.temp[2]
    
    # random postprobs
    x_ind = x_rand_ind[which_ind]
    randPPH.temp = calcPPH_star(x_grid, y_seq, x_ind, x_input, y_input,
                                indices0, indices1, type01, l01, nugget)
    rand_PPH0_mat[i + 1, j] = randPPH.temp[1]
    rand_PPH1_mat[i + 1, j] = randPPH.temp[2]
    
    # x_at1 postprobs
    x_ind = x_at1_ind[which_ind]
    at1PPH.temp = calcPPH_star(x_grid, y_seq, x_ind, x_input, y_input,
                               indices0, indices1, type01, l01, nugget)
    x_at1_PPH0_mat[i + 1, j] = at1PPH.temp[1]
    x_at1_PPH1_mat[i + 1, j] = at1PPH.temp[2]
    
    # x_diag postprobs
    x_ind = x_diag_ind[which_ind]
    diagPPH.temp = calcPPH_star(x_grid, y_seq, x_ind, x_input, y_input,
                                indices0, indices1, type01, l01, nugget)
    x_diag_PPH0_mat[i + 1, j] = diagPPH.temp[1]
    x_diag_PPH1_mat[i + 1, j] = diagPPH.temp[2]
    
    # x_sf postprobs
    x_ind = x_sf_ind[which_ind]
    sfPPH.temp = calcPPH_star(x_grid, y_seq, x_ind, x_input, y_input,
                              indices0, indices1, type01, l01, nugget)
    x_sf_PPH0_mat[i + 1, j] = sfPPH.temp[1]
    x_sf_PPH1_mat[i + 1, j] = sfPPH.temp[2]
  }
}

par(mfrow = c(1, 2))

# by mean #
smmed_PPH0_avg = apply(smmed_PPH0_mat, 1, mean, na.rm = TRUE)
mmed_PPH0_avg = apply(mmed_PPH0_mat, 1, mean, na.rm = TRUE)
rand_PPH0_avg = apply(rand_PPH0_mat, 1, mean, na.rm = TRUE)
x_at1_PPH0_avg = apply(x_at1_PPH0_mat, 1, mean, na.rm = TRUE)
x_diag_PPH0_avg = apply(x_diag_PPH0_mat, 1, mean, na.rm = TRUE)
x_sf_PPH0_avg = apply(x_sf_PPH0_mat, 1, mean, na.rm = TRUE)

smmed_PPH1_avg = apply(smmed_PPH1_mat, 1, mean, na.rm = TRUE)
mmed_PPH1_avg = apply(mmed_PPH1_mat, 1, mean, na.rm = TRUE)
rand_PPH1_avg = apply(rand_PPH1_mat, 1, mean, na.rm = TRUE)
x_at1_PPH1_avg = apply(x_at1_PPH1_mat, 1, mean, na.rm = TRUE)
x_diag_PPH1_avg = apply(x_diag_PPH1_mat, 1, mean, na.rm = TRUE)
x_sf_PPH1_avg = apply(x_sf_PPH1_mat, 1, mean, na.rm = TRUE)

# by median #
smmed_PPH0_median = apply(smmed_PPH0_mat, 1, median, na.rm = TRUE)
mmed_PPH0_median = apply(mmed_PPH0_mat, 1, median, na.rm = TRUE)
rand_PPH0_median = apply(rand_PPH0_mat, 1, median, na.rm = TRUE)
x_at1_PPH0_median = apply(x_at1_PPH0_mat, 1, median, na.rm = TRUE)
x_diag_PPH0_median = apply(x_diag_PPH0_mat, 1, median, na.rm = TRUE)
x_sf_PPH0_median = apply(x_sf_PPH0_mat, 1, median, na.rm = TRUE)

smmed_PPH1_median = apply(smmed_PPH1_mat, 1, median, na.rm = TRUE)
mmed_PPH1_median = apply(mmed_PPH1_mat, 1, median, na.rm = TRUE)
rand_PPH1_median = apply(rand_PPH1_mat, 1, median, na.rm = TRUE)
x_at1_PPH1_median = apply(x_at1_PPH1_mat, 1, median, na.rm = TRUE)
x_diag_PPH1_median = apply(x_diag_PPH1_mat, 1, median, na.rm = TRUE)
x_sf_PPH1_median = apply(x_sf_PPH1_mat, 1, median, na.rm = TRUE)

ggdata0 = data.table(
  x = 1:(N2 + 1), 
  Random = rand_PPH0_avg, 
  `X=1` = x_at1_PPH0_avg, 
  Diagonal = x_diag_PPH0_avg, 
  SpaceFilling = x_sf_PPH0_avg,
  SeqMED = smmed_PPH0_avg, 
  Hypothesis = rep("H0", idxlast)
)
ggdata1 = data.table(
  x = 1:(N2 + 1), 
  Random = rand_PPH1_avg, 
  `X=1` = x_at1_PPH1_avg,
  Diagonal = x_diag_PPH1_avg,
  SpaceFilling = x_sf_PPH1_avg,
  SeqMED = smmed_PPH1_avg, 
  Hypothesis = rep("H1", idxlast)
)
ggdata = rbind(ggdata0, ggdata1)
ggdata.melted = melt(ggdata, id = c("x", "Hypothesis"), value.name = "epph", 
                     variable.name = "Design")
ggdata.gpvs.H1true = ggdata.melted
plt2 = ggplot(ggdata.melted, aes(x = x, y = epph, color = Design, linetype = Design)) +
  facet_wrap(vars(Hypothesis)) + 
  geom_path() + 
  scale_linetype_manual(values=c(rep("dashed", 4), "solid")) + 
  geom_point(data = ggdata.melted[x == 10], aes(x = x, y = epph)) + 
  theme_bw() + 
  scale_x_continuous(breaks = 1:10) +
  theme(panel.grid.minor = element_blank(), axis.text.x = element_text(size = 5)) + 
  labs(y = "", x = "Stages") + 
  ylim(0, 1)
plt2
# ggsave("gpvsh1_epph.pdf",
#        plot = last_plot(),
#        device = "pdf",
#        path = image_path,
#        scale = 1,
#        width = 4.5,
#        height = 2,
#        units = c("in")
# )

