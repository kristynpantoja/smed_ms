################################################################################
# last updated: 12/09/20
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
which.sim = 1

# plot a seqmed
seqmed1 = seqmeds[[which.sim]]
ggdata = data.frame(x = seqmed1$D, y = seqmed1$y)
design.seqmed.includeinit = ggplot(ggdata) + 
  geom_histogram(binwidth = 0.12, closed = "right", aes(x =x, y = after_stat(density))) + 
  theme_bw(base_size = 20) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
data.seqmed.includeinit = ggplot(ggdata) + 
  geom_point(aes(x, y)) +
  stat_function(fun = fT) + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# ggarrange(design.seqmed.includeinit, data.seqmed.includeinit)
# plot without initial data
ggdata = data.frame(x = seqmed1$D[-c(1:seqN)], y = seqmed1$y[-c(1:seqN)])
design.seqmed = ggplot(ggdata) + 
  geom_histogram(binwidth = 0.12, closed = "right", aes(x =x, y = after_stat(density))) + 
  theme_bw(base_size = 20) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
data.seqmed = ggplot(ggdata) + 
  geom_point(aes(x, y)) +
  stat_function(fun = fT) + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# ggarrange(design.seqmed, data.seqmed)

# plot a bh
bh1 = bhs[[which.sim]]
bh1_input = bhs_inputs[[which.sim]]
ggdata = data.frame(
  x = c(bh1$x , bh1$x.new), 
  y = c(bh1$y, bh1$y.new))
design.bh.includeinit = ggplot(ggdata) + 
  geom_histogram(binwidth = 0.12, closed = "right", aes(x =x, y = after_stat(density))) + 
  theme_bw(base_size = 20) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
data.bh.includeinit = ggplot(ggdata) + 
  geom_point(aes(x, y)) +
  stat_function(fun = fT) + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# ggarrange(design.bh.includeinit, data.bh.includeinit)
# plot without initial data
ggdata = data.frame(x = bh1$x.new, y = bh1$y.new)
design.bh = ggplot(ggdata) + 
  geom_histogram(binwidth = 0.12, closed = "right", aes(x =x, y = after_stat(density))) + 
  theme_bw(base_size = 20) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
data.bh = ggplot(ggdata) + 
  geom_point(aes(x, y)) +
  stat_function(fun = fT) + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# ggarrange(design.bh, data.bh)

# side-by-sides of seqmed and bh
ggarrange(
  design.seqmed.includeinit, data.seqmed.includeinit,
  design.bh.includeinit, data.bh.includeinit,
  nrow = 2, ncol = 2)
ggarrange(
  design.seqmed, data.seqmed,
  design.bh, data.bh,
  nrow = 2, ncol = 2)

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

ggdata0 = data.table(
  x = 1:numSeq, 
  Dlinear = rep(epphs_dopt1[1], numSeq), 
  Dquadratic = rep(epphs_dopt2[1], numSeq), 
  SpaceFilling = rep(epphs_space[1], numSeq), 
  SeqMED = epphs_seqmed[ 1, ], 
  BH = epphs_bh[ 1, ], 
  Hypothesis = rep("H0", numSeq)
)
ggdata1 = data.table(
  x = 1:numSeq, 
  Dlinear = rep(epphs_dopt1[2], numSeq), 
  Dquadratic = rep(epphs_dopt2[2], numSeq), 
  SpaceFilling = rep(epphs_space[2], numSeq), 
  SeqMED = epphs_seqmed[ 2, ], 
  BH = epphs_bh[ 2, ], 
  Hypothesis = rep("H1", numSeq)
)
ggdataT = data.table(
  x = 1:numSeq, 
  Dlinear = rep(epphs_dopt1[3], numSeq), 
  Dquadratic = rep(epphs_dopt2[3], numSeq), 
  SpaceFilling = rep(epphs_space[3], numSeq), 
  SeqMED = epphs_seqmed[ 3, ], 
  BH = epphs_bh[ 3, ], 
  Hypothesis = rep("HT", numSeq)
)
ggdata = rbind(ggdata0, ggdata1, ggdataT)
ggdata.melted = melt(ggdata, id = c("x", "Hypothesis"), value.name = "epph", 
                     variable.name = "Design")
epph.plt = ggplot(ggdata.melted, aes(x = x, y = epph, color = Design, linetype = Design)) +
  facet_wrap(vars(Hypothesis)) + 
  geom_path(size = 1) + 
  scale_linetype_manual(values=c(rep("dashed", 3), rep("solid", 2))) + 
  geom_point(data = ggdata.melted[x == numSeq], aes(x = x, y = epph), size = 2) + 
  theme_bw(base_size = 20) + 
  theme(panel.grid.minor = element_blank()) + 
  labs(y = "", x = "Stages") + 
  ylim(0, 1)
epph.plt

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
  Designs = rep(c("Dlinear", "Dquadratic", "SpaceFilling", "SeqMED", "BH"), 4), 
  MSE = c(b0, b1, b2, b3), beta = rep(c("B0", "B1", "B2", "B3"), each = length(b0)))
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
ylimarg = c(0, 0.15)

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
  # scale_y_continuous(limits = ylimarg) + 
  # scale_x_continuous(limits = c(-1, 1)) +
  geom_path() + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  labs(y = "", x = "x")
msey.plt
