################################################################################
# last updated: 08/18/21
# purpose: to create a list of seqmed simulations
# dimT = 2:
#   dimensions (1, 2) vs dimensions (1, 2, 3)
#   where the true dimensions are (1, 2)
# dimT = 3:
#   dimensions (1, 2) vs dimensions (1, 2, 3)
#   where the true dimensions are (1, 2, 3)

dimT = 3 # 2, 3

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
Nin = 5
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
sigmasq = 0.3

# hypothesis settings
type01 = c(2, 3)
mu_full = c(0.5, 0.5, 0.5) #
indices0 = c(1, 2) #
indices1 = 1:length(mu_full)
mu0 = rep(0, length(indices0))
mu1 = rep(0, length(indices1))
sigmasq01 = 0.5
V0 = diag(rep(sigmasq01,length(mu0)))
V1 = diag(rep(sigmasq01,length(mu1)))
model0 = list(
  indices = indices0, beta.mean = mu0, beta.var = V0)
model1 = list(
  indices = indices1, beta.mean = mu1, beta.var = V1)

# seqmed settings
p = 3
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

seqmeds = readRDS(paste0(
  output_dir, "/dim", dimT, 
  "_seqmed", 
  "_Nttl", Nttl,
  "_Nin", Nin,
  "_numSeq", numSeq,
  "_seqN", seqN,
  "_seed", rng.seed,
  ".rds"))

################################################################################
# non-sequential designs
################################################################################

# random design
set.seed(12345)
x_random = matrix(runif(n = dimX * Nnew, min = xmin, max = xmax), 
                  nrow = Nnew, ncol = dimX)
y_random = simulateYvs(
  x_random[ , indicesT], Nnew, betaT, sigmasq, 1, seed = 12)

# doptimal design
dopt_pts = as.matrix(expand.grid(c(-1, 1), c(-1, 1), c(-1, 1)))
dopt_pts = dopt_pts[sample(1:dim(dopt_pts)[1], replace = FALSE), ]
x_doptimal = matrix(NA, nrow = Nnew, ncol = dimX)
for(i in 0:(Nnew - 1)){
  x_doptimal[ i + 1, ] = as.matrix(dopt_pts[ 1 + (i %% dim(dopt_pts)[1]), ])
}
y_doptimal = simulateYvs(
  x_doptimal[ , indicesT], Nnew, betaT, sigmasq, 1, seed = 12)

# 3 level factorial design
factorial_pts = as.matrix(expand.grid(c(-1, 0, 1), c(-1, 0, 1), c(-1, 0, 1)))
factorial_pts = factorial_pts[
  sample(1:dim(factorial_pts)[1], replace = FALSE), ]
x_factorial = factorial_pts
y_factorial = simulateYvs(
  x_factorial[ , indicesT], Nnew, betaT, sigmasq, 1, seed = 12)



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

numBreaks = 8
sm = seqmeds[[sim.idx]]

# fix this later
x_random = rbind(sm$x.in, x_random)
y_random = c(sm$y.in, y_random)
x_doptimal = rbind(sm$x.in, x_doptimal)
y_doptimal = c(sm$y.in, y_doptimal)
x_factorial = rbind(sm$x.in, x_factorial)
y_factorial = c(sm$y.in, y_factorial)

maxcounts = rep(NA, length(mu_full))
for(i in 1:length(mu_full)){
  marginal = i
  h = hist(sm$x.new[ , marginal], breaks = numBreaks, plot = FALSE)
  maxcounts[i] = max(h$counts)
}

marginals = matrix(NA, nrow = length((Nin + 1):Nttl), ncol = 3)
for(i in 1:(dim(marginals)[2])) {
  marginals[, i] = sm$x.new[ , i]
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



################################################################################
# plot the posterior probabilities of the hypotheses
################################################################################
source(paste(functions_dir, "/postprob_hypotheses.R", sep = ""))
library(data.table)

models = list("H0" = list(mu0, V0, NULL, indices0),
              "H1" = list(mu1, V1, NULL, indices1))
#EPPHs for other designs
EPPHs_rand = calcEPPH(x_random, Nttl, betaT, NULL, models, sigmasq, numSims, indicesT, seed = 123)
EPPHs_dopt = calcEPPH(x_doptimal, Nttl, betaT, NULL, models, sigmasq, numSims, indicesT, seed = 123)
EPPHs_fact = calcEPPH(x_factorial, Nttl, betaT, NULL, models, sigmasq, numSims, indicesT, seed = 123)
#EPPHs for smmed
EPPHs_smmed = calcEPPHseqdata(
  c(sm$y.in, sm$y.new), rbind(sm$x.in, sm$x.new), models, sigmasq, Nin, numSeq, 
  seqN)
#EPPHs for smmed sims
idxlast = numSeq + 1
EPPH0_smmedsims = matrix(NA, idxlast, numSims)
EPPH1_smmedsims = matrix(NA, idxlast, numSims)
for(i in 1:numSims){
  simH0.numSimstemp = seqmeds[[i]]
  EPPH_temp = calcEPPHseqdata(
    c(simH0.numSimstemp$y.in, simH0.numSimstemp$y.new), 
    rbind(simH0.numSimstemp$x.in, simH0.numSimstemp$x.new), 
    models, sigmasq, Nin, numSeq, seqN)
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
  geom_point(data = ggdata.melted[ggdata.melted$x == 10, ], 
             aes(x = x, y = epph)) + 
  theme_bw() + 
  theme(panel.grid.minor = element_blank()) + 
  labs(y = "", x = "Stages")
plt



################################################################################
# plot the MSE of beta-hat (posterior mean) of the hypotheses
################################################################################
source(paste(functions_dir, "/posterior_mean_mse.R", sep = ""))

if(dimT == 2){
  hyp_mu = mu0
  hyp_V = V0
  hyp_ind = indices0
} else{
  hyp_mu = mu1
  hyp_V = V1
  hyp_ind = indices1
}
mseBn_smmed = getMSEBeta(
  rbind(sm$x.in, sm$x.new), 
  Ntot, betaT, hyp_mu, hyp_V, sigmasq, NULL, hyp_ind)$MSE_postmean
mseBn_rand = getMSEBeta(
  x_random, Ntot, betaT, hyp_mu, hyp_V, sigmasq, NULL, hyp_ind)$MSE_postmean
mseBn_dopt = getMSEBeta(
  x_doptimal, Ntot, betaT, hyp_mu, hyp_V, sigmasq, NULL, hyp_ind)$MSE_postmean
mseBn_fact = getMSEBeta(
  x_factorial, Ntot, betaT, hyp_mu, hyp_V, sigmasq, NULL, hyp_ind)$MSE_postmean

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



################################################################################
# plot the MSE of y-hat (posterior mean) of the hypotheses
################################################################################
source(paste(functions_dir, "/predictive_yhat_mse.R", sep = ""))
