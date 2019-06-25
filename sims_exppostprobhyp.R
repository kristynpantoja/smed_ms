#########################################
############# Unknown Slope #############
#########################################


mean_beta0 = 1 # slope of null model
mean_beta1 = 1 / 2 # slope of alternative model
var_mean0 = 0.005; var_mean1 = var_mean0; # variance on beta
var_e = 0.025 # variance on error

xmin = 0 
xmax = 1 

f0 = function(x) mean_beta0 * x # null regression model
f1 = function(x) mean_beta1 * x # alternative regression model

N = 67
type = c(1, 1)
p = 1


library(mvtnorm)
numSims = 1000

model_evidence = function(Y, X, mean_beta, var_mean, var_e){
  # Y is a vector
  # X is a matrix
  # var_mean is a matrix
  # var_e is a scalar
  N = length(Y)
  marginaly_mean = X %*% mean_beta
  marginaly_var = diag(rep(var_e, N)) + (X %*% var_mean %*% t(X))
  return(dmvnorm(Y, mean = marginaly_mean, sigma = marginaly_var, log = FALSE))
}


simulateY = function(X, mean_beta, var_mean, var_e, numSims, plotrandsim = FALSE, dim = 1){
  Y = matrix(rep(NA, N * numSims), N, numSims) # each column is a separate simulation
  for(j in 1:numSims){
    if (dim == 1){
      beta = rnorm(n = 1, mean = mean_beta, sd = sqrt(var_mean))
    }else {
      beta = t(rmvnorm(n = 1, mean = mean_beta, sigma = var_mean))
    }
    for(i in 1:N){
      Y[i, j] = rnorm(n = 1, mean = X[i, ] %*% beta, sd = sqrt(var_e))
    }
  }
  if(plotrandsim == TRUE){
    randSim = sample(1:numSims, 1)
    plot(Y[ , randSim] ~ X)
  }
  return(Y)
}

# testing - it works fine! :D
#D_space = seq(xmin, xmax, length.out = N)
X = as.matrix(D_space)

simY = simulateY(X, mean_beta = mean_beta0, var_mean0, var_e, numSims)
randSim = sample(1:numSims, 1); plot(simY[ , randSim] ~ X, xlim = c(0, 1), ylim = c(0, 1))

simPostH0 = rep(NA, numSims)
simPostH1 = rep(NA, numSims)
simBF01 = rep(NA, numSims)

for(j in 1:numSims){
  Y = simY[, j]
  # get model evidences
  simEvidenceH0 = model_evidence(Y, X, mean_beta0, var_mean0, var_e)
  simEvidenceH1 = model_evidence(Y, X, mean_beta1, var_mean1, var_e)
  # calculate posterior probabilities of models
  simPostH0[j] = simEvidenceH0 / (simEvidenceH0 + simEvidenceH1)
  simPostH1[j] = simEvidenceH1 / (simEvidenceH0 + simEvidenceH1)
  # calculate bayes factor
  simBF01[j] = simPostH0[j] / simPostH1[j]
}
expected_postH0 = mean(simPostH0)
expected_postH1 = mean(simPostH1)
expected_BF01 = mean(simBF01)
expected_BF01

###

simY2 = simulateY2(X, mean_beta = mean_beta0, var_mean0, var_e, numSims, dim = 1)
randSim = sample(1:numSims, 1); plot(simY2[ , randSim] ~ X, xlim = c(0, 1), ylim = c(0, 1))

sim2PostH0 = rep(NA, numSims)
sim2PostH1 = rep(NA, numSims)
sim2BF01 = rep(NA, numSims)

for(j in 1:numSims){
  Y = simY2[, j]
  # get model evidences
  simEvidenceH0 = model_evidence(Y, X, mean_beta0, var_mean0, var_e)
  simEvidenceH1 = model_evidence(Y, X, mean_beta1, var_mean1, var_e)
  # calculate posterior probabilities of models
  sim2PostH0[j] = simEvidenceH0 / (simEvidenceH0 + simEvidenceH1)
  sim2PostH1[j] = simEvidenceH1 / (simEvidenceH0 + simEvidenceH1)
  # calculate bayes factor
  sim2BF01[j] = sim2PostH0[j] / sim2PostH1[j]
}
expected_postH0_2 = mean(sim2PostH0)
expected_postH1_2 = mean(sim2PostH1)
expected_BF01_2 = mean(sim2BF01)
expected_BF01_2




##################################
############# Linear #############
##################################

mean_beta0 = c(0, 1) # slope of null model
mean_beta1 = c(0, 1/2) # slope of alternative model
var_mean0 = diag(c(0.005, 0.005)); var_mean1 = var_mean0 # variance on beta
var_e = 0.025 # variance on error
xmin = 0
xmax = 1 
f0 = function(x) mean_beta0[1] + mean_beta0[2] * x # null regression model
f1 = function(x) mean_beta1[1] + mean_beta1[2] * x # alternative regression model
N = 67
type = c(2, 2)
p = 2

X = cbind(rep(1, N), D_space)

simY = simulateY(X, mean_beta = mean_beta0, var_mean = var_mean0, var_e, numSims, dim = 2)
randSim = sample(1:numSims, 1); plot(simY[ , randSim] ~ D_space, xlim = c(0, 1), ylim = c(0, 1))

simPostH0 = rep(NA, numSims)
simPostH1 = rep(NA, numSims)
simBF01 = rep(NA, numSims)

for(j in 1:numSims){
  Y = simY[, j]
  # get model evidences
  simEvidenceH0 = model_evidence(Y, X, mean_beta0, var_mean0, var_e)
  simEvidenceH1 = model_evidence(Y, X, mean_beta1, var_mean1, var_e)
  # calculate posterior probabilities of models
  simPostH0[j] = simEvidenceH0 / (simEvidenceH0 + simEvidenceH1)
  simPostH1[j] = simEvidenceH1 / (simEvidenceH0 + simEvidenceH1)
  # calculate bayes factor
  simBF01[j] = simPostH0[j] / simPostH1[j]
}
expected_postH0 = mean(simPostH0)
expected_postH1 = mean(simPostH1)
expected_BF01 = mean(simBF01)
expected_BF01

#

X2 = cbind(rep(1, N), D_space)

simY = simulateY(X2, mean_beta = mean_beta1, var_mean = var_mean1, var_e, numSims, dim = 2)
randSim = sample(1:numSims, 1); plot(simY[ , randSim] ~ D_space, xlim = c(0, 1), ylim = c(0, 1))

sim2PostH0 = rep(NA, numSims)
sim2PostH1 = rep(NA, numSims)
sim2BF01 = rep(NA, numSims)

for(j in 1:numSims){
  Y = simY[, j]
  # get model evidences
  simEvidenceH0 = model_evidence(Y, X2, mean_beta0, var_mean0, var_e)
  simEvidenceH1 = model_evidence(Y, X2, mean_beta1, var_mean1, var_e)
  # calculate posterior probabilities of models
  sim2PostH0[j] = simEvidenceH0 / (simEvidenceH0 + simEvidenceH1)
  sim2PostH1[j] = simEvidenceH1 / (simEvidenceH0 + simEvidenceH1)
  # calculate bayes factor
  sim2BF01[j] = sim2PostH0[j] / sim2PostH1[j]
}

expected_postH0_2 = mean(sim2PostH0)
expected_postH1_2 = mean(sim2PostH1)
expected_BF01_2 = mean(sim2BF01)
expected_BF01_2




hyp = c(rep("H0", numSims), rep("H1", numSims))
postH0 = c(simPostH0, sim2PostH0)
postH1 = c(simPostH1, sim2PostH1)
BF01 = c(simBF01, sim2BF01)
df = data.frame("hyp" = hyp, "postH0" = postH0, "postH1" = postH1, "BF01" = BF01)
agg_postH0 = aggregate(df$postH0, by = list(hyp = df$hyp), 
                       FUN = function(x) c(mean = mean(x), sd = sd(x), n = length(x)))
agg_postH0 = do.call(data.frame, agg_postH0)
colnames(agg_postH0) = c("hyp", "mean", "sd", "n")
agg_postH0$se = agg_postH0$sd / sqrt(agg_postH0$n)
agg_postH0$names = c("X1", "X2")
upperlim <- max(agg_postH0$mean) + agg_postH0[agg_postH0$mean == max(agg_postH0$mean), "se"] * 3

barCenters <- barplot(height = agg_postH0$mean,
                      names.arg = agg_postH0$names,
                      beside = true, las = 2,
                      ylim = c(0, upperlim),
                      cex.names = 0.75, xaxt = "n",
                      main = "E[P(Hi | Y) | Hj] for different Hj",
                      ylab = "E[P(Hi | Y) | Hj]",
                      border = "black", axes = TRUE)

# Specify the groupings. We use srt = 45 for a
# 45 degree string rotation
text(x = barCenters, y = par("usr")[3] - 1, srt = 45,
     adj = 1, labels = agg_postH0$names, xpd = TRUE)

segments(barCenters, agg_postH0$mean - agg_postH0$se * 2, barCenters,
         agg_postH0$mean + agg_postH0$se * 2, lwd = 1.5)

arrows(barCenters, agg_postH0$mean - agg_postH0$se * 2, barCenters,
       agg_postH0$mean + agg_postH0$se * 2, lwd = 1.5, angle = 90,
       code = 3, length = 0.05)
