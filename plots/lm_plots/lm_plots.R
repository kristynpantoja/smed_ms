# source files for evaluations

home = "/home/kristyn/Documents/smed_ms"

library(ggplot2)
library(ggpubr)
library(reshape2)
library(data.table)
gg_color_hue = function(n) {
  hues = seq(15, 275, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
image_path = paste0(home, "/plots/lm_plots/gg")

# libraries
library(Matrix)
library(expm)
library(matrixStats)
library(scatterplot3d)
library(knitr)
library(mvtnorm)

# --- Sources for S/MMED --- #
functions_home = paste(home, "/functions", sep="")
source(paste(functions_home, "/wasserstein_distance.R", sep = ""))
source(paste(functions_home, "/charge_function_q.R", sep = ""))
source(paste(functions_home, "/variance_marginal_y.R", sep = ""))
source(paste(functions_home, "/generate_MMED_nodata.R", sep = ""))

source(paste(functions_home, "/posterior_mean.R", sep = ""))
source(paste(functions_home, "/construct_design_matrix.R", sep = ""))
source(paste(functions_home, "/posterior_parameters.R", sep = ""))
source(paste(functions_home, "/add_MMED.R", sep = ""))
source(paste(functions_home, "/SMMED.R", sep = ""))

# --- Sources to evaluate designs : MSE(Bn), E[P(H1|Y,D)] --- #
source(paste(functions_home, "/simulate_y.R", sep = ""))
source(paste(functions_home, "/postprob_hypotheses.R", sep = ""))
source(paste(functions_home, "/posterior_mean_mse.R", sep = ""))
# source(paste(functions_home, "/plot_utils.R", sep = ""))
source(paste(functions_home, "/predictive_yhat_mse.R", sep = ""))

mu0 = c(0, 0)
mu1 = c(0, 0, 0)
typeT = 3
betaT = c(-0.2, -0.4, 0.4)
sigmasq01 = 0.25
sigmasq = 0.1

fT = function(x) betaT[1] + betaT[2] * x + betaT[3] * x^2

# MED design #
V0 = diag(rep(sigmasq01,length(mu0)))
V1 = diag(rep(sigmasq01,length(mu1)))
# function settings (including and based on prior settings above)
f0 = function(x) mu0[1] + mu0[2] * x
f1 = function(x) mu1[1] + mu1[2] * x + mu1[3] * x^2
type01 = c(2, 3)
numCandidates = 10^3
k = 4
xmin = -1
xmax = 1
p = 1
N = 100

numSeq = 10
N_seq = 10
alpha_seq = 1



# --- other designs --- #

space_filling =  seq(from = xmin, to = xmax, length.out = N)
# doptimal for linear model
dopt_linear = c(rep(1,floor(N/2)),rep(-1, N - floor(N/2)))
# doptimal for quadratic model
dopt_quadratic = c(rep(1,floor(N/3)),rep(0,ceiling(N/3)),rep(-1,N - floor(N/3) - ceiling(N/3)))

#############################
# --- case 1: quadratic --- #
#############################

smmed_data = readRDS("run_designs/smmed/designs/smmed1.rds")
smmed_data_list = readRDS("run_designs/smmed/designs/case1smmeds.rds")
numSeqMMED = length(smmed_data_list)

# par(mfrow = c(1,2))
# hist(smmed_data$D, breaks = 20, main = "", xlab = "x")
# plot(x = smmed_data$D, y = smmed_data$y, xlab = "x", ylab = "", col = 1)
# legend("bottomleft", legend = c("data", "true model"), lty = c(NA, 1), pch = c(1, NA))
# curve(fT, add = TRUE)

ggdata = data.frame(x = smmed_data$D, y = smmed_data$y)
plt1 = ggplot(ggdata) + 
  geom_histogram(binwidth = 0.12, closed = "right", aes(x =x, y = after_stat(density))) + 
  theme_bw(base_size = 20) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
plt2 = ggplot(ggdata) + 
  geom_point(aes(x, y)) +
  stat_function(fun = fT) + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggarrange(plt1, plt2)
# ggsave("seqmed_d2.pdf",
#        plot = last_plot(),
#        device = "pdf",
#        path = image_path,
#        scale = 1,
#        width = 8,
#        height = 4,
#        units = c("in")
# )
# ggsave("poster_seqmed_d2_h6.png",
#        plot = last_plot(),
#        device = "png",
#        path = image_path,
#        scale = 1,
#        width = 13.5,
#        height = 6,
#        units = c("in")
# )

ggarrange(plt2, plt1)
# ggsave("poster_seqmedv2_d2_h6.png",
#        plot = last_plot(),
#        device = "png",
#        path = image_path,
#        scale = 1,
#        width = 13.5,
#        height = 6,
#        units = c("in")
# )

# --- Fits --- #

designs = list(dopt_linear, dopt_quadratic, space_filling, smmed_data$D)
design_names = c("Dlinear", "Dquadratic", "SpaceFilling", "SeqMED")
design_col = c(4, 5, 3, 1)

# par(mfrow=c(2,2))

# col_postmeans = c(3, 5, 4)
# plot(x = 1:numSeq, y = smmed_data$postmean1[1, ], type = "l", 
#      ylim = c(-0.5, 0.5), #main = "Posterior Mean, Quadratic", 
#      xlab = "steps 1:10", ylab = "beta", col = col_postmeans[1])
# lines(x = 1:numSeq, y = smmed_data$postmean1[2, ], type = "l", col = col_postmeans[2])
# lines(x = 1:numSeq, y = smmed_data$postmean1[3, ], type = "l", col = col_postmeans[3])
# abline(h = betaT[1], lty = 2, col = col_postmeans[1])
# abline(h = betaT[2], lty = 2, col = col_postmeans[2])
# abline(h = betaT[3], lty = 2, col = col_postmeans[3])
# legend("bottomleft", c("Bn0", "Bn1", "Bn2"), lty = c(1, 1, 1), col = col_postmeans, bg = "white")

# col_postmeans = c(8, 6)
# plot(x = 1:numSeq, y = smmed_data$postmean0[1, ], type = "l", 
#      ylim = c(-0.5, 0.5), #main = "Posterior Mean, Linear", 
#      xlab = "steps 1:10", ylab = "beta", col = col_postmeans[1])
# lines(x = 1:numSeq, y = smmed_data$postmean0[2, ], type = "l", col = col_postmeans[2])
# legend("topleft", c("Bn0", "Bn1"), lty = c(1, 1), col = col_postmeans, bg = "white")

f1est = function(x) smmed_data$postmean1[1, N_seq] + 
  smmed_data$postmean1[2, N_seq] * x + smmed_data$postmean1[3, N_seq] * x^2
# curve(fT, xlim = c(-1, 1))
# curve(f1est, col = 2, add = T)
# legend("topright", c("true model", "estimated model"), lty = c(1, 1), col = c(1, 2))

f2est = function(x) smmed_data$postmean0[1, N_seq] + 
  smmed_data$postmean0[2, N_seq] * x
# curve(fT, xlim = c(-1, 1))
# curve(f2est, col = 2, add = T)

#

numseq = 1e2
x_seq = seq(from = xmin, to = xmax, length.out = numseq)
f1est_seq = sapply(x_seq, f1est)
f2est_seq = sapply(x_seq, f2est)
fT_seq = sapply(x_seq, fT)

ggdata_est = data.table::data.table(
  x = x_seq, 
  `Estimated Quadratic` = f1est_seq, 
  `Estimated Line` = f2est_seq
)
ggdata_est = data.table::melt(ggdata_est, id = c("x"), value.name = "y", variable.name = "Function")

ggdata_true = data.table::data.table(x = x_seq, y = fT_seq)

ggplot(ggdata_est) + 
  facet_wrap(facets = vars(Function)) +
  geom_path(aes(x = x, y = y, color = Function)) + 
  geom_path(data = ggdata_true, aes(x, y, color = "True Quadratic")) + 
  scale_color_manual(values = c(gg_color_hue(3)[c(1, 3)], "black")) +
  # stat_function(fun = fT) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
# ggsave("fits_d2.pdf",
#        plot = last_plot(),
#        device = "pdf",
#        path = image_path,
#        scale = 1,
#        width = 8,
#        height = 4,
#        units = c("in")
# )
# all widths are 13.5
# height = 4
# ggsave("poster_fits_d2_h4.png",
#        plot = last_plot(),
#        device = "png",
#        path = image_path,
#        scale = 1,
#        width = 13.5,
#        height = 4,
#        units = c("in")
# )

#

ggdata = data.frame(x = smmed_data$D, y = smmed_data$y)
plt1.1 = ggplot(ggdata) + 
  geom_histogram(binwidth = 0.12, closed = "right", aes(x =x, y = after_stat(density))) + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

ggdata_est = data.table::data.table(
  x = x_seq, 
  `Estimated Quadratic` = f1est_seq, 
  `Estimated Line` = f2est_seq
)
ggdata_est = data.table::melt(ggdata_est, id = c("x"), value.name = "y", variable.name = "Function")
ggdata_true = data.table::data.table(x = x_seq, y = fT_seq)
plt2.1 = ggplot(ggdata_est) + 
  facet_wrap(facets = vars(Function)) +
  geom_path(aes(x = x, y = y, color = Function)) + 
  geom_path(data = ggdata_true, aes(x, y, color = "True Quadratic")) + 
  scale_color_manual(values = c(gg_color_hue(3)[c(1, 3)], "black")) +
  # stat_function(fun = fT) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.position = "none")
ggarrange(plt1.1, plt2.1, widths = c(1, 2))
# height = 4
# ggsave("poster_seqmed_fits_d2_h4.png",
#        plot = last_plot(),
#        device = "png",
#        path = image_path,
#        scale = 1,
#        width = 13.5,
#        height = 4,
#        units = c("in")
# )


#


ggdata = data.frame(x = smmed_data$D, y = smmed_data$y)
yrange = range(smmed_data$y)
plt1 = ggplot(ggdata) + 
  geom_histogram(binwidth = 0.12, closed = "right", aes(x =x, y = after_stat(density))) + 
  theme_bw(base_size = 20) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
plt2 = ggplot(ggdata) + 
  geom_point(aes(x, y), size = 2) +
  stat_function(fun = fT, size = 2) + 
  scale_y_continuous(limits = yrange) + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggdata_est = data.table::data.table(
  x = x_seq, 
  `Estimated Quadratic` = f1est_seq, 
  `Estimated Line` = f2est_seq
)
ggdata_est = data.table::melt(ggdata_est, id = c("x"), value.name = "y", variable.name = "Function")

ggdata_true = data.table::data.table(x = x_seq, y = fT_seq)

plt3 = ggplot(ggdata_est) + 
  facet_wrap(facets = vars(Function)) +
  geom_path(aes(x = x, y = y, color = Function)) + 
  scale_y_continuous(limits = yrange) + 
  geom_path(data = ggdata_true, aes(x, y, color = "True Quadratic")) + 
  scale_color_manual(values = c(gg_color_hue(3)[c(1, 3)], "black")) +
  # stat_function(fun = fT) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
ggarrange(plt1, plt2, plt3, ncol = 3, widths = c(1, 1, 2))
# height = 4
# ggsave("poster_seqmed_data_fits_d2_h4.png",
#        plot = last_plot(),
#        device = "png",
#        path = image_path,
#        scale = 1,
#        width = 13.5,
#        height = 4,
#        units = c("in")
# )


#


ggdata = data.frame(x = smmed_data$D, y = smmed_data$y)
yrange = range(smmed_data$y)
plt1 = ggplot(ggdata) + 
  geom_histogram(binwidth = 0.12, closed = "right", aes(x =x, y = after_stat(density))) + 
  theme_bw(base_size = 20) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
plt2 = ggplot(ggdata) + 
  geom_point(aes(x, y), size = 2) +
  scale_y_continuous(limits = yrange) + 
  theme_bw(base_size = 20) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggdata_est = data.table::data.table(
  x = x_seq, 
  `Est Quadratic` = f1est_seq, 
  `Est Line` = f2est_seq, 
  `True Quadratic` = fT_seq
)
ggdata_est = data.table::melt(ggdata_est, id = c("x"), value.name = "y", variable.name = "Function")

ggdata_true = data.table::data.table(x = x_seq, y = fT_seq)

plt3 = plt2 + 
  geom_path(data = ggdata_est, mapping = aes(x = x, y = y, 
                                             color = Function), 
            inherit.aes = FALSE, size = 2) + 
  scale_color_manual(values = c(gg_color_hue(3)[c(1, 3)], 1)) +
  # stat_function(fun = fT) +
  theme_bw(base_size = 20) + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
ggarrange(plt1, plt3, ncol = 2, widths = c(1, 1.6))
# height= 5, 6
ggsave("poster_seqmed_altogether_d2_h6.png",
       plot = last_plot(),
       device = "png",
       path = image_path,
       scale = 1,
       width = 13.5,
       height = 6,
       units = c("in")
)
ggarrange(plt3, plt1, ncol = 2, widths = c(1.5, 1))
# ggsave("poster_seqmed_altogether2_d2_h6.png",
#        plot = last_plot(),
#        device = "png",
#        path = image_path,
#        scale = 1,
#        width = 13.5,
#        height = 6,
#        units = c("in")
# )


# --- wasserstein distance --- #

# par(mfrow = c(1, 1))

# testx = seq(from = xmin, to = xmax, length.out = numCandidates)
# is_already_in_D = testx %in% smmed_data$D
# testx = testx[!is_already_in_D]
# wass_testx = sapply(testx, function(x) Wasserstein_distance_postpred(x, smmed_data$postmean0[,10], smmed_data$postmean1[,10], 
#                                                                      diag(smmed_data$postvar0[,10]), diag(smmed_data$postvar1[,10]), sigmasq, type01))
# plot(x = testx, y = (wass_testx), type = "l", ylim = range(-0.5, 0.5), 
#      ylab = "wasserstein(x)", xlab = "", main = "")
# f0data = function(x) smmed_data$postmean0[,10][1] + smmed_data$postmean0[,10][2] * x
# f1data = function(x) smmed_data$postmean1[,10][1] + smmed_data$postmean1[,10][2] * x + smmed_data$postmean1[,10][3] * x^2
# fT = function(x) betaT[1] + betaT[2] * x + betaT[3] * x^2
# curve(f0data, col = 2, lwd = 5, add = T)
# curve(f1data, add = T, col = 5, lty = 2, lwd = 5)
# curve(fT, add = T, col = 1, lty = 3, lwd = 5)
# 
# # left shade
# xshade = seq(from = -1, to = -0.6, length.out = 100)
# yshade1 = sapply(xshade, FUN = f0data)
# yshade2 = sapply(xshade, FUN = f1data)
# polygon(c(xshade,rev(xshade)),c(yshade2,rev(yshade1)),col=rgb(1, 0, 0, 0.1), border = NA)
# 
# # middle shade
# xshade = seq(from = -0.6, to = 0.6, length.out = 100)
# yshade1 = sapply(xshade, FUN = f0data)
# yshade2 = sapply(xshade, FUN = f1data)
# polygon(c(xshade,rev(xshade)),c(yshade2,rev(yshade1)),col=rgb(1, 0, 0, 0.1), border = NA)
# 
# # right shade
# xshade = seq(from = 0.6, to = 1, length.out = 100)
# yshade1 = sapply(xshade, FUN = f0data)
# yshade2 = sapply(xshade, FUN = f1data)
# polygon(c(xshade,rev(xshade)),c(yshade2,rev(yshade1)),col=rgb(1, 0, 0, 0.1), border = NA)
# 
# legend("bottomleft", c("f0", "f1", "true f"), lty = c(1,2,3), lwd = 5, col = c(2, 5, 1))



w_seq = sapply(x_seq, function(x) Wasserstein_distance_postpred(x, smmed_data$postmean0[,10], smmed_data$postmean1[,10], 
                                                                     diag(smmed_data$postvar0[,10]), diag(smmed_data$postvar1[,10]), sigmasq, type01))

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
#        width = 8,
#        height = 4,
#        units = c("in")
# )


# --- EPPH --- #

exppostprobs_space = calcExpPostProbH(space_filling, N, betaT, mu0, mu1, V0, V1, sigmasq, numSims = 100, typeT, type01, seed = 123)
exppostprobs_dopt1 = calcExpPostProbH(dopt_linear, N, betaT, mu0, mu1, V0, V1, sigmasq, numSims = 100, typeT, type01, seed = 123)
exppostprobs_dopt2 = calcExpPostProbH(dopt_quadratic, N, betaT, mu0, mu1, V0, V1, sigmasq, numSims = 100, typeT, type01, seed = 123)

changing_postprobs = list()
# postprobs0 = rep(NA, numSeq)
# postprobs1 = rep(NA, numSeq)
# BF01s = rep(NA, numSeq)
# for(i in 1:numSeq){
#   changing_postprobs[[i]] = calcExpPostProbH_data(smmed_data$y[1:(N_seq * i)], 
#                                                   smmed_data$D[1:(N_seq * i)], 
#                                                   N = N_seq * i, mu0, V0, mu1, V1, sigmasq, type01)
#   postprobs0[i] = changing_postprobs[[i]][1]
#   postprobs1[i] = changing_postprobs[[i]][2]
#   BF01s[i] = changing_postprobs[[i]][3]
# }

# calculate posterior probabilities of hypotheses based on this data
postprobs0seq = matrix(NA, numSeq, numSeqMMED)
postprobs1seq = matrix(NA, numSeq, numSeqMMED)
BF01seq = matrix(NA, numSeq, numSeqMMED)
for(k in 1:numSeqMMED){
  smmed_data_k = smmed_data_list[[k]]
  for(i in 1:numSeq){
    changing_postprobs = calcExpPostProbH_data(smmed_data_k$y[1:(N_seq * i)], 
                                               smmed_data_k$D[1:(N_seq * i)], 
                                               N = N_seq * i, mu0, V0, mu1, V1, sigmasq, type01)
    postprobs0seq[i, k] = changing_postprobs[1]
    postprobs1seq[i, k] = changing_postprobs[2]
    BF01seq[i, k] = changing_postprobs[3]
  }
}

postprobs0seq_mean = apply(postprobs0seq, 1, mean)
postprobs1seq_mean = apply(postprobs1seq, 1, mean)
BF01seq_mean = apply(BF01seq, 1, mean)

# par(mfrow=c(1,2))
# 
# plot(x = 1:numSeq, y = postprobs0seq_mean, type = "l", main = "", 
#      xlab = "Number of Stages", ylab = "", ylim = c(0, 1))
# # for(k in 1:numSeqMMED){
# #   lines(x = 1:numSeq, y = postprobs0seq[,k], col=rgb(0, 0, 0, 0.1))
# # }
# # lines(x = 1:numSeq, y = postprobs0seq_mean, lty = 2)
# abline(h = exppostprobs_space[1], col = 3)
# abline(h = exppostprobs_dopt1[1], col = 4)
# abline(h = exppostprobs_dopt2[1], col = 5)
# legend("topleft", c(design_names), lty = c(rep(1, length(design_col))), col = c(design_col), bg = "white")
# 
# plot(x = 1:numSeq, y = postprobs1seq_mean, type = "l", main = "", 
#      xlab = "Number of Sequential Steps", ylab = "", ylim = c(0, 1))
# # for(k in 1:numSeqMMED){
# #   lines(x = 1:numSeq, y = postprobs1seq[,k], col=rgb(0, 0, 0, 0.1))
# # }
# # lines(x = 1:numSeq, y = postprobs1seq_mean, lty = 2)
# abline(h = exppostprobs_space[2], col = 3)
# abline(h = exppostprobs_dopt1[2], col = 4)
# abline(h = exppostprobs_dopt2[2], col = 5)

#

ggdata0 = data.table(
  x = 1:numSeq, 
  Dlinear = rep(exppostprobs_dopt1[1], numSeq), 
  Dquadratic = rep(exppostprobs_dopt2[1], numSeq), 
  SpaceFilling = rep(exppostprobs_space[1], numSeq), 
  SeqMED = postprobs0seq_mean, 
  Hypothesis = rep("H0", numSeq)
)
ggdata1 = data.table(
  x = 1:numSeq, 
  Dlinear = rep(exppostprobs_dopt1[2], numSeq), 
  Dquadratic = rep(exppostprobs_dopt2[2], numSeq), 
  SpaceFilling = rep(exppostprobs_space[2], numSeq), 
  SeqMED = postprobs1seq_mean, 
  Hypothesis = rep("H1", numSeq)
)
ggdata = rbind(ggdata0, ggdata1)
ggdata.melted = melt(ggdata, id = c("x", "Hypothesis"), value.name = "epph", 
                     variable.name = "Design")
plt6 = ggplot(ggdata.melted, aes(x = x, y = epph, color = Design, linetype = Design)) +
  facet_wrap(vars(Hypothesis)) + 
  geom_path(size = 2) + 
  scale_linetype_manual(values=c(rep("dashed", 3), "solid")) + 
  geom_point(data = ggdata.melted[x == 10], aes(x = x, y = epph), size = 2) + 
  theme_bw(base_size = 20) + 
  theme(panel.grid.minor = element_blank()) + 
  labs(y = "", x = "Stages")
plt6
# ggsave("epph_d2.pdf",
#        plot = last_plot(),
#        device = "pdf",
#        path = image_path,
#        scale = 1,
#        width = 8,
#        height = 4,
#        units = c("in")
# )



# --- MSE(Bn) --- #

postvar_doptlin = getClosedMSE(dopt_linear, N, betaT, mu1, V1, sigmasq, typeT)$MSE_postmean
postvar_doptquad = getClosedMSE(dopt_quadratic, N, betaT, mu1, V1, sigmasq, typeT)$MSE_postmean
postvar_space = getClosedMSE(space_filling, N, betaT, mu1, V1, sigmasq, typeT)$MSE_postmean
postvar_smmed = getClosedMSE(smmed_data$D, N, betaT, mu1, V1, sigmasq, typeT)$MSE_postmean

# par(mfrow = c(1,3))
# barplot(c(postvar_doptlin[1], postvar_doptquad[1], postvar_space[1], postvar_smmed[1]), 
#         names.arg = design_names, main = "B0")
# barplot(c(postvar_doptlin[2], postvar_doptquad[2], postvar_space[2], postvar_smmed[2]), 
#         names.arg = design_names, main = "B1")
# barplot(c(postvar_doptlin[3], postvar_doptquad[3], postvar_space[3], postvar_smmed[3]), 
#         names.arg = design_names, main = "B2")

#

b0 = c(postvar_doptlin[1], postvar_doptquad[1], postvar_space[1], postvar_smmed[1])
b1 = c(postvar_doptlin[2], postvar_doptquad[2], postvar_space[2], postvar_smmed[2])
b2 = c(postvar_doptlin[3], postvar_doptquad[3], postvar_space[3], postvar_smmed[3])

ggdata = data.frame(Designs = rep(c("Dlinear", "Dquadratic", "SpaceFilling", "SeqMED"), 3), 
                    MSE = c(b0, b1, b2), beta = rep(c("B0", "B1", "B2"), each = length(b0)))
ggplot(ggdata, aes(x = Designs, y = MSE)) + 
  geom_bar(stat = "identity") +
  facet_wrap(vars(beta)) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.text.x = element_text(angle = 45, vjust = 0.5)) +
  labs(y = NULL)
# ggsave("msebeta_d2.pdf",
#        plot = last_plot(),
#        device = "pdf",
#        path = image_path,
#        scale = 1,
#        width = 8,
#        height = 4,
#        units = c("in")
# )

ggdataMSEBn = data.table(
  Designs = c("Dlinear", "Dquadratic", "SpaceFilling", "SeqMED"), 
  MSE = b2
)
ggdataEPPH1 = data.table(
  x = 1:numSeq, 
  Dlinear = rep(exppostprobs_dopt1[2], numSeq), 
  Dquadratic = rep(exppostprobs_dopt2[2], numSeq), 
  SpaceFilling = rep(exppostprobs_space[2], numSeq), 
  SeqMED = postprobs1seq_mean
)
ggdataEPPH1 = melt(ggdataEPPH1, id = c("x"), value.name = "epph", 
                   variable.name = "Design")
plt1 = ggplot(ggdataEPPH1, aes(x = x, y = epph, color = Design, linetype = Design)) +
  geom_path(size = 2) + 
  scale_linetype_manual(values=c(rep("dashed", 3), "solid")) + 
  geom_point(data = ggdataEPPH1[x == 10], aes(x = x, y = epph), size = 2) + 
  theme_bw(base_size = 20) + 
  theme(panel.grid.minor = element_blank()) + 
  labs(y = "", x = "Stages")
plt2 = ggplot(ggdataMSEBn, aes(x = Designs, y = MSE)) + 
  geom_bar(stat = "identity") +
  theme_bw(base_size = 20) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.text.x = element_text(angle = 45, vjust = 0.5)) +
  labs(y = NULL)
ggarrange(plt2, plt1, widths = c(1, 1.5))
# height = 3, 4
ggsave("poster_msebn_epph1_h4.png",
       plot = last_plot(),
       device = "png",
       path = image_path,
       scale = 1,
       width = 13.5,
       height = 4,
       units = c("in")
)

ggarrange(plt2, plt6, widths = c(1, 2))
# height = 3, 4
ggsave("poster_msebn_epph_h4.png",
       plot = last_plot(),
       device = "png",
       path = image_path,
       scale = 1,
       width = 13.5,
       height = 4,
       units = c("in")
)




# --- MSE(y-hat) --- #
x_seq2 = seq(from = -1.25, to = 1.25, length.out = 1e4)
yhatmse_space = getClosedMSEyhat_seq(x_seq2, space_filling, N, betaT, typeT, mu1, V1, sigmasq, type01[2])
yhatmse_doptquad = getClosedMSEyhat_seq(x_seq2, dopt_quadratic, N, betaT, typeT, mu1, V1, sigmasq, type01[2])
yhatmse_smmed = getClosedMSEyhat_seq(x_seq2, smmed_data$D, N, betaT, typeT, mu1, V1, sigmasq, type01[2])
yhatmse_doptlin = getClosedMSEyhat_seq(x_seq2, dopt_linear, N, betaT, typeT, mu1, V1, sigmasq, type01[2])

# par(mfrow = c(1,1))
ylimarg = range(0, yhatmse_space$MSEyhat, yhatmse_doptquad$MSEyhat, yhatmse_smmed$MSEyhat)
# plot(x_seq2, yhatmse_space$MSEyhat, type = "l", col = 3, ylim = ylimarg, 
#      xlab = "x in [-1, 1]", ylab = "MSE(y-hat)")
# lines(x_seq2, yhatmse_doptquad$MSEyhat, col = 5)
# lines(x_seq2, yhatmse_smmed$MSEyhat, col = 1)
# lines(x_seq2, yhatmse_doptlin$MSEyhat, col = 4)
# legend("top", design_names, col = design_col, lty = 1)



ggdata = data.table(
  x = x_seq2, 
  Dlinear = yhatmse_doptlin$MSEyhat, 
  Dquadratic = yhatmse_doptquad$MSEyhat, 
  SpaceFilling = yhatmse_space$MSEyhat, 
  SeqMED = yhatmse_smmed$MSEyhat
)
ggdata = melt(ggdata, id = c("x"), value.name = "yhatmse", variable.name = "Design")
ggplot(ggdata, aes(x = x, y = yhatmse, color = Design)) +
  coord_cartesian(ylim = ylimarg, xlim = c(-1, 1)) +
  # scale_y_continuous(limits = ylimarg) + 
  # scale_x_continuous(limits = c(-1, 1)) +
  geom_path() + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  labs(y = "", x = "x")
# ggsave("mseyhat_d2.pdf",
#        plot = last_plot(),
#        device = "pdf",
#        path = image_path,
#        scale = 1,
#        width = 8,
#        height = 4,
#        units = c("in")
# )



#############################
# --- case 2: cubic ------- #
#############################

# MED design #
sigmasq01 = 0.25
sigmasq = 0.1

mu0 = c(0, 0)
mu1 = c(0, 0, 0)
betaT = c(0, -0.75, 0, 1)
typeT = 4
f0 = function(x) mu0[1] + mu0[2] * x
f1 = function(x) mu1[1] + mu1[2] * x + mu1[3] * x^2
xmin = -1
xmax = 1
V0 = diag(rep(sigmasq01,length(mu0)))
V1 = diag(rep(sigmasq01,length(mu1)))
type01 = c(2, 3)
numCandidates = 10^3
k = 4
p = 1
N = 100

fT = function(x) betaT[1] + betaT[2] * x + betaT[3] * x^2 + betaT[4] * x^3

numSeq = 10 
N_seq = 10
alpha_seq = 1
smmed_data2 = readRDS(paste(home, "/run_designs/smmed/designs/smmed2.rds", sep = ""))
smmed_data2_list = readRDS(paste(home, "/run_designs/smmed/designs/case2smmeds.rds", sep = ""))
numSeqMMED = length(smmed_data2_list)

# par(mfrow = c(1,2))
# hist(smmed_data2$D, breaks = 20, main = "", xlab = "x")
# plot(x = smmed_data2$D, y = smmed_data2$y, xlab = "x", ylab = "")
# legend("bottomleft", legend = c("data", "true model"), lty = c(NA, 1), pch = c(1, NA))
# curve(fT, add = TRUE)

ggdata = data.frame(x = smmed_data2$D, y = smmed_data2$y)
yrange = range(smmed_data2$y)
plt1 = ggplot(ggdata) + 
  geom_histogram(binwidth = 0.12, closed = "right", aes(x =x, y = after_stat(density))) + 
  theme_bw(base_size = 20) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
plt2 = ggplot(ggdata) + 
  geom_point(aes(x, y), size = 2) +
  stat_function(fun = fT, size = 2) + 
  theme_bw(base_size = 20) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggarrange(plt1, plt2)
# ggsave("seqmed_d3.pdf",
#        plot = last_plot(),
#        device = "pdf",
#        path = image_path,
#        scale = 1,
#        width = 8,
#        height = 4,
#        units = c("in")
# )

# --- fits --- #

# par(mfrow=c(2,3))
# 
# col_postmeans = c(1, 2, 3, 4)
smmed_data2.postmean2 = matrix(NA, length(betaT), numSeq)
for(k in 1:numSeq){
  smmed_data2.postmean2[ , k] = as.vector(postmean(smmed_data2$y[1:(N_seq * k)],
                                                   smmed_data2$D[1:(N_seq * k)], (N_seq * k),
                                                   c(0, 0, 0, 0), diag(rep(sigmasq01, 4)), sigmasq, 4))
}
# plot(x = 1:numSeq, y = smmed_data2.postmean2[1, ], type = "l", ylim = range(betaT, smmed_data2.postmean2), xlab = "steps 1:10", ylab = "beta", col = col_postmeans[1])
# lines(x = 1:numSeq, y = smmed_data2.postmean2[2, ], type = "l", col = col_postmeans[2])
# lines(x = 1:numSeq, y = smmed_data2.postmean2[3, ], type = "l", col = col_postmeans[3])
# lines(x = 1:numSeq, y = smmed_data2.postmean2[4, ], type = "l", col = col_postmeans[4])
# abline(h = betaT[1], lty = 1, col = col_postmeans[1]) # even though i said lty = 2, it gets covered.
# abline(h = betaT[2], lty = 2, col = col_postmeans[2])
# abline(h = betaT[3], lty = 2, col = col_postmeans[3])
# abline(h = betaT[4], lty = 2, col = col_postmeans[4])
# legend("bottomleft", c("Bn0", "Bn1", "Bn2", "Bn3"), lty = 1, 
#        col = col_postmeans, bg = "white")
# 
# col_postmeans = c(1, 2, 3)
# plot(x = 1:numSeq, y = smmed_data2$postmean1[1, ], type = "l", ylim = range(betaT, smmed_data2.postmean2), xlab = "steps 1:10", ylab = "beta", col = col_postmeans[1])
# lines(x = 1:numSeq, y = smmed_data2$postmean1[2, ], type = "l", col = col_postmeans[2])
# lines(x = 1:numSeq, y = smmed_data2$postmean1[3, ], type = "l", col = col_postmeans[3])
# legend("bottomleft", c("Bn0", "Bn1", "Bn2"), lty = c(1, 1, 1), col = col_postmeans, bg = "white")
# 
# col_postmeans = c(1, 2)
# plot(x = 1:numSeq, y = smmed_data2$postmean0[1, ], type = "l", ylim = range(betaT, smmed_data2.postmean2), xlab = "steps 1:10", ylab = "beta", col = col_postmeans[1])
# # would plot error bars from posterior variances, but they're soooo small... why?
# lines(x = 1:numSeq, y = smmed_data2$postmean0[2, ], type = "l", col = col_postmeans[2])
# legend("bottomleft", c("Bn0", "Bn1"), lty = c(1, 1), col = col_postmeans, bg = "white")

smmed_data2.postmean2.1 = postmean(smmed_data2$y, smmed_data2$D, N,
                                       c(0, 0, 0, 0), diag(rep(sigmasq01, 4)), sigmasq, 4)
# fEst = function(x) smmed_data2.postmean2.1[1] + smmed_data2.postmean2.1[2] * x + smmed_data2.postmean2.1[3] * x^2 + smmed_data2.postmean2.1[4] * x^3
fTest = function(x) smmed_data2.postmean2[1, N_seq] + 
  smmed_data2.postmean2[2, N_seq] * x + 
  smmed_data2.postmean2[3, N_seq] * x^2 + 
  smmed_data2.postmean2[4, N_seq] * x^3
# curve(fT, xlim = c(-1, 1))
# curve(f1est, col = 2, add = T)
# legend("bottomleft", c("true model", "estimated model"), lty = c(1, 1), col = c(1, 2))

f2est = function(x) smmed_data2$postmean1[1, N_seq] + 
  smmed_data2$postmean1[2, N_seq] * x + 
  smmed_data2$postmean1[3, N_seq] * x^2
# curve(fT, xlim = c(-1, 1))
# curve(f2est, col = 2, add = T)

f1est = function(x) smmed_data2$postmean0[1, N_seq] + 
  smmed_data2$postmean0[2, N_seq] * x
# curve(fT, xlim = c(-1, 1))
# curve(fTest, col = 2, add = T)

#

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

ggplot(ggdata_est) + 
  facet_wrap(facets = vars(Function)) +
  geom_path(aes(x = x, y = y, color = Function), size = 2) + 
  geom_path(data = ggdata_true, aes(x, y, color = "True Cubic"), size = 2) + 
  scale_color_manual(values = c(gg_color_hue(3), "black")) +
  scale_y_continuous(limits = yrange) + 
  # stat_function(fun = fT) +
  theme_bw(base_size = 20) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# ggsave("fits_d3.pdf",
#        plot = last_plot(),
#        device = "pdf",
#        path = image_path,
#        scale = 1,
#        width = 8,
#        height = 4,
#        units = c("in")
# )
# height = 3, 4
# ggsave("poster_fits_d3_h4.png",
#        plot = last_plot(),
#        device = "png",
#        path = image_path,
#        scale = 1,
#        width = 13.5,
#        height = 4,
#        units = c("in")
# )

#

plt1 = ggplot(ggdata) + 
  geom_histogram(binwidth = 0.12, closed = "right", aes(x =x, y = after_stat(density))) + 
  theme_bw(base_size = 20) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
plt2 = ggplot(ggdata) + 
  geom_point(aes(x, y)) +
  stat_function(fun = fT, size = 2) + 
  theme_bw(base_size = 20) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
plt3 = ggplot(ggdata_est) + 
  facet_wrap(facets = vars(Function)) +
  geom_path(aes(x = x, y = y, color = Function), size = 2) + 
  geom_path(data = ggdata_true, aes(x, y, color = "True Quadratic"), size = 2) + 
  scale_color_manual(values = c(gg_color_hue(3), "black")) +
  scale_y_continuous(limits = yrange) + 
  # stat_function(fun = fT) +
  theme_bw(base_size = 20) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
###
ggarrange(plt2, plt3, nrow = 1, widths = c(1, 3.5))
# hegiht = 3, 4
# ggsave("poster_data_fits_d3_h3.png",
#        plot = last_plot(),
#        device = "png",
#        path = image_path,
#        scale = 1,
#        width = 13.5,
#        height = 3,
#        units = c("in")
# )
###
plt3.1 = ggplot(ggdata_est) + 
  facet_wrap(facets = vars(Function)) +
  geom_path(aes(x = x, y = y, color = Function), size = 2) + 
  geom_path(data = ggdata_true, aes(x, y, color = "True Quadratic"), size = 2) + 
  scale_color_manual(values = c(gg_color_hue(3), "black")) +
  # stat_function(fun = fT) +
  theme_bw(base_size = 20) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggarrange(plt1, plt3.1, nrow = 1, widths = c(1,  3))
# 12, 14 in
ggsave("poster_seqmed_fits_d3_h4.png",
       plot = last_plot(),
       device = "png",
       path = image_path,
       scale = 1,
       width = 13.5,
       height = 4,
       units = c("in")
)
###

ggdata_est_plt4 = data.table::data.table(
  x = x_seq, 
  `Estimated Line` = f1est_seq, 
  `Estimated Quadratic` = f2est_seq, 
  `Estimated Cubic` = fTest_seq, 
  `True Cubic` = fT_seq
)
ggdata_est_plt4 = data.table::melt(ggdata_est_plt4, id = c("x"), value.name = "y", variable.name = "Function")

plt4 = ggplot(ggdata) + 
  geom_point(aes(x, y), size = 2) +
  geom_path(data = ggdata_est_plt4, 
            mapping = aes(x = x, y = y, color = Function), 
            inherit.aes = FALSE, size = 2) + 
  scale_color_manual(values = c(gg_color_hue(3), 1)) +
  theme_bw(base_size = 20) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggarrange(plt1, plt4, nrow = 1, ncol = 2, widths = c(1, 2))
# heght = 3, 4
# ggsave("poster_datafitsaltogether_d3_h2.png",
#        plot = last_plot(),
#        device = "png",
#        path = image_path,
#        scale = 1,
#        width = 13.5,
#        height = 2,
#        units = c("in")
# )
###

ggarrange(plt1, plt2, plt3, nrow = 1, widths = c(1, 1.1, 4))
# heght = 3, 4
# ggsave("poster_seqmeddatafits_d3_h4.png",
#        plot = last_plot(),
#        device = "png",
#        path = image_path,
#        scale = 1,
#        width = 13.5,
#        height = 4,
#        units = c("in")
# )
# --- high density areas (wasserstein distance) --- #

# par(mfrow = c(1, 1))
# testx = seq(from = xmin, to = xmax, length.out = numCandidates)
# is_already_in_D = testx %in% smmed_data$D
# testx = testx[!is_already_in_D]
# wass_testx = sapply(testx, function(x) Wasserstein_distance_postpred(x, smmed_data$postmean0[,10], smmed_data$postmean1[,10], 
#                                                                      diag(smmed_data$postvar0[,10]), diag(smmed_data$postvar1[,10]), sigmasq, type01))
# plot(x = testx, y = (wass_testx), type = "l", ylim = range(-1, 1), 
#      main = "wasserstein(x)")
# f0data = function(x) smmed_data$postmean0[,10][1] + smmed_data$postmean0[,10][2] * x
# f1data = function(x) smmed_data$postmean1[,10][1] + smmed_data$postmean1[,10][2] * x + smmed_data$postmean1[,10][3] * x^2
# fT = function(x) betaT[1] + betaT[2] * x + betaT[3] * x^2
# curve(f0data, col = 2, lwd = 5, add = T)
# curve(f1data, add = T, col = 5, lty = 2, lwd = 5)
# curve(fT, add = T, col = 1, lty = 3, lwd = 5)
# 
# # left shade
# xshade = seq(from = -1, to = -0.6, length.out = 100)
# yshade1 = sapply(xshade, FUN = f0data)
# yshade2 = sapply(xshade, FUN = f1data)
# polygon(c(xshade,rev(xshade)),c(yshade2,rev(yshade1)),col=rgb(1, 0, 0, 0.1), border = NA)
# 
# # middle shade
# xshade = seq(from = -0.6, to = 0.6, length.out = 100)
# yshade1 = sapply(xshade, FUN = f0data)
# yshade2 = sapply(xshade, FUN = f1data)
# polygon(c(xshade,rev(xshade)),c(yshade2,rev(yshade1)),col=rgb(1, 0, 0, 0.1), border = NA)
# 
# # right shade
# xshade = seq(from = 0.6, to = 1, length.out = 100)
# yshade1 = sapply(xshade, FUN = f0data)
# yshade2 = sapply(xshade, FUN = f1data)
# polygon(c(xshade,rev(xshade)),c(yshade2,rev(yshade1)),col=rgb(1, 0, 0, 0.1), border = NA)
# 
# legend("bottomright", c("f0", "f1", "true f"), lty = c(1,2,3), lwd = 5, col = c(2, 5, 1))



# --- EPPH --- #

designs = list(dopt_linear, dopt_quadratic, space_filling, smmed_data2$D)
design_names = c("Dlinear", "Dquadratic", "SpaceFilling", "SeqMED")
design_col = c(4, 5, 3, 1)

# what are expected posterior probabilities of hypotheses for space-filling design? d-optimal designs? 
models = list("H0" = list(mu0, V0, 2),
              "H1" = list(mu1, V1, 3),
              "H2" = list(betaT, diag(rep(sigmasq01, 4)), 4))
numSims = 100
exppostprobs_space = calcEPPH(space_filling, N, betaT, typeT, models, sigmasq, numSims, seed = 123)
exppostprobs_dopt1 = calcEPPH(dopt_linear,  N, betaT, typeT, models, sigmasq, numSims, seed = 123)
exppostprobs_dopt2 = calcEPPH(dopt_quadratic, N, betaT, typeT, models, sigmasq, numSims, seed = 123)

# calculate posterior probabilities of hypotheses based on this data
changing_postprobs = list() 
postprobs = matrix(NA, length(models), numSeq)
for(i in 1:numSeq){
  changing_postprobs[[i]] = calcEPPHdata(smmed_data2$y[1:(N_seq * i)],
                                                          smmed_data2$D[1:(N_seq * i)], 
                                                          N = N_seq * i, models, sigmasq)
  postprobs[ , i] = changing_postprobs[[i]]
}

postprobs0seq = matrix(NA, numSeq, numSeqMMED)
postprobs1seq = matrix(NA, numSeq, numSeqMMED)
postprobs2seq = matrix(NA, numSeq, numSeqMMED)
for(k in 1:numSeqMMED){
  smmed_data_k = smmed_data2_list[[k]]
  for(i in 1:numSeq){
    changing_postprobs = calcEPPHdata(smmed_data_k$y[1:(N_seq * i)], 
                                                       smmed_data_k$D[1:(N_seq * i)], 
                                                       N = N_seq * i, models, sigmasq)
    postprobs0seq[i, k] = changing_postprobs[1]
    postprobs1seq[i, k] = changing_postprobs[2]
    postprobs2seq[i, k] = changing_postprobs[3]
  }
}

postprobs0seq_mean = apply(postprobs0seq, 1, mean)
postprobs1seq_mean = apply(postprobs1seq, 1, mean)
postprobs2seq_mean = apply(postprobs2seq, 1, mean)


numModels = length(models)
postprobs_seq = array(NA, dim = c(numModels, numSeq, numSeqMMED))
for(k in 1:numSeqMMED){
  smmed_data_k = smmed_data2_list[[k]]
  for(i in 1:numSeq){
    changing_postprobs = calcEPPHdata(smmed_data_k$y[1:(N_seq * i)], 
                                                       smmed_data_k$D[1:(N_seq * i)], 
                                                       N = N_seq * i, models, sigmasq)
    postprobs_seq[ , i, k] = changing_postprobs
  }
}
postprobs_means = apply(postprobs_seq, c(1,2), mean)

# par(mfrow=c(1,length(models)))
# for(k in 1:length(models)){
#   plot(x = 1:numSeq, y = postprobs_means[ k, ], type = "l", 
#        xlab = "Number of Sequential Steps", ylab = paste("P(H", k - 1, "|Y)", sep = ""), 
#        ylim = range(postprobs[ k, ], exppostprobs_space[k], exppostprobs_dopt1[k], exppostprobs_dopt2[k]))
#   # for(m in 1:numSeqMMED){
#   #   lines(x = 1:numSeq, y = postprobs_seq[k, , m], col=rgb(0, 0, 0, 0.1))
#   # } 
#   # lines(x = 1:numSeq, y = postprobs_means[k, ], col = 1, lty = 2)
#   abline(h = exppostprobs_space[k], col = 3)
#   abline(h = exppostprobs_dopt1[k], col = 4)
#   abline(h = exppostprobs_dopt2[k], col = 5)
#   if(k == 1) legend("topright", c(design_names), lty = c(rep(1, length(design_col))), col = c(design_col), bg = "white")
# }

#

ggdata0 = data.table(
  x = 1:numSeq, 
  Dlinear = rep(exppostprobs_dopt1[1], numSeq), 
  Dquadratic = rep(exppostprobs_dopt2[1], numSeq), 
  SpaceFilling = rep(exppostprobs_space[1], numSeq), 
  SeqMED = postprobs_means[ 1, ], 
  Hypothesis = rep("H0", numSeq)
)
ggdata1 = data.table(
  x = 1:numSeq, 
  Dlinear = rep(exppostprobs_dopt1[2], numSeq), 
  Dquadratic = rep(exppostprobs_dopt2[2], numSeq), 
  SpaceFilling = rep(exppostprobs_space[2], numSeq), 
  SeqMED = postprobs_means[ 2, ], 
  Hypothesis = rep("H1", numSeq)
)
ggdataT = data.table(
  x = 1:numSeq, 
  Dlinear = rep(exppostprobs_dopt1[3], numSeq), 
  Dquadratic = rep(exppostprobs_dopt2[3], numSeq), 
  SpaceFilling = rep(exppostprobs_space[3], numSeq), 
  SeqMED = postprobs_means[ 3, ], 
  Hypothesis = rep("True Hypothesis", numSeq)
)
ggdata = rbind(ggdata0, ggdata1, ggdataT)
ggdata.melted = melt(ggdata, id = c("x", "Hypothesis"), value.name = "epph", 
                     variable.name = "Design")
plt5 = ggplot(ggdata.melted, aes(x = x, y = epph, color = Design, linetype = Design)) +
  facet_wrap(vars(Hypothesis)) + 
  geom_path(size = 2) + 
  scale_linetype_manual(values=c(rep("dashed", 3), "solid")) + 
  geom_point(data = ggdata.melted[x == 10], aes(x = x, y = epph), size = 2) + 
  theme_bw(base_size = 20) + 
  theme(panel.grid.minor = element_blank()) + 
  labs(y = "", x = "Stages")
plt5
# ggsave("epph_d3.pdf",
#        plot = last_plot(),
#        device = "pdf",
#        path = image_path,
#        scale = 1,
#        width = 8,
#        height = 4,
#        units = c("in")
# )
# heght = 3, 4
# ggsave("poster_epph_d3_h4.png",
#        plot = last_plot(),
#        device = "png",
#        path = image_path,
#        scale = 1,
#        width = 13.5,
#        height = 4,
#        units = c("in")
# )
ggarrange(plt4, plt5, ncol = 2, widths = c(3, 3))
# ggsave("poster_seqmed_epph_d3_h3.png",
#        plot = last_plot(),
#        device = "png",
#        path = image_path,
#        scale = 1,
#        width = 13.5,
#        height = 3,
#        units = c("in")
# )

# --- MSE(y-hat) --- #

yhatmse_space = getClosedMSEyhat_seq(x_seq, space_filling, N, betaT, typeT, 
                                     c(0, 0, 0, 0), diag(rep(sigmasq01, 4)), sigmasq, 4)
yhatmse_doptquad = getClosedMSEyhat_seq(x_seq, dopt_quadratic, N, betaT, typeT, 
                                        c(0, 0, 0, 0), diag(rep(sigmasq01, 4)), sigmasq, 4)
yhatmse_smmed = getClosedMSEyhat_seq(x_seq, smmed_data2$D, N, betaT, typeT, 
                                     c(0, 0, 0, 0), diag(rep(sigmasq01, 4)), sigmasq, 4)
yhatmse_doptlin = getClosedMSEyhat_seq(x_seq, dopt_linear, N, betaT, typeT, 
                                       c(0, 0, 0, 0), diag(rep(sigmasq01, 4)), sigmasq, 4)

# par(mfrow = c(1,1))
# ylimarg = range(0, yhatmse_space$MSEyhat, yhatmse_smmed$MSEyhat)
# plot(x_seq, yhatmse_space$MSEyhat, type = "l", col = 3, ylim = ylimarg, 
#      ylab = "", xlab = "")
# lines(x_seq, yhatmse_doptquad$MSEyhat, col = 5)
# lines(x_seq, yhatmse_smmed$MSEyhat, col = 1)
# lines(x_seq, yhatmse_doptlin$MSEyhat, col = 4)
# legend("top", design_names, col = design_col, lty = 1, bg = "white")

ggdata = data.table(
  x = x_seq, 
  Dlinear = yhatmse_doptlin$MSEyhat, 
  Dquadratic = yhatmse_doptquad$MSEyhat, 
  SpaceFilling = yhatmse_space$MSEyhat, 
  SeqMED = yhatmse_smmed$MSEyhat
)
ggdata = melt(ggdata, id = c("x"), value.name = "yhatmse", variable.name = "Design")
ggplot(ggdata, aes(x = x, y = yhatmse, color = Design)) +
  coord_cartesian(ylim = ylimarg, xlim = c(-1, 1)) +
  # scale_y_continuous(limits = ylimarg) + 
  # scale_x_continuous(limits = c(-1, 1)) +
  geom_path() + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  labs(y = "", x = "x")
# ggsave("mseyhat_d3.pdf",
#        plot = last_plot(),
#        device = "pdf",
#        path = image_path,
#        scale = 1,
#        width = 8,
#        height = 4,
#        units = c("in")
# )

#

# --- MSE(Bn) --- #

mu2 = rep(0, 4)
V2 = diag(rep(sigmasq01, length(mu2)))
postvar_doptlin = getClosedMSE(dopt_linear, N, betaT, mu2, V2, sigmasq, typeT)$MSE_postmean
postvar_doptquad = getClosedMSE(dopt_quadratic, N, betaT, mu2, V2, sigmasq, typeT)$MSE_postmean
postvar_space = getClosedMSE(space_filling, N, betaT, mu2, V2, sigmasq, typeT)$MSE_postmean
postvar_smmed = getClosedMSE(smmed_data2$D, N, betaT, mu2, V2, sigmasq, typeT)$MSE_postmean

# par(mfrow = c(1,3))
# barplot(c(postvar_doptlin[1], postvar_doptquad[1], postvar_space[1], postvar_smmed[1]), 
#         names.arg = design_names, main = "B0")
# barplot(c(postvar_doptlin[2], postvar_doptquad[2], postvar_space[2], postvar_smmed[2]), 
#         names.arg = design_names, main = "B1")
# barplot(c(postvar_doptlin[3], postvar_doptquad[3], postvar_space[3], postvar_smmed[3]), 
#         names.arg = design_names, main = "B2")

#

b0 = c(postvar_doptlin[1], postvar_doptquad[1], postvar_space[1], postvar_smmed[1])
b1 = c(postvar_doptlin[2], postvar_doptquad[2], postvar_space[2], postvar_smmed[2])
b2 = c(postvar_doptlin[3], postvar_doptquad[3], postvar_space[3], postvar_smmed[3])
b3 = c(postvar_doptlin[4], postvar_doptquad[4], postvar_space[4], postvar_smmed[4])

ggdata = data.frame(Designs = rep(c("Dlinear", "Dquadratic", "SpaceFilling", "SeqMED"), 4), 
                    MSE = c(b0, b1, b2, b3), beta = rep(c("B0", "B1", "B2", "B3"), each = length(b0)))
ggplot(ggdata, aes(x = Designs, y = MSE)) + 
  geom_bar(stat = "identity") +
  facet_wrap(vars(beta)) +
  theme_bw(base_size = 20) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.text.x = element_text(angle = 45, vjust = 0.5)) +
  labs(y = NULL)
# ggsave("msebeta_d2.pdf",
#        plot = last_plot(),
#        device = "pdf",
#        path = image_path,
#        scale = 1,
#        width = 8,
#        height = 4,
#        units = c("in")
# )

ggdataMSEBn = data.table(
  Designs = c("Dlinear", "Dquadratic", "SpaceFilling", "SeqMED"), 
  MSE = b3
)
ggdata0 = data.table(
  x = 1:numSeq, 
  Dlinear = rep(exppostprobs_dopt1[1], numSeq), 
  Dquadratic = rep(exppostprobs_dopt2[1], numSeq), 
  SpaceFilling = rep(exppostprobs_space[1], numSeq), 
  SeqMED = postprobs_means[ 1, ], 
  Hypothesis = rep("H0", numSeq)
)
ggdata1 = data.table(
  x = 1:numSeq, 
  Dlinear = rep(exppostprobs_dopt1[2], numSeq), 
  Dquadratic = rep(exppostprobs_dopt2[2], numSeq), 
  SpaceFilling = rep(exppostprobs_space[2], numSeq), 
  SeqMED = postprobs_means[ 2, ], 
  Hypothesis = rep("H1", numSeq)
)
ggdataT = data.table(
  x = 1:numSeq, 
  Dlinear = rep(exppostprobs_dopt1[3], numSeq), 
  Dquadratic = rep(exppostprobs_dopt2[3], numSeq), 
  SpaceFilling = rep(exppostprobs_space[3], numSeq), 
  SeqMED = postprobs_means[ 3, ], 
  Hypothesis = rep("True", numSeq)
)
ggdata = rbind(ggdata0, ggdata1, ggdataT)
ggdata.melted = melt(ggdata, id = c("x", "Hypothesis"), value.name = "epph", 
                     variable.name = "Design")
plt1 = ggplot(ggdata.melted, aes(x = x, y = epph, color = Design, linetype = Design)) +
  facet_wrap(vars(Hypothesis)) + 
  geom_path(size = 2) + 
  scale_linetype_manual(values=c(rep("dashed", 3), "solid")) + 
  geom_point(data = ggdata.melted[x == 10], aes(x = x, y = epph), size = 2) + 
  theme_bw(base_size = 20) + 
  theme(panel.grid.minor = element_blank()) + 
  labs(y = "", x = "Stages")

plt2 = ggplot(ggdataMSEBn, aes(x = Designs, y = MSE)) + 
  geom_bar(stat = "identity") +
  theme_bw(base_size = 20) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.text.x = element_text(angle = 45, vjust = 0.5)) +
  labs(y = NULL)
ggarrange(plt2, plt1, widths = c(1, 3))
# 12, 14 in
ggsave("poster_msebn_epph_d3_h3.png",
       plot = last_plot(),
       device = "png",
       path = image_path,
       scale = 1,
       width = 13.5,
       height = 3,
       units = c("in")
)
