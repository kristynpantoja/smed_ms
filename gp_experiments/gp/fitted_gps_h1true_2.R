# for(scenario in c(1, 2)){
scenario = 1
text_size = 12

################################################################################
# last updated: 2/27/2022
# purpose: to illustrate seqmed in scenarios 1 & 2 -- black&white-friendly

################################################################################
# Sources/Libraries
################################################################################
sims_dir = "gp_experiments/gp"
output_dir = paste0(
  sims_dir, "/modelselection_designs/scenarios_h1true/outputs")
data_dir = paste0(sims_dir, "/simulated_data/outputs")
functions_dir = "functions"

# for seqmed design
source(paste(functions_dir, "/SeqMEDgp.R", sep = ""))
source(paste(functions_dir, "/SeqMEDgp_batch.R", sep = ""))
source(paste(functions_dir, "/charge_function_q.R", sep = ""))
source(paste(functions_dir, "/covariance_functions.R", sep = ""))
source(paste(functions_dir, "/wasserstein_distance.R", sep = ""))
source(paste(functions_dir, "/gp_predictive.R", sep = ""))

# for box-hill design
source(paste(functions_dir, "/boxhill.R", sep = ""))
source(paste(functions_dir, "/boxhill_gp.R", sep = ""))
source(paste(functions_dir, "/kl_divergence.R", sep = ""))

library(mvtnorm)
# rng.seed = 123

library(ggplot2)
library(reshape2)
library(data.table)
gg_color_hue = function(n) {
  hues = seq(15, 275, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

################################################################################
# simulation settings, shared for both scenarios
################################################################################

# simulations settings
numSims = 1 ### just generate 1 for illustration
numSeq = 6
seqN = 1
Nttl = numSeq * seqN
xmin = -1
xmax = 1
numx = 10^3 + 1
x_seq = seq(from = xmin, to = xmax, length.out = numx)
sigmasq_measuremt = 1e-10
sigmasq_signal = 1

# set.seed(12345)
set.seed(1234)

################################################################################
# Scenario settings
################################################################################
if(scenario == 1){
  pT = NULL
  l01= c(0.01, 0.01) 
  type01 = c("squaredexponential", "matern")
  model0 = list(type = type01[1], l = l01[1], signal.var = sigmasq_signal, 
                measurement.var = sigmasq_measuremt)
  model1 = list(type = type01[2], l = l01[2], signal.var = sigmasq_signal, 
                measurement.var = sigmasq_measuremt)
  scenario_name = "SMM"
} else if(scenario == 2){
  pT = 0.26
  l01= c(0.3, 0.3) 
  type01 = c("matern", "periodic")
  model0 = list(type = type01[1], l = l01[1], signal.var = sigmasq_signal, 
                measurement.var = sigmasq_measuremt)
  model1 = list(type = type01[2], l = l01[2], signal.var = sigmasq_signal, 
                measurement.var = sigmasq_measuremt, p = pT)
  scenario_name = "MPP"
}
typeT = type01[2]
lT = l01[2]

################################################################################
# generate functions
null_cov = getCov(x_seq, x_seq, typeT, lT, pT, sigmasq_signal)
null_mean = rep(0, numx)

# the function values
y_seq_mat = t(rmvnorm(n = numSims, mean = null_mean, sigma = null_cov)) 
if(is.null(sigmasq_measuremt)){
  y_seq_mat = t(rmvnorm(n = numSims, mean = null_mean, sigma = null_cov))
} else{
  y_seq_mat = t(rmvnorm(n = numSims, mean = null_mean, sigma = null_cov + 
                          sigmasq_measuremt * diag(numx))) 
}

################################################################################
# generate seqmed
y_seq = y_seq_mat[ , 1]
seqmed = SeqMEDgp(
  y.in = NULL, x.in = NULL, x.in.idx = NULL, 
  candidates = x_seq, function.values = y_seq, xmin = xmin, xmax = xmax,
  model0 = model0, model1 = model1, numSeq = numSeq, seqN = seqN)

################################################################################
# make plots
################################################################################

x_train = seqmed$x.in
y_train = seqmed$y.in
newpts = seqmed$x.new
truey = seqmed$y.new
x.all = c(x_train, newpts)
y.all = c(y_train, truey)

H0_predfn = getGPPredictive(
  x_seq, x.all[-Nttl], y.all[-Nttl], model0$type, model0$l, model0$p, 
  model0$signal.var, model0$measurement.var)
H1_predfn = getGPPredictive(
  x_seq, x.all[-Nttl], y.all[-Nttl], model1$type, model1$l, model1$p, 
  model1$signal.var, model1$measurement.var)

# get w_seq
Kinv0 = solve(getCov(
  x.all[-Nttl], x.all[-Nttl], model0$type, model0$l, model0$p, 
  model0$signal.var))
Kinv1 = solve(getCov(
  x.all[-Nttl], x.all[-Nttl], model1$type, model1$l, model1$p, 
  model1$signal.var))
w_seq = sapply(x_seq, FUN = function(x1) 
  WNgp(x1, Kinv0, Kinv1, x.all[-Nttl], y.all[-Nttl], model0, model1))

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
yrange = range(
  c(ggdata$lower0, ggdata$lower1, ggdata$upper0, ggdata$upper1)
  )
yrange[1] = yrange[1] - max(ggdata$Wasserstein)
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
  x = x.all, 
  y = y.all, 
  color = c(rep(gg_color_hue(2)[2], Nttl - 1), gg_color_hue(2)[1]), 
  shape = c(rep(18, Nttl - 1), 16)
)
plt = ggplot(data = ggdata.melted, 
             aes(x = x, y = value, color = variable, linetype = variable)) +#, 
             # linetype = 1) + 
  geom_path() + 
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = variable), 
              alpha = 0.1, linetype = 0) +
  scale_linetype_manual(values = c(1, 1, 2, 3)) + 
  scale_fill_manual(values = c("#FFFFFF", "#FFFFFF", "#00BFC4", "#C77CFF")) + 
  scale_color_manual(values = c(1, "gray", "#00BFC4", "#C77CFF")) + 
  geom_point(data = ggdata_pts, mapping = aes(x = x, y = y), 
             inherit.aes = FALSE, color = ggdata_pts$color, 
             shape = ggdata_pts$shape,
             size = 2) +
  geom_point(data = ggdata_pts, mapping = aes(x = x, y = yrange[1]), 
             inherit.aes = FALSE, color = ggdata_pts$color, 
             shape = ggdata_pts$shape, 
             size = 2, alpha = 0.5) +
  scale_y_continuous(limits = yrange) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        text = element_text(size = text_size)) +
  labs(y = element_blank(), x = element_blank(), 
       fill = "Function", color = "Function", linetype = "Function") + 
  xlim(-0.25, 0.25)
plt

# manuscript plot
ggsave(
  filename = paste0(scenario_name,"_lambda", lT, "_Nttl", Nttl, "_gp_zoom.pdf"),
  plot = plt,
  width = 6, height = 2.5, units = c("in")
)
