# for(scenario in c(3, 4, 5, 6)){
scenario = 5
################################################################################
# last updated: 2/28/2022
# purpose: to test seqmedgp for scenario 3:
#   squared exponential vs. another squared exponential,
#   where the true function is matern

Ntest = 51

################################################################################
# Sources/Libraries
################################################################################
sims_dir = "gp_experiments/gp"
output_dir = paste0(
  sims_dir, "/modelselection_designs/scenarios_misspecified/outputs")
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
rng.seed = 123

library(ggplot2)
library(reshape2)
gg_color_hue = function(n) {
  hues = seq(15, 275, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

################################################################################
# simulation settings, shared for both scenarios
################################################################################

# simulations settings
numSims = 100
numSeq = 15
seqN = 1
Nttl = numSeq * seqN
xmin = -1
xmax = 1
numx = 10^3 + 1
x_seq = seq(from = xmin, to = xmax, length.out = numx)
sigmasq_measuremt = 1e-10
sigmasq_signal = 1

# shared settings
prior_probs = rep(1 / 2, 2)

################################################################################
# Scenario settings
################################################################################
if(scenario == 3){
  type01 = c("squaredexponential", "squaredexponential")
  typeT = "matern"
  l01= c(0.005, 0.01)
  lT = 0.01
  model0 = list(type = type01[1], l = l01[1], signal.var = sigmasq_signal, 
                measurement.var = sigmasq_measuremt)
  model1 = list(type = type01[2], l = l01[2], signal.var = sigmasq_signal, 
                measurement.var = sigmasq_measuremt)
  scenario_name = "SSM"
} else if(scenario == 4){
  type01 = c("matern", "squaredexponential")
  typeT = "periodic"
  l01= c(0.01, 0.01)
  lT = 0.5
  pT = 0.05 # 0.05 or 0.1
  model0 = list(type = type01[1], l = l01[1], signal.var = sigmasq_signal, 
                measurement.var = sigmasq_measuremt)
  model1 = list(type = type01[2], l = l01[2], signal.var = sigmasq_signal, 
                measurement.var = sigmasq_measuremt)
  scenario_name = "MSP"
} else if(scenario == 5){
  type01 = c("matern", "periodic")
  typeT = "squaredexponential"
  l01= c(0.01, 0.01)
  lT = 0.01
  p1 = 0.26
  model0 = list(type = type01[1], l = l01[1], signal.var = sigmasq_signal, 
                measurement.var = sigmasq_measuremt)
  model1 = list(type = type01[2], l = l01[2], signal.var = sigmasq_signal, 
                measurement.var = sigmasq_measuremt, p = p1)
  scenario_name = "MPS"
} else if(scenario == 6){
  type01 = c("squaredexponential", "periodic")
  typeT = "matern"
  l01= c(0.01, 0.01)
  lT = 0.01
  p1 = 0.26
  model0 = list(type = type01[1], l = l01[1], signal.var = sigmasq_signal, 
                measurement.var = sigmasq_measuremt)
  model1 = list(type = type01[2], l = l01[2], signal.var = sigmasq_signal, 
                measurement.var = sigmasq_measuremt, p = p1)
  scenario_name = "SPM"
} else{
  stop("invalid scenario number")
}

################################################################################
# import data
filename_append = ""
if(!is.null(sigmasq_measuremt)){
  filename_append = "_noise"
}
if(typeT == "periodic"){
  simulated_data_file = paste0(
    data_dir,
    "/", typeT,
    "_l", lT,
    "_p", pT,
    filename_append, 
    "_seed", rng.seed,
    ".rds")
} else{
  simulated_data_file = paste0(
    data_dir,
    "/", typeT,
    "_l", lT,
    filename_append, 
    "_seed", rng.seed,
    ".rds")
}
simulated.data = readRDS(simulated_data_file)
numSims = simulated.data$numSims
x_seq = simulated.data$x
numx = length(x_seq)
null_cov = simulated.data$null_cov
null_mean = simulated.data$null_mean
y_seq_mat = simulated.data$function_values_mat

################################################################################
# initial design
x_input_idx = ceiling(numx / 2)
x_input = x_seq[x_input_idx]

################################################################################
# read in the data

# filename_append.tmp for all methods alike
filename_append.tmp = paste0(
  filename_append, 
  "_seed", rng.seed,
  ".rds"
)
boxhill_sims = readRDS(paste0(
  output_dir,
  "/scenario", scenario, "_boxhill", 
  filename_append.tmp))
leaveout_sims = readRDS(paste0(
  output_dir,
  "/scenario", scenario, "_seqmed", 
  filename_append.tmp))

if(typeT == "periodic"){
  random_sims_file = paste0(
    sims_dir, 
    "/spacefilling_designs/outputs/random", 
    "_", typeT,
    "_l", lT,
    "_p", pT,
    filename_append.tmp)
  grid_sims_file = paste0(
    sims_dir,
    "/spacefilling_designs/outputs/grid", 
    "_", typeT,
    "_l", lT,
    "_p", pT,
    filename_append.tmp)
} else{
  random_sims_file = paste0(
    sims_dir, 
    "/spacefilling_designs/outputs/random", 
    "_", typeT,
    "_l", lT,
    filename_append.tmp)
  grid_sims_file = paste0(
    sims_dir,
    "/spacefilling_designs/outputs/grid", 
    "_", typeT,
    "_l", lT,
    filename_append.tmp)
}
random_sims = readRDS(random_sims_file)
grid_sims = readRDS(grid_sims_file)

################################################################################
# make RSS plots for true model
################################################################################

# models
if(typeT == "periodic"){
  modelT = list(type = typeT, l = lT, signal.var = sigmasq_signal, 
                measurement.var = sigmasq_measuremt, p = pT)
} else{
  modelT = list(type = typeT, l = lT, signal.var = sigmasq_signal, 
                measurement.var = sigmasq_measuremt)
}

getRSS = function(
  design, model, candidates, function.values, Ntest = 51
){
  # grid of Ntest points
  x.test.idx = 1 + 
    0:(Ntest - 1) * floor(length(candidates) / (Ntest - 1))
  x.test = candidates[x.test.idx]
  y.test = function.values[x.test.idx]
  
  # get yhat.test, after fitting model to given design data
  x.tmp = as.vector(na.omit(c(design$x.in, design$x.new)))
  y.tmp = as.vector(na.omit(c(design$y.in, design$y.new)))
  yhat.test = getGPPredictive(
    x.test, x.tmp, y.tmp, model$type, model$l, model$p, model$signal.var, 
    model$measurement.var)
  RSS.tmp = sum((yhat.test$pred_mean - y.test)^2)#, na.rm = TRUE)
  return(data.frame(RSS = RSS.tmp))
}

RSS.df = data.frame(
  RSS = numeric(), type = character(), sim = numeric())
for(j in 1:numSims){
  # designs at sim b
  bh = boxhill_sims[[j]]
  lo = leaveout_sims[[j]]
  r = random_sims[[j]]
  g = grid_sims[[j]]
  # sequence of PPHs for each design
  RSS.bh = getRSS(bh, modelT, x_seq, y_seq_mat[, j], Ntest)
  RSS.lo = getRSS(lo, modelT, x_seq, y_seq_mat[, j], Ntest)
  RSS.r = getRSS(r, modelT, x_seq, y_seq_mat[, j], Ntest)
  RSS.g = getRSS(g, modelT, x_seq, y_seq_mat[, j], Ntest)
  # master data frame
  RSS.bh$type = "boxhill"
  RSS.lo$type = "seqmed" # "leaveout"
  RSS.r$type = "random"
  RSS.g$type = "grid"
  # RSS.tmp = rbind(
  #   RSS.bh, RSS.qc, RSS.lo, RSS.qc2, RSS.kq, 
  #   RSS.r, RSS.g)
  RSS.tmp = rbind(
    RSS.bh, RSS.lo, #RSS.qc, RSS.lo, RSS.qc2, RSS.kq, 
    RSS.r, RSS.g)
  RSS.tmp$sim = j
  RSS = rbind(RSS.df, RSS.tmp)
  RSS.tmp$sim = j
  RSS.df = rbind(RSS.df, RSS.tmp)
}
RSSmean = aggregate(
  RSS.df$RSS, by = list(RSS.df$type), 
  FUN = function(x) mean(x, na.rm = TRUE))
names(RSSmean) = c("type", "value")

RSS.plt = ggplot(RSSmean, aes(x = type, y = value)) + 
  geom_bar(stat = "identity") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  labs(y = "E[RSS]", x = "Design") 
plot(RSS.plt)

print(RSSmean)

# slide plot
# ggsave(
#   filename = paste0("20210902_scen", scenario, "_rsst.pdf"), 
#   plot = RSS.plt, 
#   width = 6, height = 4, units = c("in")
# )

# manuscript plot
ggsave(
  filename = paste0(scenario_name, "_rsst.pdf"), 
  plot = RSS.plt, 
  width = 4.5, height = 2, units = c("in")
)

# print(paste("scenario", scenario, 
#             "################################################################"))
# }