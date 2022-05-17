################################################################################
# last updated: 2/27/2022
# purpose: to see what seqmed gp produces as n gets super large

scenario = 2 # 1, 2

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

# set up parallelization
# library(foreach)
# library(future)
# library(doFuture)
# library(parallel)
# registerDoFuture()
# nworkers = detectCores() - 2
# plan(multisession, workers = nworkers)

# library(rngtools)
# library(doRNG)
rng.seed = 123 # 123, 345
# registerDoRNG(rng.seed)
set.seed(123)

library(ggplot2)
library(reshape2)
gg_color_hue = function(n) {
  hues = seq(15, 275, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

library(microbenchmark)

################################################################################
# simulation settings, shared for both scenarios
################################################################################

# simulations settings
numSims = 100
numSeq = 200
seqN = 1
Nttl = numSeq * seqN
xmin = -1
xmax = 1
numx = 10^3 + 1
x_seq = seq(from = xmin, to = xmax, length.out = numx)
sigmasq_measuremt = 1e-10
sigmasq_signal = 1

################################################################################
# Scenario settings
################################################################################
l01= c(0.01, 0.01)
if(scenario == 1){
  type01 = c("squaredexponential", "matern")
  model0 = list(type = type01[1], l = l01[1], signal.var = sigmasq_signal, 
                measurement.var = sigmasq_measuremt)
  model1 = list(type = type01[2], l = l01[2], signal.var = sigmasq_signal, 
                measurement.var = sigmasq_measuremt)
} else if(scenario == 2){
  pT = 0.26
  type01 = c("matern", "periodic")
  model0 = list(type = type01[1], l = l01[1], signal.var = sigmasq_signal, 
                measurement.var = sigmasq_measuremt)
  model1 = list(type = type01[2], l = l01[2], signal.var = sigmasq_signal, 
                measurement.var = sigmasq_measuremt, p = pT)
}
typeT = type01[2]
lT = l01[2]

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
# generate seqmeds 

# simulations!
b = 1

y_seq = y_seq_mat[ , b]
start.time = Sys.time()
seqmed = SeqMEDgp(
  y.in = NULL, x.in = NULL, x.in.idx = NULL,
  candidates = x_seq, function.values = y_seq, xmin = xmin, xmax = xmax,
  model0 = model0, model1 = model1, numSeq = numSeq, seqN = seqN)
end.time = Sys.time()
time.diff = difftime(end.time, start.time, units = "secs")
time.diff

plot(x = c(seqmed$x.in, seqmed$x.new), y = rep(0, 100))

if(scenario == 1){ # 54 min
  # saveRDS(seqmed, file = "seqmedgp_scen1_gvm_n100.rds") 
} else if(scenario == 2){
  if(Nttl = 100){ # 60 min
    # saveRDS(seqmed, file = "seqmedgp_scen2_mpp_n100.rds") 
  } else if(Nttl = 200){
    saveRDS(seqmed, file = "seqmedgp_scen2_mpp_n200.rds") 
  }
}



