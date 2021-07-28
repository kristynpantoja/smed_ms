################################################################################
# last updated: 07/13/2021
# purpose: to test SeqMEDgpvs()

dimT = 1
seq.type = 1
lT = 0.1
typeT = "squaredexponential"

################################################################################
# Sources/Libraries
################################################################################
sims_dir = "gp_experiments/simulations_final/simulations_gpvs"
output_dir = paste0(sims_dir, "/scenarios/outputs")
data_dir = paste0(sims_dir, "/simulated_data")
functions_dir = "functions"

# for seqmed design
source(paste(functions_dir, "/SeqMEDgpvs.R", sep = ""))
source(paste(functions_dir, "/SeqMEDgpvs_batch.R", sep = ""))
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

# set up parallelization
library(foreach)
library(future)
library(doFuture)
library(parallel)
registerDoFuture()
nworkers = detectCores()
plan(multisession, workers = nworkers)

library(rngtools)
library(doRNG)
rng.seed = 123 # 123, 345
registerDoRNG(rng.seed)

library(ggplot2)
library(reshape2)
gg_color_hue = function(n) {
  hues = seq(15, 275, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

# library(expm)
# library(matrixStats)
# library(MASS)
# library(mvtnorm)
# library(knitr)

################################################################################
# simulation settings, shared for both scenarios
################################################################################

# simulations settings
# numSims = 25
Nin = 1
if(seq.type == 1){
  numSeq = 9
  seqN = 1
} else if(seq.type == 2){
  numSeq = 3
  seqN = 3
}
Nnew = numSeq * seqN
Nttl = Nin + Nnew 
xmin = 0
xmax = 1

# numx = 21
# x_seq = seq(from = xmin, to = xmax, length.out = numx)
# x_grid = expand.grid(x_seq, x_seq)

p = 2
k = 4 * p

sigmasq_measuremt = NULL
sigmasq_signal = 1

# shared settings
nugget = 1e-10
prior_probs = rep(1 / 2, 2)


################################################################################
# Scenario settings
################################################################################

l01= c(lT, lT)
type01 = c(typeT, typeT)
indices0 = c(1)
indices1 = c(1, 2)

model0 = list(type = type01[1], l = l01[1], signal.var = sigmasq_signal, 
              indices = indices0,
              measurement.var = nugget)
model1 = list(type = type01[2], l = l01[2], signal.var = sigmasq_signal, 
              indices = indices1, 
              measurement.var = nugget)

################################################################################
# import data
filename_append = ""
if(!is.null(sigmasq_measuremt)){
  filename_append = "_noise"
}
simulated_data_file = paste0(
  data_dir,
  "/", typeT,
  "_l", lT,
  "_dim", dimT,
  filename_append, 
  "_seed", rng.seed,
  ".rds")
simulated.data = readRDS(simulated_data_file)
numSims = simulated.data$numSims
x_seq = simulated.data$x_seq
x_grid = as.matrix(simulated.data$x_grid)
numx = length(x_seq)
null_cov = simulated.data$null_cov
null_mean = simulated.data$null_mean
y_seq_mat = simulated.data$function_values_mat

################################################################################
# initial design
x_input_idx = sample(1:numx, Nin)
x_input = x_grid[x_input_idx + numx * (1:Nin), , drop = FALSE]

################################################################################
# read in the designs

# filename_append.tmp for all methods alike
filename_append.tmp = paste0(
  filename_append, 
  "_seed", rng.seed,
  ".rds"
)

leaveout_sims = readRDS(paste0(
  output_dir,
  "/dim", dimT, "_seqmed", 
  "_leaveout", 
  "_seq", seq.type,
  filename_append.tmp))
cap_sims = readRDS(paste0(
  output_dir,
  "/dim", dimT, "_seqmed", 
  "_cap", 
  "_seq", seq.type,
  filename_append.tmp))
keepq_cap_sims = readRDS(paste0(
  output_dir,
  "/dim", dimT, "_seqmed", 
  "_keepq_cap", 
  "_seq", seq.type,
  filename_append.tmp))
keepq_sims = readRDS(paste0(
  output_dir,
  "/dim", dimT, "_seqmed", 
  "_keepq", 
  "_seq", seq.type,
  filename_append.tmp))

random_sims = readRDS(paste0(
  sims_dir, 
  "/fixed_designs/outputs/random", 
  "_", typeT,
  "_l", lT,
  "_dim", dimT,
  filename_append.tmp))
grid_sims = readRDS(paste0(
  sims_dir,
  "/fixed_designs/outputs/grid", 
  "_", typeT,
  "_l", lT,
  "_dim", dimT,
  filename_append.tmp))
diagonal_sims = readRDS(paste0(
  sims_dir,
  "/fixed_designs/outputs/diagonal", 
  "_", typeT,
  "_l", lT,
  "_dim", dimT,
  filename_append.tmp))
x2_sims = readRDS(paste0(
  sims_dir,
  "/fixed_designs/outputs/x2constant", 
  "_", typeT,
  "_l", lT,
  "_dim", dimT,
  filename_append.tmp))



################################################################################
# compare it to space-filling designs

sim.idx = 1

designs = list(
  cap_sims[[sim.idx]], keepq_sims[[sim.idx]], keepq_cap_sims[[sim.idx]], 
  leaveout_sims[[sim.idx]], 
  diagonal_sims[[sim.idx]], x2_sims[[sim.idx]], random_sims[[sim.idx]], 
  grid_sims[[sim.idx]]
  )
design_names = c(
  "cap q", "keepq", "keepq2", "leaveout", "diag", "x2=1", "random", "grid"
)
for(j in 1:length(designs)){
  des = designs[[j]]
  
    # fields::quilt.plot(x_grid, y_seq_mat[ , sim.idx], main = design_names[j], )
    plot(x_grid[, 1], x_grid[, 2], col = "lightgray", pch = 16, cex = 1, 
           main = design_names[j])
    points(des$x[, 1], des$x[, 2], col = "darkgreen", pch = "o", cex = 1.5)
    points(des$x.new[, 1], des$x.new[, 2], col = "purple", pch = "o", cex = 1.5)
}


if(dimT == 2){
  fields::quilt.plot(x_grid, y_seq_mat[ , sim.idx])
}
if(dimT == 1){
  fields::quilt.plot(x_grid, y_seq_mat[ , sim.idx])
  plot(x_grid[ , 1], y_seq_mat[ , sim.idx])
  order_truef = order(x_grid[,1])
  lines(x_grid[order_truef,1], y_seq_mat[order_truef, sim.idx])
}



