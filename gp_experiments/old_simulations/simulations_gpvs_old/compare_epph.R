################################################################################
# last updated: 07/13/2021
# purpose: to test SeqMEDgpvs()

lT = 0.01
typeT = "squaredexponential"
dimT = 2
seq.type = 1


################################################################################
# Sources/Libraries
################################################################################
sims_dir = "gp_experiments/simulations_gpvs"
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
Nin = 3
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

sigmasq_measuremt = 1e-10
sigmasq_signal = 1

# shared settings
nugget = sigmasq_measuremt
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
  "_cap_keepq", 
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

Evidence_gpvs = function(y, x, model){
  null_mean_vec = rep(0, length(y))
  if(is.null(model$measurement.var)){
    K_obs = getCov(
      X1 = x[, model$indices, drop = FALSE], 
      X2 = x[, model$indices, drop = FALSE], type = model$type, l = model$l, 
      p = model$p, signal.var = model$signal.var)
  } else{
    K_obs = getCov(
      X1 = x[, model$indices, drop = FALSE], 
      X2 = x[, model$indices, drop = FALSE], type = model$type, l = model$l, 
      p = model$p, signal.var = model$signal.var) + 
      model$measurement.var * diag(length(y))
  }
  evidence = dmvnorm(
    y, mean = null_mean_vec, sigma = K_obs, log = FALSE)
  return(evidence)
}

getPPHseq = function(design, model0, model1){
  len.tmp = length(as.vector(na.omit(design$y.new)))
  PPH0_seq = rep(NA, len.tmp)
  PPH1_seq = rep(NA, len.tmp)
  for(i in 1:len.tmp){
    y.tmp = c(design$y, design$y.new[1:i])
    x.tmp = rbind(design$x, design$x.new[1:i, , drop = FALSE])
    PPHs.tmp = getHypothesesPosteriors(
      prior.probs = prior_probs, 
      evidences = c(
        Evidence_gpvs(y.tmp, x.tmp, model0),
        Evidence_gpvs(y.tmp, x.tmp, model1)
      )
    )
    PPH0_seq[i] = PPHs.tmp[1]
    PPH1_seq[i] = PPHs.tmp[2]
  }
  if(length(PPH0_seq) < Nnew){
    PPH0_seq[(length(PPH0_seq) + 1):Nnew] = PPH0_seq[length(PPH0_seq)]
  }
  if(length(PPH1_seq) < Nnew){
    PPH1_seq[(length(PPH1_seq) + 1):Nnew] = PPH1_seq[length(PPH1_seq)]
  }
  return(data.frame(
    index = 1:Nnew, 
    PPH0 = PPH0_seq, 
    PPH1 = PPH1_seq
  ))
}

PPH_seq = data.frame(
  PPH0 = numeric(), PPH1 = numeric(), PPHT = numeric(), 
  type = character(), sim = numeric())
for(j in 1:numSims){
  
  # plot
  # fields::quilt.plot(x_grid, y_seq_mat[ , j])
  # plot(x_grid[ , 1], y_seq_mat[ , j])
  # order_truef = order(x_grid[,1])
  # lines(x_grid[order_truef,1], y_seq_mat[order_truef, j])
  
  # designs at sim b
  qc = cap_sims[[j]]
  kq = keepq_sims[[j]]
  kqc = keepq_cap_sims[[j]]
  lo = leaveout_sims[[j]]
  diag = diagonal_sims[[j]]
  x2 = x2_sims[[j]]
  r = random_sims[[j]]
  g = grid_sims[[j]]
  # sequence of PPHs for each design
  PPH_seq.qc = getPPHseq(qc, model0, model1)
  PPH_seq.kq = getPPHseq(kq, model0, model1)
  PPH_seq.kqc = getPPHseq(kqc, model0, model1)
  PPH_seq.lo = getPPHseq(lo, model0, model1)
  PPH_seq.diag = getPPHseq(diag, model0, model1)
  PPH_seq.x2 = getPPHseq(x2, model0, model1)
  PPH_seq.r = getPPHseq(r, model0, model1)
  PPH_seq.g = getPPHseq(g, model0, model1)
  # master data frame
  PPH_seq.qc$type = "cap q"
  PPH_seq.kq$type = "keepq"
  PPH_seq.kqc$type = "keepq2"
  PPH_seq.lo$type = "leaveout"
  PPH_seq.diag$type = "diag"
  PPH_seq.x2$type = "x2=1"
  PPH_seq.r$type = "random"
  PPH_seq.g$type = "grid"
  PPH_seq.tmp = rbind(
    PPH_seq.qc, PPH_seq.kq, PPH_seq.kqc, PPH_seq.lo, PPH_seq.diag, PPH_seq.x2, 
    PPH_seq.r, PPH_seq.g)
  PPH_seq.tmp$sim = j
  PPH_seq = rbind(PPH_seq, PPH_seq.tmp)
}

PPH0mean_seq = aggregate(PPH_seq$PPH0, by = list(PPH_seq$index, PPH_seq$type), 
                         FUN = function(x) mean(x, na.rm = TRUE))
names(PPH0mean_seq) = c("index", "type", "value")
PPH0mean_seq$Hypothesis = "H0"
PPH1mean_seq = aggregate(PPH_seq$PPH1, by = list(PPH_seq$index, PPH_seq$type), 
                         FUN = function(x) mean(x, na.rm = TRUE))
names(PPH1mean_seq) = c("index", "type", "value")
PPH1mean_seq$Hypothesis = "H1"

PPHmean_seq = rbind(PPH0mean_seq, PPH1mean_seq)
epph.plt = ggplot(PPHmean_seq, aes(x = index, y = value, color = type, 
                                   linetype = type, shape = type)) + 
  facet_wrap(~Hypothesis) + 
  geom_path() + 
  geom_point() +
  theme_bw() +
  ylim(0, 1)
plot(epph.plt)







