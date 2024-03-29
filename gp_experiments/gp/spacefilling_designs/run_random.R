for(scenario in c(1, 2, 4, 5)){
  ################################################################################
  # last updated: 05/27/2021
  # purpose: to make random design for all types of data
  
  ################################################################################
  # Sources/Libraries
  ################################################################################
  sims_dir = "gp_experiments/gp"
  output_dir = paste0(sims_dir, "/spacefilling_designs/outputs")
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
  
  ################################################################################
  # Scenario settings
  ################################################################################
  
  if(scenario %in% c(1, 3, 6)){
    typeT = "matern"
    lT = 0.01
    pT = NULL
  } else if(scenario == 2){
    typeT = "periodic"
    lT = 0.01
    pT = 0.26
  } else if(scenario == 4){
    typeT = "periodic"
    lT = 0.5
    pT = 0.05
  } else if(scenario == 5){
    typeT = "squaredexponential"
    lT = 0.01
    pT = NULL
  } else{
    stop("invalid scenario given")
  }
  
  ################################################################################
  # import data
  filename_append = ""
  if(!is.null(sigmasq_measuremt)){
    filename_append = "_noise"
  }
  if(typeT == "periodic"){
    if(is.null(pT)) pT = 0.26
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
  # generate random designs (sample x ~ U[xmin, xmax])
  
  # simulations!
  registerDoRNG(rng.seed)
  randoms = foreach(
    i = 1:numSims
  ) %dorng% {
    y_seq = y_seq_mat[ , i]
    
    x.new.idx  = sample(1:numx, Nttl)
    x.new = x_seq[x.new.idx]
    y.new = y_seq[x.new.idx]
    list(x.in = x.new[1], x.in.idx = x.new.idx[1], y.in = y.new[1], 
         x.new = x.new[-1], x.new.idx = x.new.idx[-1], y.new = y.new[-1], 
         function.values = y_seq)
  }
  
  filename_append.tmp = filename_append
  filename_append.tmp = paste0(
    filename_append.tmp, 
    "_seed", rng.seed,
    ".rds"
  )
  if(typeT == "periodic"){
    simulated_spacefilling_file = paste0(
      output_dir,
      "/random", 
      "_", typeT,
      "_l", lT,
      "_p", pT,
      filename_append.tmp)
  } else{
    simulated_spacefilling_file = paste0(
      output_dir,
      "/random", 
      "_", typeT,
      "_l", lT,
      filename_append.tmp)
  }
  saveRDS(randoms, file = simulated_spacefilling_file)
}