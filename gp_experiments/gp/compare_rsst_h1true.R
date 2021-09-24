for(scenario in c(1, 2)){
  ################################################################################
  # last updated: 09/02/2021
  # purpose: to test seqmedgp for scenario 1:
  #   squared exponential vs. matern,
  #   where the true function is matern
  
  # scenario = 1 # scenarios: 1, 2
  # input.type = 1 # 1 = extrapolation, 2 = inc spread, 3 = even coverage
  seq.type = 1 # 1 = fully sequential, 2 = stage-sequential 3x5
  
  Ntest = 51
  
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
  Nin = 1
  if(seq.type == 1){
    numSeq = 15
    seqN = 1
  } else if(seq.type == 2){
    numSeq = 3
    seqN = 5
  }
  Nnew = numSeq * seqN
  Nttl = Nin + Nnew
  xmin = 0
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
  l01= c(0.01, 0.01)
  if(scenario == 1){
    type01 = c("squaredexponential", "matern")
    model0 = list(type = type01[1], l = l01[1], signal.var = sigmasq_signal, 
                  measurement.var = sigmasq_measuremt)
    model1 = list(type = type01[2], l = l01[2], signal.var = sigmasq_signal, 
                  measurement.var = sigmasq_measuremt)
    scenario_name = "SMM"
  } else if(scenario == 2){
    pT = 0.26
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
    "_leaveout", 
    "_seq", seq.type,
    filename_append.tmp))
  qcap_sims = readRDS(paste0(
    output_dir,
    "/scenario", scenario, "_seqmed", 
    "_cap",
    "_seq", seq.type,
    filename_append.tmp))
  persist_sims = readRDS(paste0(
    output_dir,
    "/scenario", scenario, "_seqmed", 
    "_persist", 
    "_seq", seq.type,
    filename_append.tmp))
  qcap_persist_sims = readRDS(paste0(
    output_dir,
    "/scenario", scenario, "_seqmed", 
    "_cap_persist",
    "_seq", seq.type,
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
  # make plots
  ################################################################################
  
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
    qc = qcap_sims[[j]]
    lo = leaveout_sims[[j]]
    qc2 = qcap_persist_sims[[j]]
    kq = persist_sims[[j]]
    r = random_sims[[j]]
    g = grid_sims[[j]]
    # sequence of PPHs for each design
    RSS.bh = getRSS(bh, model1, x_seq, y_seq_mat[, j], Ntest)
    RSS.qc = getRSS(qc, model1, x_seq, y_seq_mat[, j], Ntest)
    RSS.lo = getRSS(lo, model1, x_seq, y_seq_mat[, j], Ntest)
    RSS.qc2 = getRSS(qc2, model1, x_seq, y_seq_mat[, j], Ntest)
    RSS.kq = getRSS(kq, model1, x_seq, y_seq_mat[, j], Ntest)
    RSS.r = getRSS(r, model1, x_seq, y_seq_mat[, j], Ntest)
    RSS.g = getRSS(g, model1, x_seq, y_seq_mat[, j], Ntest)
    # master data frame
    RSS.bh$type = "boxhill"
    RSS.qc$type = "qcap"
    RSS.lo$type = "seqmed" # "leaveout"
    RSS.qc2$type = "keepq2"
    RSS.kq$type = "keepq"
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
  
  print(paste("scenario", scenario, 
              "################################################################"))
  
}