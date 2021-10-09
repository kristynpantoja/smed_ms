for(scenario in c(1, 2)){
  ################################################################################
  # last updated: 05/25/2021
  # purpose: to test seqmedgp for scenario 1:
  #   squared exponential vs. matern,
  #   where the true function is matern
  
  # scenario = 1 # scenarios: 1, 2
  # input.type = 1 # 1 = extrapolation, 2 = inc spread, 3 = even coverage
  seq.type = 1 # 1 = fully sequential, 2 = stage-sequential 3x5
  
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
  
  # all of the designs
  idx = 1
  # designs = list(
  #   boxhill_sims[[idx]], qcap_sims[[idx]], leaveout_sims[[idx]], 
  #   qcap_persist_sims[[idx]], persist_sims[[idx]])
  # design.names = c("boxhill", "qcap", "leaveout", "keepq2", "keepq")
  # design.levels = c("leaveout", "keepq", "keepq2", "qcap", "boxhill")
  designs = list(
    "BoxHill" = boxhill_sims[[idx]], 
    "Grid" = grid_sims[[idx]], 
    "Random" = random_sims[[idx]],
    "SeqMED" = leaveout_sims[[idx]])
  design.names = names(designs)
  
  x.new.mat = matrix(NA, nrow = Nnew, ncol = length(designs))
  for(i in 1:length(designs)){
    x.new.mat[, i] = designs[[i]]$x.new
  }
  colnames(x.new.mat) = names(designs)
  data.gg = as.data.frame(x.new.mat)
  data.gg$index = as.character(1:Nnew)
  data.gg = reshape2::melt(data.gg, id.vars = "index", variable.name = "Design")
  data.gg$Design = factor(data.gg$Design, levels = design.names)
  
  data.gg0 = data.frame(
    Design = rep(design.names, each = Nin), 
    input = rep(x_input, length(designs))
  )
  text.gg = dplyr::filter(data.gg, index %in% as.character(1:Nnew))
  des.plt = ggplot() + 
    geom_point(data = data.gg0, 
               mapping = aes(x = input, y = Design), alpha = 0.25) +
    geom_point(data = data.gg, 
               mapping = aes(x = value, y = Design, color = Design), 
               inherit.aes = FALSE, alpha = 0.25) + 
    # geom_text(data = text.gg, 
    #           aes(x = value, y = Design, label = index), 
    #           vjust = -0.65 * as.numeric(paste(text.gg$index)), size = 2) +
    xlim(c(xmin, xmax)) + 
    theme_bw() + 
    labs(x = "x") +
    theme(legend.position = "none")
  if(typeT == "periodic"){
    des.plt = des.plt + geom_vline(
      xintercept = pT * (0:floor((xmax - xmin) / pT)), #color = "gray", 
      alpha = 0.125, linetype = 2)
  }
  plot(des.plt)
  
  # slide plot
  # ggsave(
  #   filename = paste0("20210815_scen", scenario, "_design.pdf"), 
  #   plot = des.plt, 
  #   width = 6, height = 4, units = c("in")
  # )
  
  # manuscript plot
  ggsave(
    filename = paste0(scenario_name, "_design.pdf"), 
    plot = des.plt, 
    width = 6.5, height = 1.75, units = c("in")
  )
  
  print(paste("scenario", scenario, 
              "################################################################"))
  
}