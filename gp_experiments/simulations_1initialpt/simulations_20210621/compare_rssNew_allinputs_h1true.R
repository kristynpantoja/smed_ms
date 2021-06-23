for(scenario in c(1, 2)){
  ################################################################################
  # last updated: 05/25/2021
  # purpose: to test seqmedgp for scenario 1:
  #   squared exponential vs. matern,
  #   where the true function is matern
  
  # scenario = 1 # scenarios: 1, 2
  seq.type = 1 # 1 = fully sequential, 2 = stage-sequential 3x5
  
  ################################################################################
  # Sources/Libraries
  ################################################################################
  sims_dir = "gp_experiments/simulations_1initialpt"
  modelsel_sims_dir = paste0(sims_dir, "/simulations_20210621")
  output_home = paste0(modelsel_sims_dir, "/scenarios_h1true/outputs")
  data_home = "gp_experiments/simulated_data"
  functions_home = "functions"
  
  # for seqmed design
  source(paste(functions_home, "/SeqMEDgp.R", sep = ""))
  source(paste(functions_home, "/SeqMEDgp_batch.R", sep = ""))
  source(paste(functions_home, "/charge_function_q.R", sep = ""))
  source(paste(functions_home, "/covariance_functions.R", sep = ""))
  source(paste(functions_home, "/wasserstein_distance.R", sep = ""))
  source(paste(functions_home, "/gp_predictive.R", sep = ""))
  
  # for box-hill design
  source(paste(functions_home, "/boxhill.R", sep = ""))
  source(paste(functions_home, "/boxhill_gp.R", sep = ""))
  source(paste(functions_home, "/kl_divergence.R", sep = ""))
  
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
  numSims = 25
  Nin = 6
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
  nugget = sigmasq_measuremt
  prior_probs = rep(1 / 2, 2)
  
  ################################################################################
  # Scenario settings
  ################################################################################
  if(scenario == 1){
    type01 = c("squaredexponential", "matern")
  } else if(scenario == 2){
    type01 = c("matern", "periodic")
  }
  typeT = type01[2]
  l01= c(0.01, 0.01)
  lT = l01[2]
  
  ################################################################################
  # import data
  filename_append = ""
  if(!is.null(sigmasq_measuremt)){
    filename_append = "_noise"
  }
  simulated.data = readRDS(paste0(
    data_home,
    "/", typeT,
    "_l", lT,
    filename_append, 
    "_seed", rng.seed,
    ".rds"))
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
  
  boxhills = list()
  # qs = list()
  # q1s = list()
  qcaps = list()
  buffers = list()
  qcap_persists = list()
  buffer_persists = list()
  randoms = list()
  spacefills = list()
  
  for(i in 1:3){
    # filename_append.tmp for all methods alike
    filename_append.tmp = paste0(
      filename_append, 
      "_input", i, 
      "_seed", rng.seed,
      ".rds"
    )
    boxhills[[i]] = readRDS(paste0(
      output_home,
      "/scenario", scenario, "_boxhill", 
      filename_append.tmp))
    buffers[[i]] = readRDS(paste0(
      output_home,
      "/scenario", scenario, "_seqmed", 
      "_buffer", 
      "_seq", seq.type,
      filename_append.tmp))
    # qs[[i]] = readRDS(paste0(
    #   output_home,
    #   "/scenario", scenario, "_seqmed", 
    #   "_q",
    #   "_seq", seq.type,
    #   filename_append.tmp))
    # q1s[[i]] = readRDS(paste0(
    #   output_home,
    #   "/scenario", scenario, "_seqmed", 
    #   "_uniform",
    #   "_seq", seq.type,
    #   filename_append.tmp))
    qcaps[[i]] = readRDS(paste0(
      output_home,
      "/scenario", scenario, "_seqmed", 
      "_cap",
      "_seq", seq.type,
      filename_append.tmp))
    buffer_persists[[i]] = readRDS(paste0(
      output_home,
      "/scenario", scenario, "_seqmed", 
      "_buffer_persist", 
      "_seq", seq.type,
      filename_append.tmp))
    qcap_persists[[i]] = readRDS(paste0(
      output_home,
      "/scenario", scenario, "_seqmed", 
      "_cap_persist",
      "_seq", seq.type,
      filename_append.tmp))
    
    randoms[[i]] = readRDS(paste0(
      "gp_experiments/spacefilling_designs/outputs/random", 
      "_", typeT,
      "_l", lT,
      filename_append.tmp))
    spacefills[[i]] = readRDS(paste0(
      "gp_experiments/spacefilling_designs/outputs/grid", 
      "_", typeT,
      "_l", lT,
      filename_append.tmp))
  }
  
  ################################################################################
  # make plots
  ################################################################################
  
  # models
  model0 = list(type = type01[1], l = l01[1], signal.var = sigmasq_signal, 
                measurement.var = nugget)
  model1 = list(type = type01[2], l = l01[2], signal.var = sigmasq_signal, 
                measurement.var = nugget)
  
  # calculate the RSSnew
  getRSS01 = function(
    design, model0, model1, candidates, function.values, Nother = 51
  ){
    x.other.idx = 1 + 0:(Nother - 1) * floor(length(candidates) / (Nother - 1))
    x.other = candidates[x.other.idx]
    y.other = function.values[x.other.idx]
    
    x.tmp = as.vector(na.omit(c(design$x, design$x.new)))
    y.tmp = as.vector(na.omit(c(design$y, design$y.new)))
    pred0.tmp = getGPPredictive(x.other, x.tmp, y.tmp, 
                                model0$type, model0$l, 
                                sigmasq_signal, model0$measurement.var)
    pred1.tmp = getGPPredictive(x.other, x.tmp, y.tmp, 
                                model1$type, model1$l, 
                                sigmasq_signal, model1$measurement.var)
    RSS0.tmp = sum((pred0.tmp$pred_mean - y.other)^2, na.rm = TRUE)  
    RSS1.tmp = sum((pred1.tmp$pred_mean - y.other)^2, na.rm = TRUE)
    RSS01.tmp = RSS0.tmp / RSS1.tmp
    return(data.frame("RSS0" = RSS0.tmp, "RSS1" = RSS1.tmp, "RSS01" = RSS01.tmp))
  }
  
  RSS.df = data.frame(RSS0 = numeric(), RSS1 = numeric(), RSS01 = numeric(), 
                      type = character(), sim = numeric(), input = numeric())
  for(k in 1:3){
    for(j in 1:numSims){
      # designs at sim b
      bh = boxhills[[k]][[j]]
      # q = qs[[k]][[j]]
      # q1 = q1s[[k]][[j]]
      qc = qcaps[[k]][[j]]
      b = buffers[[k]][[j]]
      qc2 = qcap_persists[[k]][[j]]
      b2 = buffer_persists[[k]][[j]]
      r = randoms[[k]][[j]]
      sf = spacefills[[k]][[j]]
      # sequence of PPHs for each design
      RSS.bh = getRSS01(bh, model0, model1, x_seq, y_seq_mat[, j])
      # RSS.q = getRSS01(q, model0, model1, x_seq, y_seq_mat[, j])
      # RSS.q1 = getRSS01(q1, model0, model1, x_seq, y_seq_mat[, j])
      RSS.qc = getRSS01(qc, model0, model1, x_seq, y_seq_mat[, j])
      RSS.b = getRSS01(b, model0, model1, x_seq, y_seq_mat[, j])
      RSS.qc2 = getRSS01(qc2, model0, model1, x_seq, y_seq_mat[, j])
      RSS.b2 = getRSS01(b2, model0, model1, x_seq, y_seq_mat[, j])
      RSS.r = getRSS01(r, model0, model1, x_seq, y_seq_mat[, j])
      RSS.sf = getRSS01(sf, model0, model1, x_seq, y_seq_mat[, j])
      # master data frame
      RSS.bh$type = "boxhill"
      # RSS.q$type = "q"
      # RSS.q1$type = "q=1"
      RSS.qc$type = "qcap"
      RSS.b$type = "augdist"
      RSS.qc2$type = "qcap2"
      RSS.b2$type = "augdist2"
      RSS.r$type = "random"
      RSS.sf$type = "spacefill"
      # RSS.tmp = rbind(RSS.bh, RSS.q, RSS.q1, RSS.qc, RSS.b, RSS.r, RSS.sf)
      RSS.tmp = rbind(RSS.bh, RSS.qc, RSS.b, RSS.qc2, RSS.b2, RSS.r, RSS.sf)
      RSS.tmp$sim = j
      RSS.tmp$input = k
      RSS.df = rbind(RSS.df, RSS.tmp)
    }
  }
  
  RSS0mean = aggregate(RSS.df$RSS0, by = list(RSS.df$type, RSS.df$input), 
                       FUN = function(x) mean(x, na.rm = TRUE))
  names(RSS0mean) = c("type", "input", "value")
  RSS0mean$RSS = "RSS0"
  RSS1mean = aggregate(RSS.df$RSS1, by = list(RSS.df$type, RSS.df$input), 
                       FUN = function(x) mean(x, na.rm = TRUE))
  names(RSS1mean) = c("type", "input", "value")
  RSS1mean$RSS = "RSS1"
  RSS01mean = aggregate(RSS.df$RSS01, by = list(RSS.df$type, RSS.df$input), 
                        FUN = function(x) mean(x, na.rm = TRUE))
  names(RSS01mean) = c("type", "input", "value")
  RSS01mean$RSS = "RSS01"
  
  RSSmean = rbind(RSS0mean, RSS1mean, RSS01mean)
  RSSmean$type = factor(RSSmean$type)
  RSSmean$RSS = factor(RSSmean$RSS)
  RSSmean$input = factor(RSSmean$input, 
                         labels = c("extrapolation", "inc spread", "even coverage"))
  
  # RSS1
  RSS1.plt = ggplot(dplyr::filter(RSSmean, RSS == "RSS1"), 
                    aes(x = input, y = value, group = type, color = type, 
                        linetype = type)) + 
    geom_point() + 
    
    geom_path() +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    labs(y = "RSS1", x = "Initial Data") 
  plot(RSS1.plt)
  
  ggsave(
    filename = paste0("20210615_scen", scenario, "_rsst.pdf"), 
    plot = RSS1.plt, 
    width = 6, height = 4, units = c("in")
  )
  print(paste("scenario", scenario, 
              "################################################################"))
}
