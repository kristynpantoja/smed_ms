for(scenario in c(1.1, 2.1)){
  rm(list = ls())
  ################################################################################
  # last updated: 05/25/2021
  # purpose: to test seqmedgp for scenario 1:
  #   squared exponential vs. matern,
  #   where the true function is matern
  
  # scenario = 1.1 # scenarios: 1.1, 2.1
  seq.type = 1 # 1 = fully sequential, 2 = stage-sequential 3x5
  
  ################################################################################
  # Sources/Libraries
  ################################################################################
  scenario_subtypes = unlist(strsplit(as.character(scenario), split = "\\."))
  output_home = paste0( # works for scenarios1 OR scenarios2 simulations
    "gp_experiments/scenarios", 
    scenario_subtypes[2],
    "/scenario", scenario, "/outputs")
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
  
  # shared settings
  sigmasq = 1
  if(scenario_subtypes[1] == 1){
    nuggets = c(1e-15, sigmasq_measuremt)
  } else if(scenario_subtypes[1] == 2){
    nuggets = c(1e-5, sigmasq_measuremt)
  } else{
    stop("invalid scenario")
  }
  prior_probs = rep(1 / 2, 2)
  
  ################################################################################
  # input data
  ################################################################################
  
  # 1. make space-filling design
  # space_filling = seq(from = xmin, to = xmax, length.out = Nttl)
  space_filling_idx = c(1, 1 + ((numx - 1)/(Nttl - 1)) * 1:((numx - 1) / ((numx - 1)/(Nttl - 1))))
  space_filling = x_seq[space_filling_idx]
  
  # input set 1 (extrapolation)
  x_in1_idx = space_filling_idx[1:Nin]
  x_in1 = x_seq[x_in1_idx]
  x_spacefill1_idx = space_filling_idx[-c(1:Nin)]
  x_spacefill1 = x_seq[x_spacefill1_idx]
  # all.equal(space_filling, c(x_in1, x_spacefill1))
  
  # input set 2 (increasing spread)
  x_in2_idx = space_filling_idx[c(1, 2, 4, 7, 12, 21)]
  x_in2 = x_seq[x_in2_idx]
  x_spacefill2_idx = space_filling_idx[-c(1, 2, 4, 7, 12, 21)]
  x_spacefill2 = x_seq[x_spacefill2_idx]
  # all.equal(space_filling, sort(c(x_in2, x_spacefill2)))
  
  # input set 3 (space-filling / even coverage)
  x_in3_idx = c(1, 1 + ((numx - 1)/(Nin - 1)) * 1:((numx - 1) / ((numx - 1)/(Nin - 1))))
  x_in3 = x_seq[x_in3_idx]
  x_spacefill3_idx = space_filling_idx[!(space_filling_idx %in% x_in3_idx)]
  x_spacefill3 = x_seq[x_spacefill3_idx]
  # all.equal(space_filling, sort(c(x_in3, x_spacefill3)))
  
  # input set 4 (uniform / random)
  
  ################################################################################
  # Scenario settings
  ################################################################################
  if(scenario_subtypes[1] == 1){
    type01 = c("squaredexponential", "matern")
  } else if(scenario_subtypes[1] == 2){
    type01 = c("matern", "periodic")
  }
  typeT = type01[2]
  l01= c(0.01, 0.01)
  lT = l01[2]
  
  ################################################################################
  # import matern functions
  filename_append = ""
  if(!is.null(sigmasq_measuremt)){
    filename_append = paste0(
      "_noise", strsplit(as.character(sigmasq_measuremt), "-")[[1]][2])
  }
  simulated.functions = readRDS(paste0(
    output_home,
    "/scenario", scenario, "_simulated_functions", filename_append,
    "_seed", rng.seed,
    ".rds"))
  numSims = simulated.functions$numSims
  x_seq = simulated.functions$x
  numx = length(x_seq)
  null_cov = simulated.functions$null_cov
  null_mean = simulated.functions$null_mean
  y_seq_mat = simulated.functions$function_values_mat
  
  ################################################################################
  # read in the data
  
  boxhills = list()
  qs = list()
  buffers = list()
  randoms = list()
  spacefills = list()
  seqmed.ms = list()
  
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
    qs[[i]] = readRDS(paste0(
      output_home,
      "/scenario", scenario, "_seqmed", 
      "_q",
      "_seq", seq.type,
      filename_append.tmp))
    
    randoms[[i]] = readRDS(paste0(
      output_home, 
      "/scenario", scenario, "_random", 
      filename_append.tmp))
    spacefills[[i]] = readRDS(paste0(
      output_home, 
      "/scenario", scenario, "_spacefilling", 
      filename_append.tmp))
    
    seqmed.ms[[i]] = readRDS(paste0(
      output_home,
      "/scenario", scenario, "_seqmed", 
      "_error", 
      "_seq", seq.type,
      filename_append.tmp))
  }
  
  ################################################################################
  # make plots
  ################################################################################
  
  # models
  model0 = list(type = type01[1], l = l01[1], signal.var = sigmasq,
                measurement.var = nuggets[1])
  model1 = list(type = type01[2], l = l01[2], signal.var = sigmasq, 
                measurement.var = nuggets[2])
  
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
                                sigmasq, model0$measurement.var)
    pred1.tmp = getGPPredictive(x.other, x.tmp, y.tmp, 
                                model1$type, model1$l, 
                                sigmasq, model1$measurement.var)
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
      q = qs[[k]][[j]]
      b = buffers[[k]][[j]]
      r = randoms[[k]][[j]]
      sf = spacefills[[k]][[j]]
      m = seqmed.ms[[k]][[j]]
      # sequence of PPHs for each design
      RSS.bh = getRSS01(bh, model0, model1, x_seq, y_seq_mat[, j])
      RSS.q = getRSS01(q, model0, model1, x_seq, y_seq_mat[, j])
      RSS.b = getRSS01(b, model0, model1, x_seq, y_seq_mat[, j])
      RSS.r = getRSS01(r, model0, model1, x_seq, y_seq_mat[, j])
      RSS.sf = getRSS01(sf, model0, model1, x_seq, y_seq_mat[, j])
      RSS.m = getRSS01(m, model0, model1, x_seq, y_seq_mat[, j])
      # master data frame
      RSS.bh$type = "boxhill"
      RSS.q$type = "q"
      RSS.b$type = "augdist"
      RSS.r$type = "random"
      RSS.sf$type = "spacefill"
      RSS.m$type = "measmt"
      RSS.tmp = rbind(
        RSS.bh, RSS.q, RSS.b, RSS.r, RSS.sf, RSS.m)
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
  RSS1.plt
  
  
  ggsave(
    filename = paste0("20210525_scen", scenario, "_rsst.pdf"), 
    plot = RSS1.plt, 
    width = 6, height = 4, units = c("in")
  )
  print(paste("scenario", scenario, 
              "################################################################"))
}
