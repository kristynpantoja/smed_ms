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
  } else{
    stop("invalid scenario")
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
  
  # calculate the final posterior probability
  getPPH = function(design, model0, model1){
    y.tmp = c(design$y, as.vector(na.omit(design$y.new)))
    x.tmp = c(design$x, as.vector(na.omit(design$x.new)))
    PPHs.tmp = getHypothesesPosteriors(
      prior.probs = prior_probs, 
      evidences = c(
        Evidence_gp(y.tmp, x.tmp, model0),
        Evidence_gp(y.tmp, x.tmp, model1)
      )
    )
    return(data.frame("H0" = PPHs.tmp[1], "H1" = PPHs.tmp[2]))
  }
  
  PPH = data.frame(
    PPH0 = numeric(), PPH1 = numeric(), 
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
      PPH.bh = getPPH(bh, model0, model1)
      # PPH.q = getPPH(q, model0, model1)
      # PPH.q1 = getPPH(q1, model0, model1)
      PPH.qc = getPPH(qc, model0, model1)
      PPH.b = getPPH(b, model0, model1)
      PPH.qc2 = getPPH(qc2, model0, model1)
      PPH.b2 = getPPH(b2, model0, model1)
      PPH.r = getPPH(r, model0, model1)
      PPH.sf = getPPH(sf, model0, model1)
      # master data frame
      PPH.bh$type = "boxhill"
      # PPH.q$type = "q"
      # PPH.q1$type = "q=1"
      PPH.qc$type = "qcap"
      PPH.b$type = "augdist"
      PPH.qc2$type = "qcap2"
      PPH.b2$type = "augdist2"
      PPH.r$type = "random"
      PPH.sf$type = "spacefill"
      # PPH.tmp = rbind(PPH.bh, PPH.q, PPH.q1, PPH.qc, PPH.b, PPH.r, PPH.sf)
      PPH.tmp = rbind(PPH.bh, PPH.qc, PPH.b, PPH.qc2, PPH.b2, PPH.r, PPH.sf)
      PPH.tmp$sim = j
      PPH.tmp$input = k
      PPH = rbind(PPH, PPH.tmp)
    }
  }
  
  PPH0mean = aggregate(PPH$H0, by = list(PPH$type, PPH$input), 
                       FUN = function(x) mean(x, na.rm = TRUE))
  names(PPH0mean) = c("type", "input", "value")
  PPH0mean$Hypothesis = "H0"
  PPH1mean = aggregate(PPH$H1, by = list(PPH$type, PPH$input), 
                       FUN = function(x) mean(x, na.rm = TRUE))
  names(PPH1mean) = c("type", "input", "value")
  PPH1mean$Hypothesis = "H1"
  
  PPHmean = rbind(PPH0mean, PPH1mean)
  PPHmean$type = factor(PPHmean$type)
  PPHmean$Hypothesis = factor(PPHmean$Hypothesis)
  PPHmean$input = factor(
    PPHmean$input, labels = c("extrapolation", "inc spread", "even coverage"))
  
  PPH1.plt = ggplot(dplyr::filter(PPHmean, Hypothesis == "H1"), 
                    aes(x = input, y = value, group = type, color = type, 
                        linetype = type, shape = type)) + 
    geom_path() +
    geom_point() +
    ylim(0, 1) + 
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    labs(y = "P(H1|X, Y)", x = "Initial Data")
  plot(PPH1.plt)
  
  ggsave(
    filename = paste0("20210615_scen", scenario, "_eppht.pdf"), 
    plot = PPH1.plt, 
    width = 6, height = 4, units = c("in")
  )
  print(paste("scenario", scenario,
              "################################################################"))
  
}