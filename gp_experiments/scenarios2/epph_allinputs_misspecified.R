for(scenario in c(3.2, 4.2, 5.2, 6.2)){
  ################################################################################
  # last updated: 05/25/2021
  # purpose: to test seqmedgp for scenario 3:
  #   squared exponential vs. another squared exponential,
  #   where the true function is matern
  
  # scenario = 4 # scenarios: 3.2, 4.2, 5.2, 6.2
  scenario_subtypes = unlist(strsplit(as.character(scenario), split = "\\."))
  seq.type = 1 # 1 = fully sequential, 2 = stage-sequential 3x5
  
  ################################################################################
  # Sources/Libraries
  ################################################################################
  output_home = paste0("gp_experiments/scenarios2/scenarios2_misspecified/outputs")
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
  sigmasq_signals = c(1 - 1e-10, sigmasq_signal)
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
  if(scenario_subtypes[1] == 3){
    type01 = c("squaredexponential", "squaredexponential")
    typeT = "matern"
    l01= c(0.005, 0.01)
    lT = 0.01
  } else if(scenario_subtypes[1] == 4){
    type01 = c("matern", "squaredexponential")
    typeT = "periodic"
    l01= c(0.01, 0.01)
    lT = 0.01
  } else if(scenario_subtypes[1] == 5){
    type01 = c("matern", "periodic")
    typeT = "squaredexponential"
    l01= c(0.01, 0.01)
    lT = 0.01
  } else if(scenario_subtypes[1] == 6){
    type01 = c("squaredexponential", "periodic")
    typeT = "matern"
    l01= c(0.01, 0.01)
    lT = 0.01
  }
  
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
  # read in the data
  
  boxhills = list()
  qs = list()
  q1s = list()
  seqmeds = list()
  buffers = list()
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
    qs[[i]] = readRDS(paste0(
      output_home,
      "/scenario", scenario, "_seqmed", 
      "_q",
      "_seq", seq.type,
      filename_append.tmp))
    q1s[[i]] = readRDS(paste0(
      output_home,
      "/scenario", scenario, "_seqmed", 
      "_uniform",
      "_seq", seq.type,
      filename_append.tmp))
    seqmeds[[i]] = readRDS(paste0(
      output_home,
      "/scenario", scenario, "_seqmed",
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
  model0 = list(type = type01[1], l = l01[1], signal.var = sigmasq_signals[1], 
                measurement.var = sigmasq_measuremt)
  model1 = list(type = type01[2], l = l01[2], signal.var = sigmasq_signals[2], 
                measurement.var = sigmasq_measuremt)
  modelT = list(type = typeT, l = lT, signal.var = sigmasq_signal, 
                measurement.var = sigmasq_measuremt)
  
  # calculate the final posterior probability
  getPPH = function(design, model0, model1, modelT){
    index = length(as.vector(na.omit(design$y.new)))
    y.tmp = c(design$y, as.vector(na.omit(design$y.new)))
    x.tmp = c(design$x, as.vector(na.omit(design$x.new)))
    PPHs.tmp = getHypothesesPosteriors(
      prior.probs = rep(1 / 3, 3), 
      evidences = c(
        Evidence_gp(y.tmp, x.tmp, model0),
        Evidence_gp(y.tmp, x.tmp, model1),
        Evidence_gp(y.tmp, x.tmp, modelT)
      )
    )
    return(data.frame(
      index = index, 
      "H0" = PPHs.tmp[1], "H1" = PPHs.tmp[2], "HT" = PPHs.tmp[3]))
  }
  
  PPH = data.frame(
    PPH0 = numeric(), PPH1 = numeric(), PPHT = numeric(), 
    type = character(), sim = numeric(), input = numeric())
  for(k in 1:3){
    for(j in 1:numSims){
      # designs at sim b
      bh = boxhills[[k]][[j]]
      q = qs[[k]][[j]]
      q1 = q1s[[k]][[j]]
      sm = seqmeds[[k]][[j]]
      b = buffers[[k]][[j]]
      r = randoms[[k]][[j]]
      sf = spacefills[[k]][[j]]
      # sequence of PPHs for each design
      PPH.bh = getPPH(bh, model0, model1, modelT)
      PPH.q = getPPH(q, model0, model1, modelT)
      PPH.q1 = getPPH(q1, model0, model1, modelT)
      PPH.sm = getPPH(sm, model0, model1, modelT)
      PPH.b = getPPH(b, model0, model1, modelT)
      PPH.r = getPPH(r, model0, model1, modelT)
      PPH.sf = getPPH(sf, model0, model1, modelT)
      # master data frame
      PPH.bh$type = "boxhill"
      PPH.q$type = "q"
      PPH.q1$type = "seqmed,q1"
      PPH.sm$type = "seqmed"
      PPH.b$type = "augdist"
      PPH.r$type = "random"
      PPH.sf$type = "spacefill"
      PPH.tmp = rbind(
        PPH.bh, PPH.q, PPH.q1, PPH.sm, PPH.b, PPH.r, PPH.sf)
      PPH.tmp$sim = j
      PPH.tmp$input = k
      PPH = rbind(PPH, PPH.tmp)
    }
  }
  
  ### average over all, regardless of index of last point
  
  PPH0mean = aggregate(PPH$H0, by = list(PPH$type, PPH$input), 
                       FUN = function(x) mean(x, na.rm = TRUE))
  names(PPH0mean) = c("type", "input", "value")
  PPH0mean$Hypothesis = "H0"
  PPH1mean = aggregate(PPH$H1, by = list(PPH$type, PPH$input), 
                       FUN = function(x) mean(x, na.rm = TRUE))
  names(PPH1mean) = c("type", "input", "value")
  PPH1mean$Hypothesis = "H1"
  PPHTmean = aggregate(PPH$HT, by = list(PPH$type, PPH$input), 
                       FUN = function(x) mean(x, na.rm = TRUE))
  names(PPHTmean) = c("type", "input", "value")
  PPHTmean$Hypothesis = "HT"
  
  PPHmean = rbind(PPH0mean, PPH1mean, PPHTmean)
  PPHmean$type = factor(PPHmean$type)
  PPHmean$Hypothesis = factor(PPHmean$Hypothesis)
  PPHmean$input = factor(PPHmean$input, 
                         labels = c("extrapolation", "inc spread", "even coverage"))
  
  PPHT.plt = ggplot(dplyr::filter(PPHmean, Hypothesis == "HT"), 
                    aes(x = input, y = value, group = type, color = type, 
                        linetype = type, shape = type)) + 
    geom_path() +
    geom_point() +
    ylim(0, 1) + 
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    labs(y = "E[P(HT|X, Y)|X]", x = "Initial Data")
  PPHT.plt
  
  ggsave(
    filename = paste0("20210604_scen", scenario, "_eppht.pdf"), 
    plot = PPHT.plt, 
    width = 6, height = 4, units = c("in")
  )
  print(paste("scenario", scenario, 
              "################################################################"))
  
}