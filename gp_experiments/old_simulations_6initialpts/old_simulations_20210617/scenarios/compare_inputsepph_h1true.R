for(scenario in c(1, 2)){
  for(input.type in 1:3){
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
    output_home = paste0("gp_experiments/scenarios/scenarios_h1true/outputs")
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
    # read in the data
    
    boxhills = list()
    qs = list()
    q1s = list()
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
    
    # input set
    bh.in = boxhills[[input.type]]
    q.in = qs[[input.type]]
    q1.in = q1s[[input.type]]
    buf.in = buffers[[input.type]]
    ran.in = randoms[[input.type]]
    sf.in = spacefills[[input.type]]
    
    # input set
    if(input.type == 1){
      x_input = x_in1
      x_input_idx = x_in1_idx
    } else if(input.type == 2){
      x_input = x_in2
      x_input_idx = x_in2_idx
    } else if(input.type == 3){
      x_input = x_in3
      x_input_idx = x_in3_idx
    }
    
    getPPHinit = function(y, x, model0, model1){
      PPHs.tmp = getHypothesesPosteriors(
        prior.probs = prior_probs, 
        evidences = c(
          Evidence_gp(y, x, model0),
          Evidence_gp(y, x, model1)
        )
      )
      return(data.frame(H0 = PPHs.tmp[1], H1 = PPHs.tmp[2]))
    }
    
    PPHs = data.frame(
      H0 = numeric(), H1 = numeric())
    for(j in 1:numSims){
      y_seq = y_seq_mat[, j]
      y_input = y_seq[x_input_idx]
      PPH.tmp = getPPHinit(y_input, x_input, model0, model1)
      PPHs = rbind(PPHs, PPH.tmp)
    }
    PPHmeans = colMeans(PPHs)
    
    print(paste0("scenario", scenario, ", input ", input.type, ", EPPH = ", paste(as.character(round(PPHmeans, 3)), collapse = ", "), 
                "################################################################"))
    
  }
}