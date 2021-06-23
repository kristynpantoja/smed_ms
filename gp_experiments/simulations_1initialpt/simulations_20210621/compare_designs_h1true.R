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
    # if(scenario == 1){
    #   nuggets = c(1e-15, sigmasq_measuremt)
    # } else if(scenario == 2){
    #   nuggets = c(1e-5, sigmasq_measuremt)
    # } else{
    #   stop("invalid scenario")
    # }
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
    
      # filename_append.tmp for all methods alike
      filename_append.tmp = paste0(
        filename_append, 
        "_seed", rng.seed,
        ".rds"
      )
      boxhill_sims = readRDS(paste0(
        output_home,
        "/scenario", scenario, "_boxhill", 
        filename_append.tmp))
      leaveout_sims = readRDS(paste0(
        output_home,
        "/scenario", scenario, "_seqmed", 
        "_leaveout", 
        "_seq", seq.type,
        filename_append.tmp))
      qcap_sims = readRDS(paste0(
        output_home,
        "/scenario", scenario, "_seqmed", 
        "_cap",
        "_seq", seq.type,
        filename_append.tmp))
      leaveout_persist_sims = readRDS(paste0(
        output_home,
        "/scenario", scenario, "_seqmed", 
        "_leaveout_persist", 
        "_seq", seq.type,
        filename_append.tmp))
      qcap_persist_sims = readRDS(paste0(
        output_home,
        "/scenario", scenario, "_seqmed", 
        "_cap_persist",
        "_seq", seq.type,
        filename_append.tmp))
      
      random_sims = readRDS(paste0(
        sims_dir, 
        "/spacefilling_designs/outputs/random", 
        "_", typeT,
        "_l", lT,
        filename_append.tmp))
      spacefill_sims = readRDS(paste0(
        sims_dir,
        "/spacefilling_designs/outputs/grid", 
        "_", typeT,
        "_l", lT,
        filename_append.tmp))
    
    ################################################################################
    # make plots
    ################################################################################
    
    # input set
    bh.in = boxhills[[input.type]]
    # q.in = qs[[input.type]]
    # q1.in = q1s[[input.type]]
    qcap.in = qcaps[[input.type]]
    buf.in = buffers[[input.type]]
    qcap2.in = qcap_persists[[input.type]]
    buf2.in = buffer_persists[[input.type]]
    ran.in = randoms[[input.type]]
    sf.in = spacefills[[input.type]]
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
    
    # all 6 designs
    idx = 1
    # designs = list(bh.in[[idx]], q.in[[idx]], q1.in[[idx]], qcap.in[[idx]], buf.in[[idx]])
    # design.names = c("boxhill", "q", "q=1", "qcap", "augdist")
    # design.levels = c("augdist", "qcap", "q=1", "q", "boxhill")
    designs = list(bh.in[[idx]], qcap.in[[idx]], buf.in[[idx]], qcap2.in[[idx]], buf2.in[[idx]])
    design.names = c("boxhill", "qcap", "augdist", "qcap2", "augdist2")
    design.levels = c("augdist", "augdist2", "qcap", "qcap2", "boxhill")
    
    x.new.mat = matrix(NA, nrow = Nnew, ncol = length(designs))
    for(i in 1:length(designs)){
      x.new.mat[, i] = designs[[i]]$x.new
    }
    
    data.gg = data.frame(
      index = as.character(rep(1:Nnew, length(designs))), 
      type = factor(rep(design.names, each = Nnew), levels = design.levels), 
      value = as.vector(x.new.mat)
    )
    data.gg0 = data.frame(
      type = factor(rep(design.names, each = Nin), levels = design.levels), 
      input = rep(x_input, length(designs))
    )
    text.gg = dplyr::filter(data.gg, index %in% as.character(1:Nnew))
    des.plt = ggplot() + 
      geom_point(data = data.gg0, 
                 mapping = aes(x = input, y = type)) +
      geom_point(data = data.gg, 
                 mapping = aes(x = value, y = type, color = type), 
                 inherit.aes = FALSE) + 
      geom_text(data = text.gg, 
                aes(x = value, y = type, label = index), 
                vjust = -0.65 * as.numeric(paste(text.gg$index)), size = 2) +
      xlim(c(xmin, xmax))
    plot(des.plt)
    
    ggsave(
      filename = paste0("20210615_scen", scenario, "_in", input.type, "_design.pdf"), 
      plot = des.plt, 
      width = 6, height = 4, units = c("in")
    )
    print(paste("scenario", scenario, ", input ", input.type, 
                "################################################################"))
    
}