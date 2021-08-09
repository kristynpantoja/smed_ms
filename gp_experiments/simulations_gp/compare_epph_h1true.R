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
    sims_dir = "gp_experiments/simulations_gp"
    output_dir = paste0(
      sims_dir, "/modelselection_designs/scenarios_misspecified/outputs")
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
    nugget = sigmasq_measuremt
    prior_probs = rep(1 / 2, 2)
    
  
    ################################################################################
    # Scenario settings
    ################################################################################
    l01= c(0.01, 0.01)
    if(scenario == 1){
      type01 = c("squaredexponential", "matern")
      model0 = list(type = type01[1], l = l01[1], signal.var = sigmasq_signal, 
                    measurement.var = nugget)
      model1 = list(type = type01[2], l = l01[2], signal.var = sigmasq_signal, 
                    measurement.var = nugget)
    } else if(scenario == 2){
      pT = 0.26
      type01 = c("matern", "periodic")
      model0 = list(type = type01[1], l = l01[1], signal.var = sigmasq_signal, 
                    measurement.var = nugget)
      model1 = list(type = type01[2], l = l01[2], signal.var = sigmasq_signal, 
                    measurement.var = nugget, p = pT)
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
    
    getPPHseq = function(design, model0, model1){
      PPH0_seq = rep(NA, length(as.vector(na.omit(design$y.new))))
      PPH1_seq = rep(NA, length(as.vector(na.omit(design$y.new))))
      for(i in 1:length(as.vector(na.omit(design$y.new)))){
        y.tmp = c(design$y, as.vector(na.omit(design$y.new))[1:i])
        x.tmp = c(design$x, as.vector(na.omit(design$x.new))[1:i])
        PPHs.tmp = getHypothesesPosteriors(
          prior.probs = prior_probs, 
          evidences = c(
            Evidence_gp(y.tmp, x.tmp, model0),
            Evidence_gp(y.tmp, x.tmp, model1)
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
      # designs at sim b
      bh = boxhill_sims[[j]]
      qc = qcap_sims[[j]]
      lo = leaveout_sims[[j]]
      qc2 = qcap_persist_sims[[j]]
      kq = persist_sims[[j]]
      r = random_sims[[j]]
      g = grid_sims[[j]]
      # sequence of PPHs for each design
      PPH_seq.bh = getPPHseq(bh, model0, model1)
      PPH_seq.qc = getPPHseq(qc, model0, model1)
      PPH_seq.lo = getPPHseq(lo, model0, model1)
      PPH_seq.qc2 = getPPHseq(qc2, model0, model1)
      PPH_seq.kq = getPPHseq(kq, model0, model1)
      PPH_seq.r = getPPHseq(r, model0, model1)
      PPH_seq.g = getPPHseq(g, model0, model1)
      # master data frame
      PPH_seq.bh$type = "boxhill"
      PPH_seq.qc$type = "qcap"
      PPH_seq.lo$type = "leaveout"
      PPH_seq.qc2$type = "keepq2"
      PPH_seq.kq$type = "keepq"
      PPH_seq.r$type = "random"
      PPH_seq.g$type = "grid"
      PPH_seq.tmp = rbind(
        PPH_seq.bh, PPH_seq.qc, PPH_seq.lo, PPH_seq.qc2, PPH_seq.kq, 
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
    
    ggsave(
      filename = paste0("20210708_scen", scenario, "_epph.pdf"), 
      plot = epph.plt, 
      width = 6, height = 4, units = c("in")
    )
    print(paste("scenario", scenario, 
                "################################################################"))
    
}