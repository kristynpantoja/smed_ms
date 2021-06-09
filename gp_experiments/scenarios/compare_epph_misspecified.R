for(scenario in c(3, 4, 5, 6)){
  for(input.type in 1:3){
    ################################################################################
    # last updated: 05/25/2021
    # purpose: to test seqmedgp for scenario 3:
    #   squared exponential vs. another squared exponential,
    #   where the true function is matern
    
    # scenario = 4 # scenarios: 3, 4, 5, 6
    # input.type = 2 # 1 = extrapolation, 2 = inc spread, 3 = even coverage
    seq.type = 1 # 1 = fully sequential, 2 = stage-sequential 3x5
    
    ################################################################################
    # Sources/Libraries
    ################################################################################
    output_home = paste0("gp_experiments/scenarios/scenarios_misspecified/outputs")
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
    if(scenario == 3){
      type01 = c("squaredexponential", "squaredexponential")
      typeT = "matern"
      l01= c(0.005, 0.01)
      lT = 0.01
    } else if(scenario == 4){
      type01 = c("matern", "squaredexponential")
      typeT = "periodic"
      l01= c(0.01, 0.01)
      lT = 0.01
    } else if(scenario == 5){
      type01 = c("matern", "periodic")
      typeT = "squaredexponential"
      l01= c(0.01, 0.01)
      lT = 0.01
    } else if(scenario == 6){
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
    buffers = list()
    q1s = list()
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
    # make sequential EPPH plots
    ################################################################################
    
    # models
    model0 = list(type = type01[1], l = l01[1], signal.var = sigmasq_signal, 
                  measurement.var = nugget)
    model1 = list(type = type01[2], l = l01[2], signal.var = sigmasq_signal, 
                  measurement.var = nugget)
    modelT = list(type = typeT, l = lT, signal.var = sigmasq_signal, 
                  measurement.var = sigmasq_measuremt)
    
    # input set
    bh.in = boxhills[[input.type]]
    q.in = qs[[input.type]]
    q1.in = q1s[[input.type]]
    buf.in = buffers[[input.type]]
    ran.in = randoms[[input.type]]
    sf.in = spacefills[[input.type]]
    
    getPPHseq = function(design, model0, model1, modelT){
      PPH0_seq = rep(NA, length(as.vector(na.omit(design$y.new))))
      PPH1_seq = rep(NA, length(as.vector(na.omit(design$y.new))))
      PPHT_seq = rep(NA, length(as.vector(na.omit(design$y.new))))
      for(i in 1:length(as.vector(na.omit(design$y.new)))){
        y.tmp = c(design$y, as.vector(na.omit(design$y.new))[1:i])
        x.tmp = c(design$x, as.vector(na.omit(design$x.new))[1:i])
        PPHs.tmp = getHypothesesPosteriors(
          prior.probs = rep(1 / 3, 3), 
          evidences = c(
            Evidence_gp(y.tmp, x.tmp, model0),
            Evidence_gp(y.tmp, x.tmp, model1), 
            Evidence_gp(y.tmp, x.tmp, modelT)
          )
        )
        PPH0_seq[i] = PPHs.tmp[1]
        PPH1_seq[i] = PPHs.tmp[2]
        PPHT_seq[i] = PPHs.tmp[3]
      }
      if(length(PPH0_seq) < Nnew){
        PPH0_seq[(length(PPH0_seq) + 1):Nnew] = NA
        PPH1_seq[(length(PPH1_seq) + 1):Nnew] = NA
        PPHT_seq[(length(PPHT_seq) + 1):Nnew] = NA
      }
      return(data.frame(
        index = 1:Nnew, 
        "H0" = PPH0_seq, 
        "H1" = PPH1_seq, 
        "HT" = PPHT_seq
      ))
    }
    
    PPH_seq = data.frame(
      PPH0 = numeric(), PPH1 = numeric(), PPHT = numeric(), 
      type = character(), sim = numeric())
    for(j in 1:numSims){
      # designs at sim b
      bh = bh.in[[j]]
      q = q.in[[j]]
      q1 = q1.in[[j]]
      b = buf.in[[j]]
      r = ran.in[[j]]
      sf = sf.in[[j]]
      # sequence of PPHs for each design
      PPH_seq.bh = getPPHseq(bh, model0, model1, modelT)
      PPH_seq.q = getPPHseq(q, model0, model1, modelT)
      PPH_seq.q1 = getPPHseq(q1, model0, model1, modelT)
      PPH_seq.b = getPPHseq(b, model0, model1, modelT)
      PPH_seq.r = getPPHseq(r, model0, model1, modelT)
      PPH_seq.sf = getPPHseq(sf, model0, model1, modelT)
      # master data frame
      PPH_seq.bh$type = "boxhill"
      PPH_seq.q$type = "q"
      PPH_seq.q1$type = "seqmed,q1"
      PPH_seq.b$type = "augdist"
      PPH_seq.r$type = "random"
      PPH_seq.sf$type = "spacefill"
      PPH_seq.tmp = rbind(
        PPH_seq.bh, PPH_seq.q, PPH_seq.q1, PPH_seq.b, PPH_seq.r, PPH_seq.sf)
      PPH_seq.tmp$sim = j
      PPH_seq = rbind(PPH_seq, PPH_seq.tmp)
    }
    
    PPH0mean_seq = aggregate(PPH_seq$H0, by = list(PPH_seq$index, PPH_seq$type), 
                             FUN = function(x) mean(x, na.rm = TRUE))
    names(PPH0mean_seq) = c("index", "type", "value")
    PPH0mean_seq$Hypothesis = "H0"
    PPH1mean_seq = aggregate(PPH_seq$H1, by = list(PPH_seq$index, PPH_seq$type), 
                             FUN = function(x) mean(x, na.rm = TRUE))
    names(PPH1mean_seq) = c("index", "type", "value")
    PPH1mean_seq$Hypothesis = "H1"
    PPHTmean_seq = aggregate(PPH_seq$HT, by = list(PPH_seq$index, PPH_seq$type), 
                             FUN = function(x) mean(x, na.rm = TRUE))
    names(PPHTmean_seq) = c("index", "type", "value")
    PPHTmean_seq$Hypothesis = "HT"
    
    PPHmean_seq = rbind(PPH0mean_seq, PPH1mean_seq, PPHTmean_seq)
    epph.plt = ggplot(PPHmean_seq, aes(x = index, y = value, color = type, 
                                       linetype = type, shape = type)) + 
      facet_wrap(~Hypothesis) + 
      geom_path() + 
      geom_point() +
      theme_bw() +
      ylim(0, 1)
    epph.plt
    
    ggsave(
      filename = paste0("20210530_scen", scenario, "_in", input.type, "_epph.pdf"), 
      plot = epph.plt, 
      width = 6, height = 4, units = c("in")
    )
    
    print(paste("scenario", scenario, ", input ", input.type, 
                "################################################################"))
  }
}