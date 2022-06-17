
scenario = 1
# scenario = 1 : SE vs Matern, f ~ Matern
# scenario = 2 : Matern vs Periodic, f ~ Periodic
text_size = 12

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
# simulation settings
################################################################################

numSims = 100
numSeq = 15
seqN = 1
Nttl = numSeq * seqN
xmin = -1
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
################################################################################
################################################################################
# plots!!!
################################################################################
################################################################################
################################################################################



################################################################################
# plot the posterior probabilities of the hypotheses
################################################################################

# useful functions #############################################################

getPPHseq = function(design, model0, model1, n, randomize.order = FALSE){
  n.new = n - 1
  x.new.idx = design$x.new.idx
  x.new = design$x.new
  y.new = design$y.new
  len.tmp = length(as.vector(na.omit(y.new)))
  if(randomize.order){
    new.order = sample(1:len.tmp, len.tmp, replace = FALSE)
    x.new.idx = x.new.idx[new.order]
    x.new = x.new[new.order]
    y.new = y.new[new.order]
  }
  # calculate posterior probs for each new point
  PPH0_seq = rep(NA, len.tmp)
  PPH1_seq = rep(NA, len.tmp)
  for(i in 1:len.tmp){
    y.tmp = c(design$y.in, y.new[1:i])
    x.tmp = c(design$x.in, x.new[1:i])
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
  if(length(PPH0_seq) < n.new){
    PPH0_seq[(length(PPH0_seq) + 1):n.new] = PPH0_seq[length(PPH0_seq)]
  }
  if(length(PPH1_seq) < n.new){
    PPH1_seq[(length(PPH1_seq) + 1):n.new] = PPH1_seq[length(PPH1_seq)]
  }
  # include posterior probs for initial point(s)
  len.init = length(design$y.in)
  PPHs.init = getHypothesesPosteriors(
    prior.probs = prior_probs, 
    evidences = c(
      Evidence_gp(design$y.in, design$x.in, model0),
      Evidence_gp(design$y.in, design$x.in, model1)
    )
  )
  PPH0_seq = c(PPHs.init[1], PPH0_seq)
  PPH1_seq = c(PPHs.init[2], PPH1_seq)
  return(data.frame(
    index = 1:n, 
    PPH0 = PPH0_seq, 
    PPH1 = PPH1_seq
  ))
}

# design model probs ###########################################################

PPH_seq = data.frame(
  PPH0 = numeric(), PPH1 = numeric(), PPHT = numeric(), 
  Design = character(), sim = numeric())
for(j in 1:numSims){
  # designs at sim b
  bh = boxhill_sims[[j]]
  lo = leaveout_sims[[j]]
  r = random_sims[[j]]
  g = grid_sims[[j]]
  # sequence of PPHs for each design
  PPH_seq.bh = getPPHseq(bh, model0, model1, Nttl)
  PPH_seq.lo = getPPHseq(lo, model0, model1, Nttl)
  PPH_seq.r = getPPHseq(r, model0, model1, Nttl)
  PPH_seq.g = getPPHseq(g, model0, model1, Nttl, randomize.order = TRUE)
  # master data frame
  PPH_seq.bh$Design = "(B) BoxHill"
  PPH_seq.lo$Design = "(A) SeqMED" 
  PPH_seq.r$Design = "(D) Random"
  PPH_seq.g$Design = "(C) Grid"
  PPH_seq.tmp = rbind(
    PPH_seq.bh, PPH_seq.lo,
    PPH_seq.r, PPH_seq.g)
  PPH_seq.tmp$sim = j
  PPH_seq = rbind(PPH_seq, PPH_seq.tmp)
}

PPH0mean_seq = aggregate(PPH_seq$PPH0, by = list(PPH_seq$index, PPH_seq$Design), 
                         FUN = function(x) mean(x, na.rm = TRUE))
names(PPH0mean_seq) = c("index", "Design", "value")
PPH0mean_seq$Hypothesis = "H0"
PPH1mean_seq = aggregate(PPH_seq$PPH1, by = list(PPH_seq$index, PPH_seq$Design), 
                         FUN = function(x) mean(x, na.rm = TRUE))
names(PPH1mean_seq) = c("index", "Design", "value")
PPH1mean_seq$Hypothesis = "H1"
PPH_seq = rbind(PPH0mean_seq, PPH1mean_seq) %>% 
  filter(Hypothesis != "H0")
epph.plt = ggplot(
  PPH_seq, aes(
    x = index, y = value, color = Design, linetype = Design, shape = Design)) + 
  facet_wrap(~Hypothesis) + 
  geom_path() + 
  geom_point() +
  ylim(0, 1) +
  labs(x = element_blank(), y = element_blank()) + 
  scale_x_continuous(breaks = c(5, 10, 15)) +
  theme_bw() + 
  theme(text = element_text(size = text_size), 
        panel.background = element_rect(fill = "transparent"), 
        plot.background = element_rect(
          fill = "transparent", color = NA), 
        legend.background = element_rect(fill = "transparent"), 
        legend.box.background = element_rect(
          fill = "transparent", color = NA))

epph.plt
ggsave(
  filename = paste0(scenario_name, "_epph.pdf"),
  plot = epph.plt,
  width = 4, height = 2.5, units = c("in")
)

################################################################################
# plot the designs
################################################################################
# all of the designs
idx = 1
designs = list(
  "(B) BoxHill" = boxhill_sims[[idx]], 
  "(C) Grid" = grid_sims[[idx]], 
  "(D) Random" = random_sims[[idx]],
  "(A) SeqMED" = leaveout_sims[[idx]])
design.names = sort(names(designs), decreasing = TRUE)

x.mat = matrix(NA, nrow = Nttl, ncol = length(designs))
for(i in 1:length(designs)){
  x.mat[, i] = c(designs[[i]]$x.in, designs[[i]]$x.new)
}
colnames(x.mat) = names(designs)
data.gg = as.data.frame(x.mat)
data.gg$index = as.character(1:Nttl)
data.gg = reshape2::melt(data.gg, id.vars = "index", variable.name = "Design")
data.gg$Design = factor(data.gg$Design, levels = design.names)

text.gg = dplyr::filter(data.gg, index %in% as.character(1:Nttl))
des.plt = ggplot() + 
  geom_point(data = data.gg, 
             mapping = aes(x = value, y = Design, color = Design), 
             inherit.aes = FALSE, alpha = 0.25) + 
  xlim(c(xmin, xmax)) + 
  labs(x = element_blank(), y = element_blank()) +
  theme_bw() + 
  theme(legend.position = "none", text = element_text(size = text_size), 
        panel.background = element_rect(fill = "transparent"), 
        plot.background = element_rect(
          fill = "transparent", color = NA), 
        legend.background = element_rect(fill = "transparent"), 
        legend.box.background = element_rect(
          fill = "transparent", color = NA))
if(typeT == "periodic"){
  des.plt = des.plt + geom_vline(
    xintercept = pT * (0:floor((xmax - xmin) / pT)), #color = "gray", 
    alpha = 0.125, linetype = 2)
}

des.plt
ggsave(
  filename = paste0(scenario_name, "_design.pdf"),
  plot = des.plt,
  width = 5.5, height = 2, units = c("in")
)
