# for(scenario in c(1, 2)){
scenario = 1
################################################################################
# last updated: 2/27/2022
# purpose: to test seqmedgp for scenario 1:
#   squared exponential vs. matern,
#   where the true function is matern

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
numSeq = 15
seqN = 1
Nttl = numSeq * seqN
xmin = -1
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
  # geom_text(data = text.gg, 
  #           aes(x = value, y = Design, label = index), 
  #           vjust = -0.65 * as.numeric(paste(text.gg$index)), size = 2) +
  xlim(c(xmin, xmax)) + 
  theme_bw() + 
  labs(x = element_blank(), y = element_blank()) +
  theme(legend.position = "none")
if(typeT == "periodic"){
  des.plt = des.plt + geom_vline(
    xintercept = pT * (0:floor((xmax - xmin) / pT)), #color = "gray", 
    alpha = 0.125, linetype = 2)
}
plot(des.plt)

# # slide plot
# ggsave(
#   filename = paste0("scen", scenario, "_", scenario_name, "_design.pdf"),
#   plot = des.plt,
#   width = 6, height = 4, units = c("in")
# )

# manuscript plot
ggsave(
  filename = paste0(scenario_name, "_design.pdf"),
  plot = des.plt,
  width = 4, height = 1.75, units = c("in")
)

# }
