
dimT = 1
text_size = 12

################################################################################
# Sources/Libraries
################################################################################
sims_dir = "gp_experiments/gpvs"
output_dir = paste0(sims_dir, "/modelselection_designs/outputs")
data_dir = paste0(sims_dir, "/simulated_data/outputs")
functions_dir = "functions"

# for seqmed design
source(paste(functions_dir, "/SeqMEDgpvs.R", sep = ""))
source(paste(functions_dir, "/SeqMEDgpvs_batch.R", sep = ""))
source(paste(functions_dir, "/charge_function_q.R", sep = ""))
source(paste(functions_dir, "/covariance_functions.R", sep = ""))
source(paste(functions_dir, "/wasserstein_distance.R", sep = ""))
source(paste(functions_dir, "/gp_predictive.R", sep = ""))

# for box-hill design
source(paste(functions_dir, "/boxhill.R", sep = ""))
source(paste(functions_dir, "/boxhill_gp.R", sep = ""))
source(paste(functions_dir, "/boxhill_gpvs.R", sep = ""))
source(paste(functions_dir, "/kl_divergence.R", sep = ""))

library(mvtnorm)
rng.seed = 123

# parallelization settings
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

seq.type = 1
lT = 0.01
typeT = "squaredexponential"

numSims = 100
Nin = 1
if(seq.type == 1){
  numSeq = 9
  seqN = 1
} else if(seq.type == 2){
  numSeq = 3
  seqN = 3
}
Nnew = numSeq * seqN
Nttl = Nin + Nnew 
xmin = 0
xmax = 1
p = 2
k = 4 * p
sigmasq_measuremt = 1e-10
sigmasq_signal = 1

# shared settings
prior_probs = rep(1 / 2, 2)


################################################################################
# Scenario settings
################################################################################

l01= c(lT, lT)
type01 = c(typeT, typeT)
indices0 = c(1)
indices1 = c(1, 2)

model0 = list(type = type01[1], l = l01[1], signal.var = sigmasq_signal, 
              indices = indices0,
              measurement.var = sigmasq_measuremt)
model1 = list(type = type01[2], l = l01[2], signal.var = sigmasq_signal, 
              indices = indices1, 
              measurement.var = sigmasq_measuremt)

################################################################################
# import data
filename_append = ""
if(!is.null(sigmasq_measuremt)){
  filename_append = "_noise"
}
simulated_data_file = paste0(
  data_dir,
  "/", typeT,
  "_l", lT,
  "_dim", dimT,
  filename_append, 
  "_seed", rng.seed,
  ".rds")
simulated.data = readRDS(simulated_data_file)
numSims = simulated.data$numSims
x_seq = simulated.data$x_seq
x_grid = as.matrix(simulated.data$x_grid)
numx = length(x_seq)
null_cov = simulated.data$null_cov
null_mean = simulated.data$null_mean
y_seq_mat = simulated.data$function_values_mat

################################################################################
# initial design
x_input_idx = sample(1:numx, Nin)
x_input = x_grid[x_input_idx + numx * (1:Nin), , drop = FALSE]

################################################################################
# read in the designs

# filename_append.tmp for all methods alike
filename_append.tmp = paste0(
  filename_append, 
  "_seed", rng.seed,
  ".rds"
)

leaveout_sims = readRDS(paste0(
  output_dir,
  "/seqmed", 
  "_leaveout", 
  "_seq", seq.type,
  "_l", lT,
  "_dim", dimT,
  filename_append.tmp))
cap_sims = readRDS(paste0(
  output_dir,
  "/seqmed", 
  "_cap", 
  "_seq", seq.type,
  "_l", lT,
  "_dim", dimT,
  filename_append.tmp))
keepq_cap_sims = readRDS(paste0(
  output_dir,
  "/seqmed", 
  "_keepq_cap", 
  "_seq", seq.type,
  "_l", lT,
  "_dim", dimT,
  filename_append.tmp))
keepq_sims = readRDS(paste0(
  output_dir,
  "/seqmed", 
  "_keepq", 
  "_seq", seq.type,
  "_l", lT,
  "_dim", dimT,
  filename_append.tmp))
boxhill_sims = readRDS(paste0(
  output_dir,
  "/boxhill", 
  "_seq", seq.type,
  "_l", lT,
  "_dim", dimT,
  filename_append.tmp))

random_sims = readRDS(paste0(
  sims_dir, 
  "/other_designs/outputs/random", 
  "_", typeT,
  "_l", lT,
  "_dim", dimT,
  filename_append.tmp))
grid_sims = readRDS(paste0(
  sims_dir,
  "/other_designs/outputs/grid", 
  "_", typeT,
  "_l", lT,
  "_dim", dimT,
  filename_append.tmp))
diagonal_sims = readRDS(paste0(
  sims_dir,
  "/other_designs/outputs/diagonal", 
  "_", typeT,
  "_l", lT,
  "_dim", dimT,
  filename_append.tmp))
x2_sims = readRDS(paste0(
  sims_dir,
  "/other_designs/outputs/x2constant", 
  "_", typeT,
  "_l", lT,
  "_dim", dimT,
  filename_append.tmp))



################################################################################
# compare it to space-filling designs

getPPHseq = function(design, model0, model1, n, randomize.order = FALSE){
  x.new.idx = design$x.new.idx
  x.new = design$x.new
  y.new = design$y.new
  len.tmp = length(as.vector(na.omit(y.new)))
  if(randomize.order){
    new.order = sample(1:len.tmp, len.tmp, replace = FALSE)
    x.new.idx = x.new.idx[new.order]
    x.new = x.new[new.order, ]
    y.new = y.new[new.order]
  }
  # calculate posterior probs for each new point
  PPH0_seq = rep(NA, len.tmp)
  PPH1_seq = rep(NA, len.tmp)
  for(i in 1:len.tmp){
    y.tmp = c(design$y.in, y.new[1:i])
    x.tmp = rbind(design$x.in, x.new[1:i, , drop = FALSE])
    PPHs.tmp = getHypothesesPosteriors(
      prior.probs = prior_probs, 
      evidences = c(
        Evidence_gpvs(y.tmp, x.tmp, model0),
        Evidence_gpvs(y.tmp, x.tmp, model1)
      )
    )
    PPH0_seq[i] = PPHs.tmp[1]
    PPH1_seq[i] = PPHs.tmp[2]
  }
  if(length(PPH0_seq) < n){
    PPH0_seq[(length(PPH0_seq) + 1):n] = PPH0_seq[length(PPH0_seq)]
  }
  if(length(PPH1_seq) < n){
    PPH1_seq[(length(PPH1_seq) + 1):n] = PPH1_seq[length(PPH1_seq)]
  }
  # include posterior probs for initial point(s)
  len.init = length(design$y.in)
  PPHs.init = getHypothesesPosteriors(
    prior.probs = prior_probs, 
    evidences = c(
      Evidence_gpvs(design$y.in, design$x.in, model0),
      Evidence_gpvs(design$y.in, design$x.in, model1)
    )
  )
  PPH0_seq = c(PPHs.init[1], PPH0_seq)
  PPH1_seq = c(PPHs.init[2], PPH1_seq)
  return(data.frame(
    index = 0:((len.init - 1) + n), 
    PPH0 = PPH0_seq, 
    PPH1 = PPH1_seq
  ))
}

PPH_seq = data.frame(
  PPH0 = numeric(), PPH1 = numeric(), PPHT = numeric(), 
  type = character(), sim = numeric())
for(j in 1:numSims){
  
  # designs at sim b
  qc = cap_sims[[j]]
  kq = keepq_sims[[j]]
  kqc = keepq_cap_sims[[j]]
  lo = leaveout_sims[[j]]
  bh = boxhill_sims[[j]]
  diag = diagonal_sims[[j]]
  x2 = x2_sims[[j]]
  r = random_sims[[j]]
  g = grid_sims[[j]]
  # sequence of PPHs for each design
  PPH_seq.qc = getPPHseq(qc, model0, model1, Nnew)
  PPH_seq.kq = getPPHseq(kq, model0, model1, Nnew)
  PPH_seq.kqc = getPPHseq(kqc, model0, model1, Nnew)
  PPH_seq.lo = getPPHseq(lo, model0, model1, Nnew)
  PPH_seq.bh = getPPHseq(bh, model0, model1, Nnew)
  PPH_seq.diag = getPPHseq(diag, model0, model1, Nnew, randomize.order = TRUE)
  PPH_seq.x2 = getPPHseq(x2, model0, model1, Nnew, randomize.order = TRUE)
  PPH_seq.r = getPPHseq(r, model0, model1, Nnew)
  PPH_seq.g = getPPHseq(g, model0, model1, Nnew, randomize.order = TRUE)
  # master data frame
  PPH_seq.qc$Design = "cap q"
  PPH_seq.kq$Design = "keepq"
  PPH_seq.kqc$Design = "keepq2"
  PPH_seq.lo$Design = "SeqMED" 
  PPH_seq.bh$Design = "BoxHill"
  PPH_seq.diag$Design = "Diagonal"
  PPH_seq.x2$Design = "x2=1"
  PPH_seq.r$Design = "Random"
  PPH_seq.g$Design = "Grid"
  PPH_seq.tmp = rbind(
    PPH_seq.lo, PPH_seq.bh, PPH_seq.diag,
    PPH_seq.x2, PPH_seq.r, PPH_seq.g)
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

PPHmean_seq = rbind(PPH0mean_seq, PPH1mean_seq)
epph.plt = ggplot(PPHmean_seq, aes(x = index, y = value, color = Design, 
                                   linetype = Design, shape = Design)) + 
  facet_wrap(~Hypothesis) + 
  geom_path() + 
  geom_point() +
  ylim(0, 1) + 
  scale_x_continuous(breaks = c(0, 3, 6, 9)) +
  labs(x = element_blank(), y = element_blank()) +
  theme_bw() + 
  theme(text = element_text(size = text_size), 
        panel.background = element_rect(fill = "transparent"), 
        plot.background = element_rect(
          fill = "transparent", color = NA), 
        legend.background = element_rect(fill = "transparent"), 
        legend.box.background = element_rect(
          fill = "transparent", color = NA))
plot(epph.plt)

ggsave(
  filename = paste0("gpvs_dim", dimT, "_epph.pdf"),
  plot = epph.plt,
  width = 6, height = 2.5, units = c("in")
)
