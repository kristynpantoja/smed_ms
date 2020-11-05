################################################################################
### --- Gaussian Process Regression ---------------------------------------- ###
################################################################################

# --- Working Directory --- #
home = "/home/kristyn/Documents/research/seqmed/smed_ms"

# --- Sources/Libraries --- #
functions_home = paste(home, "/functions", sep="")
source(paste(functions_home, "/add_MMEDgp.R", sep = ""))
source(paste(functions_home, "/wasserstein_distance.R", sep = ""))
source(paste(functions_home, "/charge_function_q.R", sep = ""))
source(paste(functions_home, "/covariance_functions.R", sep = ""))
source(paste(functions_home, "/gp_predictive.R", sep = ""))
source(paste(functions_home, "/SMMEDgp.R", sep = ""))
###
source(paste(functions_home, "/boxhill.R", sep = ""))
source(paste(functions_home, "/boxhill_gp.R", sep = ""))

library(expm)
library(matrixStats)
library(MASS)
library(mvtnorm)
library(fields)
library(knitr)

# for plots
library(ggplot2)
library(ggpubr)
library(reshape2)
library(data.table)
gg_color_hue = function(n) {
  hues = seq(15, 275, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
image_path = "/home/kristyn/Pictures"

### --- shared settings for both 1d f & 2d f cases --- ###

# --- simulations --- #

# settings: Squared Exponential vs. Matern
N = 6
l01= c(0.1, 0.1)
type01 = c(1, 4)
numCandidates = 1001
nugget = 1e-10
N.new = 5

# x_seq, grid over which to generate subsequent functions
xmin = 0; xmax = 1
numx = 1001
x_seq = seq(from = xmin, to = xmax, length.out = numx) # set input points

# generate matern function (H1 is true)
set.seed(12)
null_cov = getCov(x_seq, x_seq, type01[2], l01[2])
null_mean = rep(0, numx)
y_seq = as.vector(rmvnorm(n = 1, mean = null_mean, sigma = null_cov)) # the function values


################################################################################
# input data
################################################################################

# input points : randomly-selected points
set.seed(1997)
x_input_idx = sample(1:numx, N)
x_input = x_seq[x_input_idx]
y_input = y_seq[x_input_idx]

################################################################################
# box and hill method
model0 = list(type = type01[1], l = l01[1])
model1 = list(type = type01[2], l = l01[2])

# calculate prior probabilities using preliminary data (input data)
prior_probs = rep(1 / 2, 2)

BHres = getBHgp_m2(y_input, prior_probs, x_input, model0, model1, x_seq, y_seq, N.new, nugget)
x_new_idx = BHres$x.new.idx
x_new = BHres$x.new
y_new = y_seq[x_new_idx]

# plot
x_input.gg = x_input
y_input.gg = y_input
x_new.gg = x_new
y_new.gg = y_new
H0_predfn = getGPPredictive(x_seq, x_input.gg, y_input.gg, type01[1], l01[1],
                            nugget = NULL)
H1_predfn = getGPPredictive(x_seq, x_input.gg, y_input.gg, type01[2], l01[2],
                            nugget = NULL)
err0 = 2 * sqrt(diag(H0_predfn$pred_var))
err1 = 2 * sqrt(diag(H1_predfn$pred_var))
ggdata = data.table(
  x = x_seq, 
  `True Function` = y_seq, 
  `H0 Predictive Mean` = H0_predfn$pred_mean, 
  `H1 Predictive Mean` = H1_predfn$pred_mean,
  lower0 = H0_predfn$pred_mean - err0, 
  lower1 = H1_predfn$pred_mean - err1, 
  upper0 = H0_predfn$pred_mean + err0,
  upper1 = H1_predfn$pred_mean + err1
)
yrange = range(ggdata$lower0, ggdata$lower1, 
               ggdata$upper0, ggdata$upper1, 
               na.rm = TRUE)
yrange[1] = yrange[1] - 1
ggdata$zero1 = NA
ggdata.melted = melt(ggdata, id.vars = c("x"), 
                     measure.vars = c("True Function", "H0 Predictive Mean", "H1 Predictive Mean"))
ggdata.lower = melt(ggdata, id.vars = c("x"), 
                    measure.vars = c("zero1", "lower0", "lower1"))
ggdata.upper = melt(ggdata, id.vars = c("x"), 
                    measure.vars = c("zero1", "upper0", "upper1"))
ggdata.melted = cbind(ggdata.melted, 
                      lower = ggdata.lower$value, 
                      upper = ggdata.upper$value)
ggdata_pts = data.table(
  x = c(x_input.gg, x_new.gg), 
  y = c(y_input.gg, y_new.gg), 
  color = c(rep(gg_color_hue(2)[2], length(x_input.gg)), 
            rep(gg_color_hue(2)[1], length(x_new.gg))), 
  shape = c(rep(8, length(x_input.gg)), 
            rep(16, length(x_new.gg)))
)
ggplot(data = ggdata.melted, aes(x = x, y =value, color = variable), 
       linetype = 1) + 
  geom_path() + 
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = variable), 
              alpha = 0.1, linetype = 0) +
  scale_linetype_manual(values = c(1, 1, 2, 2)) + 
  scale_fill_manual(values = c(NA, "#00BFC4", "#C77CFF")) + 
  scale_color_manual(values = c(1, "#00BFC4", "#C77CFF")) + 
  geom_point(data = ggdata_pts, mapping = aes(x = x, y = y), 
             inherit.aes = FALSE, color = ggdata_pts$color, 
             shape = ggdata_pts$shape, 
             size = 3, alpha = 0.25) +
  geom_label(label = c(rep(NA, length(x_input.gg)), c(1:length(x_new.gg))), 
             data = ggdata_pts, mapping = aes(x = x, y = y + 0.3), 
             inherit.aes = FALSE, alpha = 0.25) +
  geom_point(data = ggdata_pts, mapping = aes(x = x, y = yrange[1]), 
             inherit.aes = FALSE, color = ggdata_pts$color, 
             shape = ggdata_pts$shape, 
             size = 3, alpha = 0.25) + 
  scale_y_continuous(limits = yrange) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  labs(y = "y", x = "x", fill = "Function", color = "Function")
# ggsave("BHgp_110520.pdf",
#        plot = last_plot(),
#        device = "pdf",
#        path = image_path,
#        scale = 1,
#        width = 6,
#        height = 4,
#        units = c("in")
# )

post.probs.gg = data.frame(
  x = 1:dim(BHres$post.probs)[1],
  H0 = BHres$post.probs[, 1], 
  H1 = BHres$post.probs[, 2]
)
post.probs.ggm = reshape2::melt(
  post.probs.gg, id.vars = c("x"))
ggplot(post.probs.ggm, aes(x = x, y = value)) + 
  facet_wrap(vars(variable)) + 
  geom_path() + 
  theme_bw() + 
  ylab("posterior probability of hypothesis")
# ggsave("BHgp_pph_110520.pdf",
#        plot = last_plot(),
#        device = "pdf",
#        path = image_path,
#        scale = 1,
#        width = 6,
#        height = 4,
#        units = c("in")
# )

# plot criterion for first point
BHcrit1 = BHDgp_m2(y_input, post_probs0, x_input, x_seq, model0, model1, nugget)
BHD.gg = data.frame(
  x = x_seq, 
  y = BHcrit1
)
ggplot(BHD.gg, aes(x = x, y = y)) + 
  geom_path() + 
  theme_bw() + 
  ylab("posterior probability of hypothesis")
ggplot(data = ggdata.melted, aes(x = x, y =value, color = variable), 
       linetype = 1) + 
  geom_path() + 
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = variable), 
              alpha = 0.1, linetype = 0) +
  scale_linetype_manual(values = c(1, 1, 2, 2)) + 
  scale_fill_manual(values = c(NA, "#00BFC4", "#C77CFF")) + 
  scale_color_manual(values = c(1, "#00BFC4", "#C77CFF")) + 
  geom_point(data = ggdata_pts, mapping = aes(x = x, y = y), 
             inherit.aes = FALSE, color = ggdata_pts$color, 
             shape = ggdata_pts$shape, 
             size = 3, alpha = 0.25) +
  geom_label(label = c(rep(NA, length(x_input.gg)), c(1:length(x_new.gg))), 
             data = ggdata_pts, mapping = aes(x = x, y = y + 0.3), 
             inherit.aes = FALSE, alpha = 0.25) +
  geom_point(data = ggdata_pts, mapping = aes(x = x, y = yrange[1]), 
             inherit.aes = FALSE, color = ggdata_pts$color, 
             shape = ggdata_pts$shape, 
             size = 3, alpha = 0.25) + 
  geom_path(data = BHD.gg, aes(x = x, y = (1/1000) * y - 3), inherit.aes = FALSE, color = 3) + 
  scale_y_continuous(limits = yrange) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  labs(y = "y", x = "x", fill = "Function", color = "Function")
# ggsave("BHDgp_110520.pdf",
#        plot = last_plot(),
#        device = "pdf",
#        path = image_path,
#        scale = 1,
#        width = 6,
#        height = 4,
#        units = c("in")
# )
