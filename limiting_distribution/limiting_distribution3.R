###########################
# Making W a Distribution #
###########################

# --- Directory --- #
# Kristyn
home = "/Users/kristyn/Documents/research/smed_ms/"

# --- Sources/Libraries --- #
functions_home = paste(home, "/functions", sep="")
source(paste(functions_home, "/wasserstein_distance.R", sep = ""))
source(paste(functions_home, "/charge_function_q.R", sep = ""))
source(paste(functions_home, "/variance_marginal_y.R", sep = ""))
source(paste(functions_home, "/generate_MED_oneatatime.R", sep = ""))
source(paste(functions_home, "/generate_MED_fast.R", sep = ""))
library(expm)
library(matrixStats)


#####

mu0 = c(0, 0)
mu1 = c(0, 0, 0)
typeT = 3
betaT = c(-0.2, -0.4, 0.4)
sigmasq01 = 0.25
sigmasq = 0.1

# MED design #
typeT = 3
xmin = -1
xmax = 1
# MED design #
V0 = diag(rep(sigmasq01,length(mu0)))
V1 = diag(rep(sigmasq01,length(mu1)))
# function settings (including and based on prior settings above)
f0 = function(x) mu0[1] + mu0[2] * x
f1 = function(x) mu1[1] + mu1[2] * x + mu1[3] * x^2
type01 = c(2, 3)
numCandidates = 10^3
k = 4
xmin = -1
xmax = 1
p = 2
N = 100

fT = function(x) betaT[1] + betaT[2] * x + betaT[3] * x^2
curve(f0, col = 2, lwd = 5, xlim = c(xmin, xmax), xlab = "x", ylab = "f")
curve(f1, add = T, col = 3, lty = 2, lwd = 5)
curve(fT, add = T, col = 1, lty = 3, lwd = 5)
legend("bottom", c("f0", "f1", "true f"), lty = c(1,2,3), lwd = 5, col = c(2, 3, 1))

# marginal distribution of y, variance
x_seq = seq(from = -1, to = 1, length.out = 1e3)
# H0:
var_margy0 = sapply(x_seq, FUN = function(x) var_marginaly(x, V0, sigmasq, type = 2))
f0_seq = sapply(x_seq, f0)
y0_seq_lower = f0_seq - var_margy0
y0_seq_upper = f0_seq + var_margy0
polygon(c(x_seq,rev(x_seq)),c(y0_seq_lower,rev(y0_seq_upper)),col=rgb(1, 0, 0, 0.1), border = NA)
# lines(x_seq, y0_seq_lower, lty = 3, col = 2, lwd = 5)
# lines(x_seq, y0_seq_upper, lty = 3, col = 2, lwd = 5)
# H1:
var_margy1 = sapply(x_seq, FUN = function(x) var_marginaly(x, V1, sigmasq, type = 3))
f1_seq = sapply(x_seq, f1)
y1_seq_lower = f1_seq - var_margy1
y1_seq_upper = f1_seq + var_margy1
polygon(c(x_seq,rev(x_seq)),c(y1_seq_lower,rev(y1_seq_upper)),col=rgb(0, 1, 0, 0.1), border = NA)
# lines(x_seq, y1_seq_lower, lty = 3, col = 3, lwd = 5)
# lines(x_seq, y1_seq_upper, lty = 3, col = 3, lwd = 5)

# MED design #
V0 = diag(rep(sigmasq01,length(mu0)))
V1 = diag(rep(sigmasq01,length(mu1)))
# function settings (including and based on prior settings above)
f0 = function(x) mu0[1] + mu0[2] * x
f1 = function(x) mu1[1] + mu1[2] * x + mu1[3] * x^2
type01 = c(2, 3)
numCandidates = 1000 #5000
k = 4
xmin = -1
xmax = 1
p = 2
N = 300

# alpha = 1, i.e. original, more space-filling-design
mmedv2 = MED_ms_oneatatime(mu0, mu1, V0, V1, sigmasq, f0, f1, type01, N, 
                           numCandidates, k, xmin, xmax, p = p, alpha = 1)




###
f_dens = function(x, mean_beta0, mean_beta1, var_mean0, var_mean1, var_e, f0, f1, type){
  if(length(type) != 2) stop("type should be vector with length == 2")
  mu1 = f0(x) # mean of marginal dist of y | H0
  mu2 = f1(x) # mean of marginal dist of y | H1
  var1 = var_marginaly(x, var_mean0, var_e, type = type[1]) # variance of marginal dist of y | H0
  var2 = var_marginaly(x, var_mean1, var_e, type = type[2]) # variance of marginal dist of y | H1
  Wass_dist = Wasserstein_distance(mu1, mu2, var1, var2)
  return(Wass_dist)
}

unnormalized_f_dens_curve = function(x) f_dens(x, mu0, mu1, V0, V1, sigmasq, f0, f1, type01)
integratefn_int_const = integrate(f = unnormalized_f_dens_curve, lower = -1, upper = 1)$value

f_dens_curve2 = function(x) (1/integratefn_int_const) * unnormalized_f_dens_curve(x)

###
hist(mmedv2, probability = T, ylim = c(0, 2.5), breaks = 20)
curve(f_dens_curve2, add = T)

saveRDS(mmedv2, "mmedv2_N300.rds")










