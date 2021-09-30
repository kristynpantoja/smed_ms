# candidates
dimX = 3
numCandidates = 5000
xmin=-1; xmax=1
x_seq = seq(from = xmin, to = xmax, length.out = floor((numCandidates)^(1 / 3)))
candidates = as.matrix(expand.grid(x_seq, x_seq, x_seq))
colnames(candidates) = paste("x", 1:dimX, sep = "")

################################################################################
library(OptimalDesign)

# nxm matrix F, denoted Fx, of all candidate regressors
# the result is equivalent to candidates, but with named variables
Fx <- Fx_cube(
  ~-1 + x1 + x2 + x3, lower = rep(-1, dimX), upper = rep(1, dimX),
  n.levels = rep(floor((numCandidates)^(1 / 3)), dimX))

# different algorithms for computing optimal or nearly-optimal approximate or 
#   exact experimental design

# od_REX computes optimal approx design under standard (size) constraint using
#   one of 3 methods (alg.AA = "REX", "MUL", "VDM")
# res_REX = od_REX(candidates, crit = "D")
# design_REX = as.data.frame(candidates)
# design_REX$freq = round(res_REX$w.best, 2)
# design_REX = design_REX[design_REX$freq > 0, ]
# design_REX

# od_KL computes optimal or near-optimal exact design of experiments under std
#   constraint
res_KL = od_KL(Fx, N = 24, crit = "D")
design_KL = as.data.frame(Fx)
design_KL$freq = res_KL$w.best
design_KL = design_KL[design_KL$freq > 0, ]
design_KL

################################################################################
library(AlgDesign)

res_Fed = optFederov(
  ~-1+x1+x2+x3, data = candidates, nTrials = 24, approximate = FALSE, 
  criterion = "D")
