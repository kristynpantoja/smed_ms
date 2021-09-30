# candidates
numCandidates = 5000
xmin=-1; xmax=1
x_seq = seq(from = xmin, to = xmax, length.out = floor((numCandidates)^(1 / 3)))
candidates = as.matrix(expand.grid(x_seq, x_seq, x_seq))
dimX = ncol(candidates)
colnames(candidates) = paste("x", 1:dimX, sep = "")

library(AlgDesign)

# consider linear model y = b1 x1 + b2 x2 + b3 x3
res_Fed_Doptlin = optFederov(
  ~-1+x1+x2+x3, data = candidates, approximate = TRUE, criterion = "D")
res_Fed_Doptlin
# consider linear model with 1st order interactions y = x1+x2+x3+x1x2+x1x3+x2x3
res_Fed_Doptquad = optFederov(
  ~-1+x1+x2+x3+x1*x2+x1*x3+x2*x3, data = candidates, approximate = TRUE, 
  criterion = "D")
res_Fed_Doptquad
