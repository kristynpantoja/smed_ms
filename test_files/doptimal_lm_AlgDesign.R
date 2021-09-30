# candidates
numCandidates = 101 # size doesn't matter when using approximate = TRUE
xmin=-1; xmax=1
x_seq = seq(from = xmin, to = xmax, length.out = numCandidates)

library(AlgDesign)
candidates = data.frame(x = x_seq)

## approximate = TRUE ###
# outputs proportions for each support point in the design

# consider linear model, y = b0 + b1 x
res_Fed_Doptlin = optFederov(
  ~1+x, data = candidates, approximate = TRUE, criterion = "D")
res_Fed_Doptlin
# consider quadratic model, y = b0 + b1 x + b2 x^2
res_Fed_Doptquad = optFederov(
  ~1+x+I(x^2), data = candidates, approximate = TRUE, criterion = "D")
res_Fed_Doptquad

## approximate = FALSE ###
# outputs sample size for each support point in the design
# doesn't work great because doesn't replace a candidate point after choosing it

# # consider linear model, y = b0 + b1 x
# res_Fed_Doptlin_N30 = optFederov(
#   ~1+x, data = candidates, approximate = FALSE, criterion = "D", 
#   nTrials = 30)
# res_Fed_Doptlin_N30
# # consider quadratic model, y = b0 + b1 x + b2 x^2
# res_Fed_Doptquad_N30 = optFederov(
#   ~1+x+I(x^2), data = candidates, approximate = FALSE, criterion = "D",
#   nTrials = 30)
# res_Fed_Doptquad_N30
