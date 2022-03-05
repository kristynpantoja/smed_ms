d_seq = seq(0, 0.1, length.out = 101)
lam = 0.01
period = 0.05

se_fn = function(d) exp(-d^2 / (2 * lam^2))
ma_fn = function(d) (1 + (sqrt(3) * d / lam)) * exp(-(sqrt(3) * d) / lam)
pd_fn = function(d) exp(-2 * (sin(pi * d / period))^2 / lam^2)
dat = data.frame(
  d = d_seq, 
  SE = se_fn(d_seq), 
  Matern = ma_fn(d_seq), 
  Periodic = pd_fn(d_seq)
)
dat.mlt = reshape2::melt(dat, id.vars = "d")
library(ggplot2)
ggplot(dat.mlt, aes(x = d, y = value, color = variable)) + 
  geom_path() + 
  labs(x = "d", y = element_blank())

