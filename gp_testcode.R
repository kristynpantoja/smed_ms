library(MASS)
library(reshape2)
library(ggplot2)

library(transport)
library(mined)



# --- Begin Helper Functions --- #

# Covariance function
#  here, used radial aka squared exponential aka gaussian
C_fn_elementwise = function(Xi,Xj, l) exp(-0.5 * (Xi - Xj) ^ 2 / l ^ 2) 
#C_fn = function(X, Y) outer(X, Y, FUN = function(X, Y) C_fn_elementwise(X, Y, l))
C_fn = function(X, Y, l) outer(X, Y, FUN = C_fn_elementwise, l)
#C_fn2 = function(X, Y) outer(X, Y, C_fn_elementwise, l)

# --- End Helper Functions --- #





# --- Begin GP Notes from Example --- #


# GP with Squared Exponential C
# resource for getting started: 
#  https://www.carlboettiger.info/2012/10/17/basic-regression-in-gaussian-processes.html



# GP Parameter Specifications
l = 0.2 # scale parameter for radial covariance function
# observation: the smaller the "l" the more wiggly the plot.



## With Data

# In general we aren’t interested in drawing from the prior,
#  but want to include information from the data as well.
#  Here, we assume no noise.

# data:
obs = data.frame(x = c(-4, -3, -1,  0,  2), y = c(-2,  0,  1,  2, -1))

# new design points
xstar = seq(-5,5,len=50)
# parameter
l = 1
#C_fn = function(X, Y, l) outer(X, Y, FUN = C_fn_elementwise, l)



## draws from null distribution (prior at those points xstar)

# Get covariance matrix from Covariance function
covariance_xstar = C_fn(xstar, xstar, l)

# Generate a number of functions from the GP, mean = 0 for now.
numFns = 2
set.seed(1234)
values <- mvrnorm(numFns, rep(0, length=length(xstar)), covariance_xstar)

# R eshape the data into long (tidy) form, listing x value, y value, and sample number
data <- data.frame(x=xstar, t(values))
data <- melt(data, id="x")
head(data)

# Plot
fig1 <- ggplot(data,aes(x=x,y=value)) +
  #geom_rect(xmin=-Inf, xmax=Inf, ymin=-2, ymax=2, fill="grey80") +
  geom_line(aes(group=variable)) +   theme_bw() +
  scale_y_continuous(lim=c(-3,3), name="output, f(x)") +
  xlab("input, x")
fig1



# Assume y ~ N(0, cov(X, X))

# Then the conditional probability of observing our data is
#  y* | X*, X, y ~ N(K(X*, X) K(X, X)^(-1) y, 
#                    K(X*, X*) - K(X*, X) K(X, X)^(-1) K(X, X*))

# Calculate posterior mean and covariance:
cov_xx_inv <- solve(C_fn(obs$x, obs$x, l))
Ef <- C_fn(xstar, obs$x, l) %*% cov_xx_inv %*% obs$y
Cf <- C_fn(xstar, xstar, l) - C_fn(xstar, obs$x, l)  %*% cov_xx_inv %*% C_fn(obs$x, xstar, l)

# Draw random samples from posterior
set.seed(1234)
values <- mvrnorm(3, Ef, Cf)

dat = data.frame(x=xstar, t(values))
dat = melt(dat, id="x")

dat2 = data.frame(xstar, Ef)

fig2 = ggplot(dat,aes(x=dat$x,y=dat$value)) +
  geom_ribbon(data=dat2, 
              aes(x=xstar, y=Ef, ymin=(Ef-2*sqrt(diag(Cf))), 
                  ymax=(Ef+2*sqrt(diag(Cf)))),
              fill="grey80") +
  geom_line(aes(color=variable)) + #REPLICATES
  geom_line(data=dat2,aes(x=xstar,y=Ef), size=1) + #MEAN
  geom_point(data=obs,aes(x=x,y=y)) +  #OBSERVED DATA
  scale_y_continuous(lim=c(-3,3), name="output, f(x)") +
  xlab("input, x")
fig2



## Noise in Data

# y = f(x) + epsilon (i.e. some process noise)
# epsilon ~ Normal(0, sigma_eps)

sigma_eps <- 0.25
cov_xx_inv_noise <- solve(C_fn(obs$x, obs$x, l) + sigma_eps^2 * diag(1, length(obs$x)))
E_f_noise <- C_fn(xstar, obs$x, l) %*% cov_xx_inv_noise %*% obs$y
C_f_noise <- C_fn(xstar, xstar, l) - C_fn(xstar, obs$x, l)  %*% cov_xx_inv_noise %*% C_fn(obs$x, xstar, l)

# 3 samples from posterior 
set.seed(1234)
values_noise <- mvrnorm(3, E_f_noise, C_f_noise)

# plot

dat_noise <- data.frame(x=xstar, t(values_noise))
dat_noise <- melt(dat_noise, id="x")

dat_noise2 = data.frame(xstar, E_f_noise)

fig3 <- ggplot(dat_noise,aes(x=x,y=value)) +
  geom_ribbon(data=dat_noise2, 
              aes(x=xstar, y=E_f_noise, ymin=(E_f_noise-2*sqrt(diag(C_f_noise))), 
                  ymax=(E_f_noise+2*sqrt(diag(C_f_noise)))),
              fill="grey80") + # Var
  geom_line(aes(color=variable)) + #REPLICATES
  geom_line(data=dat_noise2,aes(x=xstar,y=E_f_noise), size=1) + #MEAN
  geom_point(data=obs,aes(x=x,y=y)) +  #OBSERVED DATA
  scale_y_continuous(lim=c(-3,3), name="output, f(x)") +
  xlab("input, x")
fig3

# --- End GP Notes from Example --- #





# --- Begin GP Other Notes / Playing Around --- #



## choose appropriate parameter l and other stuff for points X in [0, 1]

############
## Drawing from GP Prior (i.e. null distribution), for X between 0 and 1

# Try some X
X = seq(from = 0, to = 1, length.out = 51)

# GP Parameter Specifications
l = 0.2 # scale parameter for radial covariance function
# observation: the smaller the "l" the more wiggly the plot.

# Covariance function
#  here, used radial aka squared exponential aka gaussian
C_fn_elementwise = function(Xi,Xj, l) exp(-0.5 * (Xi - Xj) ^ 2 / l ^ 2) 
#C_fn = function(X, Y) outer(X, Y, FUN = function(X, Y) C_fn_elementwise(X, Y, l))
C_fn = function(X, Y, l) outer(X, Y, FUN = C_fn_elementwise, l)
#C_fn2 = function(X, Y) outer(X, Y, C_fn_elementwise, l)

# Get covariance matrix from Covariance function
covariance = C_fn(X, X, l)

# Generate a number of functions from the GP, mean = 0 for now.
numFns = 2
set.seed(1234)
values <- mvrnorm(numFns, rep(0, length=length(X)), covariance)

# R eshape the data into long (tidy) form, listing x value, y value, and sample number
data <- data.frame(x=X, t(values))
data <- melt(data, id="x")
head(data)

# Plot
fig1 <- ggplot(data,aes(x=x,y=value)) +
  #geom_rect(xmin=-Inf, xmax=Inf, ymin=-2, ymax=2, fill="grey80") +
  geom_line(aes(group=variable)) +   theme_bw() +
  scale_y_continuous(lim=c(-3,3), name="output, f(x)") +
  xlab("input, x")
fig1

# this is a plot of draws from the null distribution.

#data1 = subset(data, variable == "X1")
#data2 = subset(data, variable == "X2")
#plot(x = data1$x, data1$value, type = "l")
############


# ...







## GPs Different l 


# Try some X
X = seq(from = 0, to = 1, length.out = 51)

# GP Parameter Specifications
l1 = 0.2 # scale parameter for radial covariance function
l2 = 0.3
# observation: the smaller the "l" the more wiggly the plot.

# Covariance function
# still using radial aka squared exponential aka gaussian
C_fn_elementwise = function(Xi,Xj, l) exp(-0.5 * (Xi - Xj) ^ 2 / l ^ 2) 
C_fn = function(X, Y, l) outer(X, Y, FUN = C_fn_elementwise, l)

# Get covariance matrix from Covariance function
cov1 = C_fn(X, X, l1)
cov2 = C_fn(X, X, l2)

# Generate a fn from each GP
set.seed(1234)
values1 <- mvrnorm(1, rep(0, length=length(X)), cov1)
values2 <- mvrnorm(1, rep(0, length=length(X)), cov2)

# R eshape the data into long (tidy) form, listing x value, y value, and sample number
data1 = data.frame(x = X, values = values1)
data2 = data.frame(x = X, values = values2)

# Plot
ymin = min(c(data1$values, data2$values))
ymax = max(c(data1$values, data2$values))
plot(data1, type = "l", ylim = c(ymin, ymax))
lines(data2)

















