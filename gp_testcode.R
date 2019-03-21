library(MASS)
library(reshape2)
library(ggplot2)

library(transport)
library(mined)


# --- GP with Squared Exponential C ---  #

# resource for getting started: https://www.carlboettiger.info/2012/10/17/basic-regression-in-gaussian-processes.html

## Drawing from GP Prior

# Try some X
X = seq(from = 0, to = 1, length.out = 51)

# GP Parameter Specifications
l = 0.2 # scale parameter for radial covariance function
# observation: the smaller the "l" the more wiggly the plot.

# Covariance function
C_fn_elementwise = function(Xi,Xj, l) exp(-0.5 * (Xi - Xj) ^ 2 / l ^ 2) # radial aka squared exponential aka gaussian
C_fn = function(X, Y) outer(X, Y, FUN = function(X, Y) C_fn_elementwise(X, Y, l))
#C_fn2 = function(X, Y) outer(X, Y, C_fn_elementwise, l)

# Get covariance matrix from Covariance function
covariance = C_fn(X, X)

# Generate a number of functions from the GP, mean = 0 for now.
numFns = 2
set.seed(1234)
values <- mvrnorm(numFns, rep(0, length=length(X)), covariance)

# R eshape the data into long (tidy) form, listing x value, y value, and sample number
data <- data.frame(x=X, t(values))
data <- melt(data, id="x")
head(data)

# Plot
fig2a <- ggplot(data,aes(x=x,y=value)) +
  #geom_rect(xmin=-Inf, xmax=Inf, ymin=-2, ymax=2, fill="grey80") +
  geom_line(aes(group=variable)) +   theme_bw() +
  scale_y_continuous(lim=c(-3,3), name="output, f(x)") +
  xlab("input, x")
fig2a

# this is a plot of draws from the null distribution.

#data1 = subset(data, variable == "X1")
#data2 = subset(data, variable == "X2")
#plot(x = data1$x, data1$value, type = "l")

## Incorporating Data

# In general we aren’t interested in drawing from the prior,
#  but want to include information from the data as well.

# data:
obs = data.frame(x = c(-4, -3, -1,  0,  2), y = c(-2,  0,  1,  2, -1))

# new design points
xstar = seq(-5,5,len=50)
# parameter
l = 1

# Assume y ~ N(0, cov(X, X))

# Then the conditional probability of observing our data is
#  y* | X*, X, y ~ N(K(X*, X) K(X, X)^(-1) y, K(X*, X*) - K(X*, X) K(X, X)^(-1) K(X, X*))

# Calculate posterior mean and covariance:
cov_xx_inv <- solve(C_fn(obs$x, obs$x, l))
Ef <- C_fn(xstar, obs$x, l) %*% cov_xx_inv %*% obs$y
Cf <- C_fn(xstar, xstar, l) - C_fn(xstar, obs$x, l)  %*% cov_xx_inv %*% C_fn(obs$x, xstar, l)

# Draw rndom samples from posterior
values <- mvrnorm(3, Ef, Cf)

dat = data.frame(x=xstar, t(values))
dat = melt(dat, id="x")

dat2 = data.frame(xstar, Ef)

fig2b = ggplot(dat,aes(x=dat$x,y=dat$value)) +
  geom_ribbon(data=dat2, 
              aes(x=xstar, y=Ef, ymin=(Ef-2*sqrt(diag(Cf))), ymax=(Ef+2*sqrt(diag(Cf)))),
              fill="grey80") +
  geom_line(aes(color=variable)) + #REPLICATES
  geom_line(data=dat2,aes(x=xstar,y=Ef), size=1) + #MEAN
  geom_point(data=obs,aes(x=x,y=y)) +  #OBSERVED DATA
  scale_y_continuous(lim=c(-3,3), name="output, f(x)") +
  xlab("input, x")
fig2b



# --- GPs Different l --- #


# Try some X
X = seq(from = 0, to = 1, length.out = 51)

# GP Parameter Specifications
l1 = 0.2 # scale parameter for radial covariance function
l2 = 0.3
# observation: the smaller the "l" the more wiggly the plot.

# Covariance function
C_fn_elementwise = function(Xi,Xj, l) exp(-0.5 * (Xi - Xj) ^ 2 / l ^ 2) # radial aka squared exponential aka gaussian
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



# now the question is how to pick our x's in X to distinguish these 2 different GPs

#













