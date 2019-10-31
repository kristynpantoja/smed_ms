library(rodd)
library(mvtnorm)

sigmasq = 0.5
xmin = -1
xmax = 1
N = 50

betaT = c(-0.2, 0, 0.5) # like ex2.1
typeT = 3

numSims = 500
space_filling = seq(from = xmin, to = xmax, length.out = N)
ysims_space = simulateY(space_filling, N, betaT, sigmasq, numSims = numSims, typeT, seed = 123)

lm_lin_sims = matrix(NA, 2, numSims)
lm_quad_sims = matrix(NA, 3, numSims)
lm_cub_sims = matrix(NA, 4, numSims)
for(i in 1:numSims){
  lm_lin_sims[ , i] = lm(ysims_space[ , i] ~ space_filling)$coefficients
  lm_quad_sims[ , i] = lm(ysims_space[ , i] ~ space_filling + I(space_filling^2))$coefficients
  lm_cub_sims[ , i] = lm(ysims_space[ , i] ~ space_filling + I(space_filling^2) + 
                           I(space_filling^3))$coefficients
}
lm_lin = apply(lm_lin_sims, 1, mean)
lm_quad = apply(lm_quad_sims, 1, mean)
lm_cub = apply(lm_cub_sims, 1, mean)





# https://cran.r-project.org/web/packages/rodd/index.html
# install.packages("rodd")
library(rodd)
?rodd

### Auxiliary libraries for examples
library(mvtnorm)

#List of models
eta0 <- function(x, theta0) 
  theta0[1] + theta0[2] * x

eta1 <- function(x, theta1) 
  theta1[1] + theta1[2] * x + theta1[3] * x^2

eta <- list(eta0, eta1)

#List of fixed parameters
theta0 <- lm_lin
theta1 <- lm_quad
theta.fix <- list(theta0, theta1)

#Comparison table, holding quadratic (i = 2) model parameters fixed while letting linear model parameters vary
p <- matrix(c(0, 0, 
              1, 0), c(2, 2), byrow = TRUE)
#Design estimation - is the D-Optimal quadratic design... why?
res <- tpopt(x = space_filling, eta = eta, theta.fix = theta.fix, p = p, 
             x.lb = -1, x.rb = 1)


#Comparison table, holding linear (i = 1) model parameters fixed while letting quadratic model parameters vary
p2 <- matrix(c(0, 1, 
              0, 0), c(2, 2), byrow = TRUE)
#Design estimation - is close to D-Optimal quadratic design, but with different weights (more in middle), and slightly off-center
res2 <- tpopt(x = space_filling, eta = eta, theta.fix = theta.fix, p = p2, 
             x.lb = -1, x.rb = 1)








# https://cran.r-project.org/web/packages/rodd/index.html
# install.packages("rodd")
library(rodd)

### Auxiliary libraries for examples
library(mvtnorm)

#List of models
eta0 <- function(x, theta0) 
  theta0[1] + theta0[2] * x

eta1 <- function(x, theta1) 
  theta1[1] + theta1[2] * x + theta1[3] * x^2

eta2 <- function(x, theta2) 
  theta2[1] + theta2[2] * x + theta2[3] * x^2 + theta2[4] * x^3

eta <- list(eta0, eta1, eta2)

#List of fixed parameters
theta0 <- lm_lin
theta1 <- lm_quad
theta2 <- lm_cub
theta.fix <- list(theta0, theta1, theta2)

#Comparison table, holding quadratic (i = 2) model parameters fixed while letting cubic model parameters vary
p <- matrix(c(0, 0, 0,
              0, 0, 1,
              0, 0, 0), c(3, 3), byrow = TRUE)
#Design estimation - is the D-Optimal quadratic design... why?
res <- tpopt(x = space_filling, eta = eta, theta.fix = theta.fix, p = p, 
             x.lb = -1, x.rb = 1)


#Comparison table, holding cubic (i = 1) model parameters fixed while letting quadratic model parameters vary
p2 <- matrix(c(0, 0, 0,
              0, 0, 0,
              0, 1, 0), c(3, 3), byrow = TRUE)
#Design estimation - is the D-Optimal quadratic design... why?
res2 <- tpopt(x = space_filling, eta = eta, theta.fix = theta.fix, p = p2, 
             x.lb = -1, x.rb = 1)


















#######################
#######################

# https://cran.r-project.org/web/packages/rodd/index.html
# install.packages("rodd")
library(rodd)
?rodd

### Auxiliary libraries for examples
library(mvtnorm)

#List of models
eta0 <- function(x, theta0) 
  theta0[1] + theta0[2] * x

eta1 <- function(x, theta1) 
  theta1[1] + theta1[2] * x + theta1[3] * x^2

eta <- list(eta0, eta1)

#List of fixed parameters
theta0 <- c(1, 1)
theta1 <- c(1, 1, 1)
theta.fix <- list(theta0, theta1)

#Comparison table
p <- matrix(c(0, 0, 
              1, 0), c(2, 2), byrow = TRUE)

#Design estimation - is the D-Optimal design... why?
res <- tpopt(x = c(-1, 1), eta = eta, theta.fix = theta.fix, p = p, 
             x.lb = -1, x.rb = 1)
res <- tpopt(x = c(-1, 0, 1), eta = eta, theta.fix = theta.fix, p = p, 
             x.lb = -1, x.rb = 1)

plot(res)
summary(res)











