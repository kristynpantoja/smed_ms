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
<<<<<<< HEAD
#Design estimation - is the D-Optimal design... why?
res <- tpopt(x = c(-1, 1), eta = eta, theta.fix = theta.fix, p = p, 
             x.lb = -1, x.rb = 1)
res <- tpopt(x = c(-1, 0, 1), eta = eta, theta.fix = theta.fix, p = p, 
             x.lb = -1, x.rb = 1)




#Comparison table
p <- matrix(c(0, 1, 
              0, 0), c(2, 2), byrow = TRUE)
=======

>>>>>>> 215423400ca954fbc720cc0dbcdc4b6a76225f70
#Design estimation - is the D-Optimal design... why?
res <- tpopt(x = c(-1, 1), eta = eta, theta.fix = theta.fix, p = p, 
             x.lb = -1, x.rb = 1)
res <- tpopt(x = c(-1, 0, 1), eta = eta, theta.fix = theta.fix, p = p, 
             x.lb = -1, x.rb = 1)

plot(res)
summary(res)

















### EMAX vs MM
#List of models: model1 is larger than model 2
eta.1 <- function(x, theta.1) 
  theta.1[1] + theta.1[2] * x / (x + theta.1[3])

eta.2 <- function(x, theta.2) 
  theta.2[1] * x / (x + theta.2[2])

eta <- list(eta.1, eta.2)

#List of fixed parameters
theta.1 <- c(1, 1, 1)
theta.2 <- c(1, 1)
theta.fix <- list(theta.1, theta.2)

#Comparison table
p <- matrix(
  c(
    0, 1,
    0, 0
  ), c(2, 2), byrow = TRUE)
# here model 1 is fixed, model 2 is allowed to vary
# since the fixed model, model 1, is smaller, 
#   must also fix values for shared parameter
#Design estimation
res <- tpopt(x = c(1.2, 1.5, 1.7), eta = eta, theta.fix = theta.fix, p = p, 
             x.lb = 1, x.rb = 2)



#Comparison table
p <- matrix(
  c(
    0, 0,
    1, 0
  ), c(2, 2), byrow = TRUE)
# here model 2 is fixed, model 1 is allowed to vary
#Design estimation
res <- tpopt(x = c(1.2, 1.5, 1.7), eta = eta, theta.fix = theta.fix, p = p, 
             x.lb = 1, x.rb = 2)
