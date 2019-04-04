############################################################
## Based on Fast Algorithm, Joseph et. al. 2018
############################################################

### --- Linear Regression Different Slopes (both intercept at 0, same error variance) --- ###

## Like bayes_linear_regression.R, but I try to implement some of the techniques in
##  "Deterministic Sampling of Expensive Posteriors Using Minimum Energy Designs,"
##  Joseph 2018

library(transport)
library(mined)
library(expm)

source("smed_ms_functions.R")



##########################



### Testing



##########################



## Pick Parameters

mean_beta0 = 1 # slope of null model
mean_beta1 = 1 / 2 # slope of alternative model
var_mean = 0.05
var_e = 0.1 # same variance

p = 1
xmin = 0
xmax = 1

# parameters we change sometimes:
N = 51 # in paper, n = numCandidates - not true, it's numCandidates generated for each x_i^k at each step
numCandidates = 7 # largest prime number less than 100 + 5p = 103
K = 4 # ceiling(4* sqrt(p))



############################
## Running Algorithm, K = 4
############################

#### N = 11


N = 11
X_test = SMED_ms_fast(mean_beta0, mean_beta1, var_e, var_mean, N, xmin, xmax, K, p)
#X_test_D_sorted = apply(X_test$D, 2, sort)
#X_test$candidates

f0 = function(x) mean_beta0 * x # null regression model
f1 = function(x) mean_beta1 * x # alternative regression model

test_k = 4
X_k = sort(X_test$D[ ,test_k])
curve(f0, col = 1, from = xmin, to = xmax, xlab = "design points", ylab = "f1, f2")
curve(f1, col = 1, add = TRUE)
text(X_test$D[ ,test_k], f0(X_test$D[ ,test_k]), c(1:N), col=4)
text(X_test$D[ ,test_k], f1(X_test$D[ ,test_k]), c(1:N), col=4)
points(X_k, rep(0, N), col = 2)

#dev.copy(png, 'fast_N11.png')
#dev.off()

X_k = sort(X_test$D[ ,test_k])
curve(f0, col = 1, from = xmin, to = xmax, xlab = "design points", ylab = "f1, f2")
curve(f1, col = 1, add = TRUE)
for(i in 1:N){
  text(X_test$D[i ,test_k], f1(X_test$D[i ,test_k]) + i * 0.01, i, col=4)
}
points(X_k, rep(0, N), col = 2)

#dev.copy(png, 'fast_pattern_N11.png')
#dev.off()


#### N = 67
# does not vary by much

N = 67
X_test = SMED_ms_fast(mean_beta0, mean_beta1, var_e, var_mean, N, xmin, xmax, K, p)

test_k = 4
X_k = sort(X_test$D[ ,test_k])
curve(f0, col = 1, from = xmin, to = xmax, xlab = "design points", ylab = "f1, f2")
curve(f1, col = 1, add = TRUE)
text(X_test$D[ ,test_k], f0(X_test$D[ ,test_k]), c(1:N), col=4)
text(X_test$D[ ,test_k], f1(X_test$D[ ,test_k]), c(1:N), col=4)
points(X_k, rep(0, N), col = 2)

#dev.copy(png, 'fast_N67.png')
#dev.off()

X_k = sort(X_test$D[ ,test_k])
curve(f0, col = 1, from = xmin, to = xmax, xlab = "design points", ylab = "f1, f2")
curve(f1, col = 1, add = TRUE)
for(i in 1:N){
  text(X_test$D[i ,test_k], f1(X_test$D[i ,test_k]) + i * 0.01, i, col=4)
}
points(X_k, rep(0, N), col = 2)

#dev.copy(png, 'fast_pattern_N67.png')
#dev.off()


##########################


# What if we use uniform, instead of Lattice?
# high variability!!! why???

N = 67
X_test = SMED_ms_fast2(mean_beta0, mean_beta1, var_e, var_mean, N, xmin, xmax, K, p)

test_k = 4
X_k = sort(X_test$D[ ,test_k])
curve(f0, col = 1, from = xmin, to = xmax, xlab = "design points", ylab = "f1, f2")
curve(f1, col = 1, add = TRUE)
text(X_test$D[ ,test_k], f0(X_test$D[ ,test_k]), c(1:N), col=4)
text(X_test$D[ ,test_k], f1(X_test$D[ ,test_k]), c(1:N), col=4)
points(X_k, rep(0, N), col = 2)

#dev.copy(png, 'fastunif_N67.png')
#dev.off()

X_k = sort(X_test$D[ ,test_k])
curve(f0, col = 1, from = xmin, to = xmax, xlab = "design points", ylab = "f1, f2")
curve(f1, col = 1, add = TRUE)
for(i in 1:N){
  text(X_test$D[i ,test_k], f1(X_test$D[i ,test_k]) + i * 0.01, i, col=4)
}
points(X_k, rep(0, N), col = 2)

#dev.copy(png, 'fastunif_pattern_N67.png')
#dev.off()










############################
## Running Algorithm, K = 20
############################



#### N = 20

K = 20
N = 11
X_test = SMED_ms_fast(mean_beta0, mean_beta1, var_e, var_mean, N, xmin, xmax, K, p)
#X_test_D_sorted = apply(X_test$D, 2, sort)
#X_test$candidates

f0 = function(x) mean_beta0 * x # null regression model
f1 = function(x) mean_beta1 * x # alternative regression model

test_k = K
X_k = sort(X_test$D[ ,test_k])
curve(f0, col = 1, from = xmin, to = xmax, xlab = "design points", ylab = "f1, f2")
curve(f1, col = 1, add = TRUE)
text(X_test$D[ ,test_k], f0(X_test$D[ ,test_k]), c(1:N), col=4)
text(X_test$D[ ,test_k], f1(X_test$D[ ,test_k]), c(1:N), col=4)
points(X_k, rep(0, N), col = 2)

#dev.copy(png, 'fast_N11_K20.png')
#dev.off()



#### N = 67

K = 20
N = 67
X_test = SMED_ms_fast(mean_beta0, mean_beta1, var_e, var_mean, N, xmin, xmax, K, p)

test_k = K

curve(f0, col = 1, from = xmin, to = xmax, xlab = "design points", ylab = "f1, f2", ylim = c(0, 3))
curve(f1, col = 1, add = TRUE)
y = rep(NA, N)
for(i in 1:N){
  y[i] = i * 0.04
  text(X_test$D[i ,test_k], y[i], i, col=4)
}
points(X_test$D[ ,test_k], rep(0, N), col = 2)
lines(X_test$D[ ,test_k], y, col = 3)

#dev.copy(png, 'fast_pattern_N67_K20.png')
#dev.off()











############################
## Running Algorithm, K = 40
############################



#### N = 11

K = 100
N = 11
X_test = SMED_ms_fast(mean_beta0, mean_beta1, var_e, var_mean, N, xmin, xmax, K, p)
#X_test_D_sorted = apply(X_test$D, 2, sort)
#X_test$candidates

f0 = function(x) mean_beta0 * x # null regression model
f1 = function(x) mean_beta1 * x # alternative regression model

test_k = K
X_k = sort(X_test$D[ ,test_k])
curve(f0, col = 1, from = xmin, to = xmax, xlab = "design points", ylab = "f1, f2")
curve(f1, col = 1, add = TRUE)
text(X_test$D[ ,test_k], f0(X_test$D[ ,test_k]), c(1:N), col=4)
text(X_test$D[ ,test_k], f1(X_test$D[ ,test_k]), c(1:N), col=4)
points(X_k, rep(0, N), col = 2)

#dev.copy(png, 'fast_N11_K100.png')
#dev.off()



#### N = 67


K = 100
N = 67
X_test = SMED_ms_fast(mean_beta0, mean_beta1, var_e, var_mean, N, xmin, xmax, K, p)

test_k = 100

X_k = sort(X_test$D[ ,test_k])
curve(f0, col = 1, from = xmin, to = xmax, xlab = "design points", ylab = "f1, f2", ylim = c(0.0, 2.5))
curve(f1, col = 1, add = TRUE)
for(i in 1:N){
  text(X_test$D[i ,test_k], f1(X_test$D[i ,test_k]) + i * 0.03, i, col=4)
}
points(X_k, rep(0, N), col = 2)

ks = seq(from = 0, to = 100, by = 10)
ks[1] = 1
for(i in 1:length(ks)){
  name = paste("fast_pattern_N67_K100_k", ks[i], ".png", sep = "")
  
  test_k = ks[i]
  X_k = sort(X_test$D[ ,test_k])
  curve(f0, col = 1, from = xmin, to = xmax, xlab = "design points", ylab = "f1, f2", 
        ylim = c(0.0, 2.5))
  curve(f1, col = 1, add = TRUE)
  for(i in 1:N){
    text(X_test$D[i ,test_k], f1(X_test$D[i ,test_k]) + i * 0.03, i, col=4)
  }
  points(X_k, rep(0, N), col = 2)
  dev.copy(png, name)
  dev.off()
}









#############################################
## Running Algorithm, K = 40, Lattice Version
#############################################



#### N = 11

K = 100
N = 11
X_test = SMED_ms_fast3(mean_beta0, mean_beta1, var_e, var_mean, N, xmin, xmax, K, p)
#X_test_D_sorted = apply(X_test$D, 2, sort)
#X_test$candidates

f0 = function(x) mean_beta0 * x # null regression model
f1 = function(x) mean_beta1 * x # alternative regression model

test_k = K
X_k = sort(X_test$D[ ,test_k])
curve(f0, col = 1, from = xmin, to = xmax, xlab = "design points", ylab = "f1, f2")
curve(f1, col = 1, add = TRUE)
text(X_test$D[ ,test_k], f0(X_test$D[ ,test_k]), c(1:N), col=4)
text(X_test$D[ ,test_k], f1(X_test$D[ ,test_k]), c(1:N), col=4)
points(X_k, rep(0, N), col = 2)

#dev.copy(png, 'fastlattice_N11_K100.png')
#dev.off()



#### N = 67


K = 40
N = 67
X_test = SMED_ms_fast3(mean_beta0, mean_beta1, var_e, var_mean, N, xmin, xmax, K, p)

test_k = K

X_k = sort(X_test$D[ ,test_k])
curve(f0, col = 1, from = xmin, to = xmax, xlab = "design points", ylab = "f1, f2", ylim = c(0.0, 2.5))
curve(f1, col = 1, add = TRUE)
for(i in 1:N){
  text(X_test$D[i ,test_k], f1(X_test$D[i ,test_k]) + i * 0.03, i, col=4)
}
points(X_k, rep(0, N), col = 2)

# looks similar to seq version

ks = seq(from = 0, to = 100, by = 10)
ks[1] = 1
for(i in 1:length(ks)){
  name = paste("fast_pattern_N67_K100_k", ks[i], ".png", sep = "")
  
  test_k = ks[i]
  X_k = sort(X_test$D[ ,test_k])
  curve(f0, col = 1, from = xmin, to = xmax, xlab = "design points", ylab = "f1, f2", 
        ylim = c(0.0, 2.5))
  curve(f1, col = 1, add = TRUE)
  for(i in 1:N){
    text(X_test$D[i ,test_k], f1(X_test$D[i ,test_k]) + i * 0.03, i, col=4)
  }
  points(X_k, rep(0, N), col = 2)
  dev.copy(png, name)
  dev.off()
}




















##########################

# What if we use uniform, instead of Lattice?
# high variability!!! why???

# with K = 4, 20, variability is high. What about 40? Takes a long time.
K = 40 # didn't run --- took too long!
N = 67
X_test = SMED_ms_fast2(mean_beta0, mean_beta1, var_e, var_mean, N, xmin, xmax, K, p)

X_k = sort(X_test$D[ ,test_k])
curve(f0, col = 1, from = xmin, to = xmax, xlab = "design points", ylab = "f1, f2")
curve(f1, col = 1, add = TRUE)
for(i in 1:N){
  text(X_test$D[i ,test_k], f1(X_test$D[i ,test_k]) + i * 0.01, i, col=4)
}
points(X_k, rep(0, N), col = 2)

#dev.copy(png, 'fastunif_pattern_N67_K20.png')
#dev.off()