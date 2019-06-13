






mean_beta0 = c(0, 1) # slope of null model
mean_beta1 = c(0, 1/2) # slope of alternative model
var_mean = diag(c(0.001, 0.001)) # variance on beta
var_e = 0.01 # variance on error

xmin = 0 
xmax = 1 

f0 = function(x) mean_beta0[1] + mean_beta0[2] * x # null regression model
f1 = function(x) mean_beta1[1] + mean_beta1[2] * x # alternative regression model

N = 67
# for one-at-a-time algorithm:
numCandidates = 1000 # suggested 10^5

type = c(2, 2)

###############################
# One-at-a-Time Algorithm k = 4
k = 4
X_1atatime2 = SMED_ms(mean_beta0, mean_beta1, var_mean0, var_mean1, var_e, 
                      f0 = f0, f1 = f1, type = type, var_margy0 = NULL, var_margy1 = NULL, 
                      N = N, numCandidates = numCandidates, k = k, xmin = xmin, xmax = xmax, 
                      genCandidates = 1, initialpt = 1)

# design points locations
curve(f0, col = 1, from = xmin, to = xmax, xlab = "", ylab = "", ylim = c(0, 3), axes = F, main = "1atT k=4")
curve(f1, col = 1, add = TRUE)
axis(1)
y = rep(NA, N)
for(i in 1:N){
  y[i] = i * 0.04
  text(X_1atatime2[i], y[i], i, col=4)
}
points(X_1atatime2, rep(0, N), col = 2, pch = "*")
lines(X_1atatime2, y, col = 3)
hist(X_1atatime2, freq = T, main = "1atT k=4")
###############################
# Fast Algorithm, K = 20
K = 20

X_fast2 = SMED_ms_fast(mean_beta0, mean_beta1, var_mean0, var_mean1, var_e, 
                       f0 = f0, f1 = f1, type = type, var_margy0 = NULL, var_margy1 = NULL, 
                       N = N, K = K, p = 2, xmin = xmin, xmax = xmax, 
                       genCandidates = 1, initialpt = 1)

# design points locations
X_fast2_K = X_fast2$D[ , K]
curve(f0, col = 1, from = xmin, to = xmax, xlab = "", ylab = "", ylim = c(0, 3), axes = F, main = "1atT k=4")
curve(f1, col = 1, add = TRUE)
axis(1)
y = rep(NA, N)
for(i in 1:N){
  y[i] = i * 0.04
  text(X_fast2_K[i], y[i], i, col=4)
}
points(X_fast2_K, rep(0, N), col = 2, pch = "*")
lines(X_fast2_K, y, col = 3)
hist(X_fast2_K, freq = T, main = "1atT k=4")




mean_beta0 = c(0, 1) # slope of null model
mean_beta1 = c(1, 4, -4) # slope of alternative model
var_mean0 = diag(c(0.001, 0.001)) # variance on beta0
var_mean1 = diag(c(0.001, 0.001, 0.001)) # variance on beta1
var_e = 0.01 # variance on error

xmin = 0 
xmax = 1 

f0 = function(x) mean_beta0[1] + mean_beta0[2] * x # null regression model
f1 = function(x) mean_beta1[1] + mean_beta1[2] * x + mean_beta1[3] * x^2 # alternative regression model

N = 67
# for one-at-a-time algorithm:
numCandidates = 1000 # suggested 10^5

type = c(2, 3)
###############################
# One-at-a-Time Algorithm k = 4
k = 4
X_1atatime3 = SMED_ms(mean_beta0, mean_beta1, var_mean0, var_mean1, var_e, 
                      f0 = f0, f1 = f1, type = type, var_margy0 = NULL, var_margy1 = NULL, 
                      N = N, numCandidates = numCandidates, k = k, xmin = xmin, xmax = xmax, 
                      genCandidates = 1, initialpt = 1)
# design points locations
curve(f0, col = 1, from = xmin, to = xmax, xlab = "", ylab = "", ylim = c(0, 3), axes = F, main = "1atT k=4")
curve(f1, col = 1, add = TRUE)
axis(1)
y = rep(NA, N)
for(i in 1:N){
  y[i] = i * 0.04
  text(X_1atatime3[i], y[i], i, col=4)
}
points(X_1atatime3, rep(0, N), col = 2, pch = "*")
lines(X_1atatime3, y, col = 3)
hist(X_1atatime3, freq = T, main = "1atT k=4")
###############################
# Fast Algorithm, K = 20
K = 20

X_fast3 = SMED_ms_fast(mean_beta0, mean_beta1, var_mean0, var_mean1, var_e, 
                       f0 = f0, f1 = f1, type = type, var_margy0 = NULL, var_margy1 = NULL, 
                       N = N, K = K, p = 2, xmin = xmin, xmax = xmax, 
                       genCandidates = 1, initialpt = 1)

# design points locations
X_fast3_K = X_fast3$D[ , K]
curve(f0, col = 1, from = xmin, to = xmax, xlab = "", ylab = "", ylim = c(0, 3), axes = F, main = "1atT k=4")
curve(f1, col = 1, add = TRUE)
axis(1)
y = rep(NA, N)
for(i in 1:N){
  y[i] = i * 0.04
  text(X_fast3_K[i], y[i], i, col=4)
}
points(X_fast3_K, rep(0, N), col = 2, pch = "*")
lines(X_fast3_K, y, col = 3)
hist(X_fast3_K, freq = T, main = "1atT k=4")


