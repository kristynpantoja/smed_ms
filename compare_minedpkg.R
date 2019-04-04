#########################################################
# My Implementation of Fast SMED, 1-dimensional functions
#########################################################

# charge function at design point x
q = function(f, x, p = 1){ # p = 1 since it's the 1-dimensional case rn.
  1.0 / ((f(x))^(1/(2 * p)))
}

f_min_fast = function(f, candidate_jk, D_k, gamma_k){
  q(f, candidate_jk)^gamma_k * max(sapply(D_k, function(x_i) (q(f, x_i)^gamma_k / abs(x_i - candidate_jk))))
}

fast_smed_1d = function(f, N = 11, xmin = 0, xmax = 1, K, p = 1){
  
  # -- Make D_1 -- #
  # check that n >= 3
  if(N < 3) stop("not enough samples - need at least 3.")
  
  # generate candidate points, C1. for first design, C1 = D1 = lattice over [0, 1]^p
  C1 = seq(from = xmin, to = xmax, length.out = N)
  D1 = C1
  
  # -- If K = 1, return the space-filling design -- #
  if(K == 1){
    D = D1
    C = C1
    return(list("beta0" = mean_beta0, "beta1" = mean_beta1, "D" = D, "candidates" = C))
  }
  
  # -- If K > 1, choose new design -- #
  D = matrix(rep(D1, K), nrow = N, ncol = K)
  gammas = c(1:K) / (K - 1)
  # save candidates for each K
  C <- list()
  for (j in 1:N){
    C[[j]] = D1
  }
  
  # at index k, determine the next design k + 1
  for(k in 1:(K - 1)){
    
    ## For j = 1, i.e. 1st point in design k + 1:
    
    # Get candidates in neighborhood L1k = (lower, upper):
    # for j = 1
    # get candidates in L_1k
    R1k = min(abs(D[-1, k] - D[1, k])) # radius of L1k # can do it like this bc sorted D1, which was used to initialize D
    L1k_lower = max(D[1, k] - R1k, 0) # is this necessary, to max with 0?
    L1k_upper = min(D[1, k] + R1k, 1) # HERE IT IS BECAUSE o/w GET NaNs in q evaluation!
    # candidates from space-filling design, tildeD1_kplus1
    tildeD1_kplus1 = seq(from = L1k_lower, to = L1k_upper, length.out = N)
    # save the candidates to be used in future designs
    #candidates_1k = c(candidates_1k, candidates)
    C[[1]] = c(C[[1]], tildeD1_kplus1)
    # criterion to choose first candidate from candidate set: the point at which f1 and f2 are most different
    q_evals = sapply(C[[1]], FUN = function(x) q(f, x))
    xinitind = which.min(q_evals)
    #xinitind <- which.max(f(C[[1]])) # get the point that has the highest likelihood
    
    D[1, k + 1] = C[[1]][xinitind] # x1, first element of set of design points, D
    
    # for j = 2:n
    for(j in 2:N){
      # get candidates in neighborhood L_jk = (lower, upper)
      if(j == N){
        #R_jk = D[j, k] - D[j - 1, k]
        R_jk = min(abs(D[-j, k] - D[j, k]))
        lower = min(D[j, k] + R_jk, 1)
        upper = min(D[j, k] + R_jk, 1)
        tildeDj_kplus1 = seq(from = lower, to = upper, length.out = N)
        C[[j]] = c(C[[j]], tildeDj_kplus1) # This is now C_j^{k+1}
        
        
        f_min_candidates = sapply(C[[j]], function(x) f_min_fast(f, x, D[1:(j - 1), k + 1], gammas[k]))
        #choose that which has largest evaluation of criterion
        chosen_cand = which.min(f_min_candidates)
        D[j, k + 1] = C[[j]][chosen_cand]
      } else{
        R_jk = min(abs(D[-j, k] - D[j, k])) #which.min(c(D[j, k] - D[j - 1, k], D[j + 1, k] - D[j, k]))
        lower = max(D[j, k] - R_jk, 0) 
        upper = min(D[j, k] + R_jk, 1)
        tildeDjk = seq(from = lower, to = upper, length.out = N)
        C[[j]] = c(C[[j]], tildeDjk) # This is now C_j^{k+1}
        
        f_min_candidates = sapply(C[[j]], function(x) f_min_fast(f, x, D[1:(j - 1), k + 1], gammas[k]))
        #choose that which has largest evaluation of criterion
        chosen_cand = which.min(f_min_candidates)
        D[j, k + 1] = C[[j]][chosen_cand]
      }
      
    }
  }
  
  return(list("D" = D, "candidates" = C))
}






##### test, f1 is normal(1/2, 1/9) in [0, 1]

# my version
K = 20
N = 11
xmin=0; xmax=1; p = 1
f1 = function(x) exp(-(x - 1/2)^2 / (2 * (1/9)^2))
logf1 = function(x) -(x - 1/2)^2 / (2 * (1/9)^2)
xmin=0; xmax=1
f1_smedtest = fast_smed_1d(f1, N = N, xmin=xmin, xmax=xmax, K = K)
# plot
Xtest1 = f1_smedtest$D[,K]
fXtest1 = f1(Xtest1)
curve(f1,from=xmin,to=xmax)
for(i in 1:N){
  text(Xtest1[i],fXtest1[i],i,col=4)
  points(Xtest1[i],0,col=2)
}

# mined version
library(mined)
initial = matrix(seq(from = xmin, to = xmax, length.out = N), N, p)
f1_smed = mined(initial, logf1)
# plot
Xmined_f1 = f1_smed$points
fXmined_f1 = f1(Xmined_f1)
curve(unnormalized_normal,from=xmin,to=xmax)
for(i in 1:N){
  text(Xmined_f1[i],fXmined_f1[i],i,col=4)
  points(Xmined_f1[i],0,col=2)
}


##### test, f2 is beta(a = 1/2, b = 1/2)

# just to see what it looks like
a = 2
b = 5
beta_dens = function(x) dbeta(x, a, b)
curve(beta_dens)

# my version
K = 4 # for higher K, it's like spacing becomes more important.
N = 11
xmin = 0; xmax = 1; p = 1
f2 = function(x) x^(a - 1) * (1 - x)^(b - 1)
logf2 = function(x) (a - 1) * log(x) + (b - 1) * log(1 - x)
f2_smedtest = fast_smed_1d(f2, N = N, xmin=xmin, xmax=xmax, K = K)
# plot
Xtest2 = f2_smedtest$D[,4]
fXtest2 = f2(Xtest2)
curve(f2,from=xmin,to=xmax)
for(i in 1:N){
  text(Xtest2[i],fXtest2[i],i,col=4)
  points(Xtest2[i],0,col=2)
}

# mined version
initial = matrix(seq(from = xmin, to = xmax, length.out = N), N, p)
f2_smed = mined(initial, logf2)
# plot
Xmined_f2 = f2_smed$points
fXmined_f2 = f2(Xmined_f2)
curve(f2,from=xmin,to=xmax)
for(i in 1:N){
  text(Xmined_f2[i],fXmined_f2[i],i,col=4)
  points(Xmined_f2[i],0,col=2)
}

