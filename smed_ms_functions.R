# --- Functions Used in both One-At-a-Time and Fast Algorithms --- #

# Wasserstein distance betwen two (univariate) normals, N(mu1, var1) and N(mu2, var2)
Wasserstein_distance = function(mu1, mu2, var1, var2){
  return(sqrt(mu1 - mu2)^2 + var1 + var2 - 2 * sqrt(var1 * var2))
}

var_marginaly = function(x, var_e, var_mean) var_e + x^2 * var_mean

# charge function at design point x
q = function(x, mean_beta0, mean_beta1, var_e, var_mean){
  mu1 = mean_beta0 * x # mean of marginal dist of y | H0
  mu2 = mean_beta1 * x # mean of marginal dist of y | H1
  var = var_marginaly(x, var_e, var_mean) # variance of marginal dist of y | H1
  Wass_dist = Wasserstein_distance(mu1, mu2, var, var)
  return(1.0 / Wass_dist^(1/2))
}