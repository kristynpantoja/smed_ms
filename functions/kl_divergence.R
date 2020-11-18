# KL between two multivariate normals (of same dimension)
KLMVN = function(
  mu1, Sigma1, 
  mu2, Sigma2
){
  if(length(mu1) != length(mu2)) stop("Error in getKLMVN : mu1 and mu2 are not same length!")
  if(length(mu1) != dim(Sigma1)[1]) stop("Error in getKLMVN : mu1 and Sigma1 dimensions don't match!")
  if(length(mu2) != dim(Sigma2)[1]) stop("Error in getKLMVN : mu1 and Sigma1 dimensions don't match!")
  d = length(mu1)
  log.det.term = log(det(Sigma2)) - log(det(Sigma1))
  trace.term = sum(diag(solve(Sigma2, Sigma1))) 
  quadratic.term = t(mu2 - mu1) %*% solve(Sigma2, mu2 - mu1)
  KL = 0.5 * (log.det.term - d + trace.term + quadratic.term)
  return(KL)
}

# KL between two univariate normals
KLN1 = function(
  mu1, sigma1, 
  mu2, sigma2
){
  KL = log(sigma2) - log(sigma1) + ((sigma1^2 + (mu1 - mu2)^2) / (2 * sigma2^2)) - 0.5
}

# KL between two normals (wrapper for univariate and multivariate)
KLN = function(
  mu1, sigma1,
  mu2, sigma2, 
  dim = 1
){
  if(dim == 1){
    KLN1(mu1, sigma1, mu2, sigma2)
  } else{
    KLMVN(mu1, sigma1, mu2, sigma2)
  }
}