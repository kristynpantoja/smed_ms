# scaled prediction variance
getSPV = function(X, N){
  XtX = t(X) %*% X
  if(is.null(dim(X)) | (dim(X)[2] == 1)){ # if X has one dimension
    SPV = sapply(X, FUN = function(x) N * t(as.matrix(x)) %*% solve(XtX, x))
  } else{ # if X has multiple dimensions (i.e. beta is vector, not scalar)
    XtX = t(X) %*% X
    SPV = apply(X, 1, FUN = function(x) N * t(as.matrix(x)) %*% solve(XtX, x))
  }
  return(SPV)
}

# variance dispersion graphs (VDGs) to evaluate the prediction variance properties for experimental designs.
# allows the practitioner to gain an impression regarding the stability of the scaled prediction variance

# the ratio of the prediction variance to the error variance
# is not a function of the error variance. This ratio, called the relative variance of prediction, depends only on
# the design and the factor setting and can be calculated before acquiring the data.
# The prediction variance profile plots the relative variance of prediction as a function of 
# each factor at fixed values of the other factors 



### --- Compute Criterion --- ###







evalD = function(D, N, p) det(crossprod(D))

evalDe = function(D, N, p) det(crossprod(D))^(1/p) / N

evalA = function(D, N, p) {
  Mi = solve(crossprod(D) / N)
  return(sum(diag(Mi))/p)
}

evalI = function(X, D, N, p) {
  Mi = solve(crossprod(D) / N)
  XMiX = rep(NA, dim(X)[1])
  for(i in 1:dim(X)[1]){
    XMiX[i] = X[i, ] %*% Mi %*% X[i, ]
  }
  return(mean(XMiX))
}

#evalGe = function(X, D, N, p){
#  Mi = solve(crossprod(D) / N)
#  XMiX = rep(NA, dim(X)[1])
#  for(i in 1:dim(X)[1]){
#    XMiX[i] = X[i, ] %*% Mi %*% X[i, ]
#  }
#  return(p / max(XMiX))
#}



### 2d
# Scaled Prediction Variance (SPV) : N V[y-hat(x_0)] / sigma^2 = N x_0' (X'X)^(-1) x_0

getSPV_2d = function(D, N, type){
  numGridPts = 10000
  gridpts_x1 = seq(from = 0, to = 1, length.out = sqrt(numGridPts))
  gridpts_x2 = seq(from = 0, to = 1, length.out = sqrt(numGridPts))
  gridpts = cbind(rep(gridpts_x1, each = sqrt(numGridPts)), 
                  rep(gridpts_x2, times = sqrt(numGridPts))) # 10k gridpts
  X = NULL
  if(type == 4){
    N = dim(D)[1]
    X = cbind(rep(1, N), D)
    gridpts = cbind(rep(1, numGridPts), gridpts)
  }
  if(type == 5){
    N = dim(D)[1]
    X = cbind(rep(1, N), D[,1], D[,1]^2, D[,2], D[,2]^2)
    gridpts = cbind(rep(1, numGridPts), gridpts[,1], gridpts[,1]^2, gridpts[,2], gridpts[,2]^2)
  }
  XtX = crossprod(X)
  spv = N * apply(gridpts, 1, function(x) x %*% solve(XtX, x))
  spv = sort(spv)
  return(spv)
}


# Covariance function
#  here, used radial aka squared exponential aka gaussian
C_fn_elementwise_sqexp2d = function(Xi,Xj, l) exp(-0.5 * sum(((Xi - Xj) / l)^2)) 
C_fn_sqexp2d = function(X, Y, l){
  lenX = dim(X)[1]
  lenY = dim(Y)[1]
  mat = matrix(rep(NA, lenX*lenY), lenX, lenY)
  for(i in 1:lenX){
    for(j in 1:lenY){
      mat[i, j] = C_fn_elementwise_sqexp2d(X[i,], Y[j,], l)
    }
  }
  return(mat)
}

#  here is product exponential
C_fn_elementwise_prodexp2d = function(Xi,Xj, theta) exp(-sum(theta * (Xi - Xj)^2)) 
C_fn_prodexp2d = function(X, Y, theta){
  lenX = dim(X)[1]
  lenY = dim(Y)[1]
  mat = matrix(rep(NA, lenX*lenY), lenX, lenY)
  for(i in 1:lenX){
    for(j in 1:lenY){
      mat[i, j] = C_fn_elementwise_prodexp(X[i,], Y[j,], l)
    }
  }
  return(mat)
}

# doesn't work... det(solve(K_X)) = 0
getSPV_GASP_2d = function(D, N, hyperparam = NULL, type, covtype = 2){
  D = unique(D)
  if(covtype == 1){ # doesn't work yet
    C_fn_elementwise = C_fn_elementwise_prodexp2d
    C_fn = C_fn_prodexp2d
  } else if(covtype == 2){
    C_fn_elementwise = C_fn_elementwise_sqexp2d
    C_fn = C_fn_sqexp2d
  }
  numGridPts = 10000
  gridpts_x1 = seq(from = 0, to = 1, length.out = sqrt(numGridPts))
  gridpts_x2 = seq(from = 0, to = 1, length.out = sqrt(numGridPts))
  gridpts = cbind(rep(gridpts_x1, each = sqrt(numGridPts)), 
                  rep(gridpts_x2, times = sqrt(numGridPts))) # 10k gridpts
  X = NULL
  if(type[2] == 4){
    N = dim(D)[1]
    X = cbind(rep(1, N), D)
    gridpts = cbind(rep(1, numGridPts), gridpts)
    if(is.null(hyperparam)) hyperparam = rep(0.1, 3)
  }
  if(type[2] == 5){
    N = dim(D)[1]
    X = cbind(rep(1, N), D[,1], D[,1]^2, D[,2], D[,2]^2)
    gridpts = cbind(rep(1, numGridPts), gridpts[,1], gridpts[,1]^2, gridpts[,2], gridpts[,2]^2)
    if(is.null(hyperparam)) hyperparam = rep(0.1, 5)
  }
  K_X = C_fn(X, X, hyperparam)
  ones = rep(1, N)
  spv = rep(NA, numGridPts)
  for(i in 1:numGridPts){
    newx = gridpts[i,]
    k_newx = apply(X, 1, function(x) C_fn_elementwise(x, newx, hyperparam))
    term2 = t(k_newx) %*% solve(K_X, k_newx)
    term3 = (1 - t(ones) %*% solve(K_X, k_newx))^2 / (t(ones) %*% solve(K_X, ones))
    spv[i] = 1 - term2 + term3
  }
  return(spv)
}




# RMSE for surface h1
RMSE = function(pred, obs){
  sqrt(mean((pred - obs)^2))
}

getRMSE_2d = function(D, N, f_ell){
  # create hypothetical response variable based on test function, here we use f_ell
  fake_y = apply(D, 1, function(x) f_ell(x))
  # analyze design using 2nd order polynomial
  fake_data = as.data.frame(cbind(D, fake_y))
  colnames(fake_data) = c("x1", "x2", "y")
  polymodel_deg1 = lm(fake_y ~ 1 + x1 + x2, data = fake_data)
  polymodel_deg2 = lm(fake_y ~ 1 + x1 + x2 + I(x1^2) + I(x2^2), data = fake_data)
  polymodel_deg3 = lm(fake_y ~ 1 + x1 + x2 + I(x1^2) + I(x2^2) + I(x1^3) + I(x2^3), data = fake_data)
  # predict responses across grid
  numGridPts = 10000
  gridpts_x1 = seq(from = 0, to = 1, length.out = sqrt(numGridPts))
  gridpts_x2 = seq(from = 0, to = 1, length.out = sqrt(numGridPts))
  gridpts = cbind(rep(gridpts_x1, each = sqrt(numGridPts)), 
                  rep(gridpts_x2, times = sqrt(numGridPts))) # 10k gridpts
  new_data = as.data.frame(gridpts)
  colnames(new_data) = c("x1", "x2")
  polymodel_deg1_predy = predict(polymodel_deg1, newdata = new_data)
  polymodel_deg2_predy = predict(polymodel_deg2, newdata = new_data)
  polymodel_deg3_predy = predict(polymodel_deg3, newdata = new_data)
  # calculate RMSE
  obs = apply(gridpts, 1, function(x) f_ell(x))
  polymodel_deg1_RMSE = RMSE(polymodel_deg1_predy, obs)
  polymodel_deg2_RMSE = RMSE(polymodel_deg2_predy, obs)
  polymodel_deg3_RMSE = RMSE(polymodel_deg3_predy, obs)
  return(list("RMSE_deg1" = polymodel_deg1_RMSE, "RMSE_deg2" = polymodel_deg2_RMSE, "RMSE_deg3" = polymodel_deg3_RMSE))
}




