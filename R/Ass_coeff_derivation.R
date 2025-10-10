##################### FUNCTIONS TO DERIVE THE ASSOCIATION COEFFICIENT IN JMM AND JSM #######################################

#------------------- DERIVATION OF THE ASSOCIATION COEFFICIENT IN JMM WITH THE CLOSED FORMULA (NORMAL CASE) ------------------- 

# INPUT:
# model: JMM model
# time: time points for which I want the estimation (number/vector)

coeff_fun <- function(model, time) {
  # resampling
  temp <- inla.hyperpar.sample(10^5, model)
  
  # variance of the endogenous variable
  var_x <- lapply(time, 
                  function(x) 1/temp[, 1] + 1/temp[, 3] + (x^2) * 1/temp[, 5] + 2*x*temp[, 8]*sqrt(1/temp[, 3])*sqrt(1/temp[, 5]))
  
  # covariance between response and endogenous variables
  cov_xy <- lapply(time, function(x) temp[, 7]*sqrt(1/temp[, 3])*sqrt(1/temp[, 4]) + 
                     x * (temp[, 9]*sqrt(1/temp[, 3])*sqrt(1/temp[, 6])) +
                     x * (temp[, 10]*sqrt(1/temp[, 4])*sqrt(1/temp[, 5])) +
                     (x^2) * (temp[, 12]*sqrt(1/temp[, 5])*sqrt(1/temp[, 6])))
  
  # closed formula
  coeff_sampled <- mapply(function(x, y) x/y, cov_xy, var_x)
  coeff <- cbind(time, t(apply(coeff_sampled, 2, function(x) quantile(x, probs = c(0.025, 0.5, 0.975), na.rm = TRUE))))
  
  return(coeff)
}

# OUTPUT:
# coeff: matrix of association coefficients (and its credible intervals) estimated for each time point



#------------------- DERIVATION OF THE ASSOCIATION COEFFICIENT WITH THE INTEGRATION IN JMM -------------------

# INPUT: 
# model: JMM model
# fix_eff: fixed_effects of JMM (vector)
# R: number of Montecarlo samples (number)
# n: number of hyperparameters resamples (number)
# a: value of x that we are interesting on (number)
# t: for which time points I want the estimations (number/vector)
# cores: number of cores available (number)
# link_function: name of the link function used in the model (character)

marginal_jmm <- function(model, fix_eff, R, n, a, t, cores, link_function){
  # resampling
  temp <- inla.hyperpar.sample(n, model)
  
  # marginal of the random effects b
  mean_b <- rep(0, 4)
  var_b <- list(NA)
  for(i in 1:n){
    cor <- matrix(1, nrow=4, ncol=4)
    cor[lower.tri(cor)] <- temp[i, 7:12]
    cor[upper.tri(cor)] <- t(cor)[upper.tri(cor)]
    sd <- c(sqrt(1/temp[i, 3]), sqrt(1/temp[i, 4]), sqrt(1/temp[i, 5]), sqrt(1/temp[i, 6]))
    var_b[[i]] <- diag(sd) %*% cor %*% diag(sd)
  }
  
  # marginal of x
  beta_x <- fix_eff[1:2]
  mean_x <- c(1, t) %*% beta_x
  sd_x <- apply(temp, 1, function(x) sqrt(1/x[3]+(1/x[5])*(t^2)+2*t*(x[8]*sqrt(1/x[5])*sqrt(1/x[3]))+1/x[1]))
  
  # posterior distribution of the random effects b|x
  vector_bx <- lapply(var_b, function(x) x[1, ] + t*x[, 3])
  mean_bx <- var_bx <- r_bx <- mean_bx2 <- r_bx2 <- list(NA)
  for(i in 1:n){
    mean_bx[[i]] <- mean_b + vector_bx[[i]]*sd_x[i]^(-2)*as.vector(a-mean_x)
    var_bx[[i]] <- round(var_b[[i]] - (vector_bx[[i]] * sd_x[i]^(-2)) %*% t(vector_bx[[i]]), 6)
    r_bx[[i]] <- mvtnorm::rmvnorm(R, mean_bx[[i]], var_bx[[i]], method = "svd")
    mean_bx2[[i]] <- mean_b + vector_bx[[i]]*sd_x[i]^(-2)*as.vector((a+1)-mean_x)
    r_bx2[[i]] <- mvtnorm::rmvnorm(R, mean_bx2[[i]], var_bx[[i]], method = "svd")
  }
  
  # conditional mean of y|b 
  beta_y <- fix_eff[3:4]
  nu <- mclapply(r_bx, function(x) apply(x, 1, function(x) c(1, t)%*%beta_y+x[2]+t*x[4]), mc.cores = cores)
  nu2 <- mclapply(r_bx2, function(x) apply(x, 1, function(x) c(1, t)%*%beta_y+x[2]+t*x[4]), mc.cores = cores)
  mu_cond_y <- mclapply(nu, link_function, mc.cores = cores)
  mu_cond_y2 <- mclapply(nu2, link_function, mc.cores = cores)
  
  # averaging
  average <- unlist(mclapply(mu_cond_y, mean, na.rm = TRUE, mc.cores = cores))
  average2 <- unlist(mclapply(mu_cond_y2, mean, na.rm = TRUE, mc.cores = cores))
  mu <- mean(average2 - average)
  result <- mean(mu)
  var_tot <- var(average2 - average)
  
  return(c(result = result, variance = var_tot))
}

# OUTPUT:
# table with first row the association coefficient estimations at each time point, and second row the variances



#------------------- DERIVATION OF THE ASSOCIATION COEFFICIENT IN JSM (first formulation) WITH THE CLOSED FORMULA (NORMAL CASE) ------------------- 

# INPUT:
# model: JSM model
# time: time points for which I want the estimation (number/vector)

coeff_jsm_old <- function(model, time) {
  # resampling
  temp <- inla.hyperpar.sample(10^5, model)
  
  var_x <- lapply(time, function(x) 1/temp[, 1] + 1/temp[, 5] + (x^2 * 1/temp[, 6]) + (2*x*temp[, 7]*sqrt(1/temp[, 5])*sqrt(1/temp[,6])))
  var_m <- lapply(time, function(x) 1/temp[, 5] + (x^2 * 1/temp[, 6]) + (2*x*temp[, 7]*sqrt(1/temp[, 5])*sqrt(1/temp[,6])))
  gamma <- temp[, 11]
  
  # closed formula
  coeff_sampled <- mapply(function(y, z) gamma*(y/z), var_m, var_x)
  coeff <- cbind(time, t(apply(coeff_sampled, 2, function(x) quantile(x, probs = c(0.025, 0.5, 0.975), na.rm = TRUE))))
  
  return(coeff)
}

# OUTPUT:
# coeff: matrix of association coefficients (and its credible intervals) estimated for each time point



#-------------------- DERIVATION OF THE ASSOCIATION COEFFICIENT WITH THE INTEGRATION IN JSM (first formulation) ----------------------

# INPUT:
# model: JSM model
# fix_eff: fixed_effects of JSM (vector)
# R: number of Montecarlo samples (number)
# n: number of hyperparameters resamples (number)
# a: value of x that we are interesting on (number)
# t: for which time points I want the estimations (number/vector)
# cores: number of cores available (number)
# gamma: gamma parameter of JSM (number)
# link_function: name of the link function used in the model (character)

marginal_jsm_old <- function(model, fix_eff, R, n, a, t, cores, gamma, link_function){
  temp <- inla.hyperpar.sample(n, model)
  
  # marginal of the random effects b
  mean_b <- rep(0, 4)
  var_b <- list(NA)
  for(i in 1:n){
    cov1 <- matrix(c(1/temp[i, 5], temp[i, 7]*sqrt(1/temp[i, 5]*1/temp[i, 6]), temp[i, 7]*sqrt(1/temp[i, 5]*1/temp[i, 6]), 1/temp[i, 6]), nrow=2, ncol=2)
    cov2 <- matrix(c(1/temp[i, 8], temp[i, 10]*sqrt(1/temp[i, 8]*1/temp[i, 9]), temp[i, 10]*sqrt(1/temp[i, 8]*1/temp[i, 9]), 1/temp[i, 9]), nrow=2, ncol=2)
    var_b[[i]] <- as.matrix(bdiag(cov1, cov2)) 
  }
  
  # marginal of x
  beta_x <- fix_eff[1:2]
  mean_x <- c(1, t) %*% beta_x
  sd_x <- apply(temp, 1, function(x) sqrt(1/x[5]+(1/x[6])*(t^2)+2*t*(x[7]*sqrt(1/x[5])*sqrt(1/x[6]))+1/x[1]))
  
  # posterior distribution of the random effects b|x
  vector_bx <- lapply(var_b, function(x) x[1, ] + t*x[, 2])
  mean_bx <- var_bx <- r_bx <- mean_bx2 <- r_bx2 <- list(NA)
  for(i in 1:n){
    mean_bx[[i]] <- mean_b + vector_bx[[i]]*sd_x[i]^(-2)*as.vector(a-mean_x)
    var_bx[[i]] <- as.matrix(round(var_b[[i]] - (vector_bx[[i]] * sd_x[i]^(-2)) %*% t(vector_bx[[i]]), 6))
    r_bx[[i]] <- mvtnorm::rmvnorm(R, mean_bx[[i]], var_bx[[i]], method = "svd")
    mean_bx2[[i]] <- mean_b + vector_bx[[i]]*sd_x[i]^(-2)*as.vector((a+1)-mean_x)
    r_bx2[[i]] <- mvtnorm::rmvnorm(R, mean_bx2[[i]], var_bx[[i]], method = "svd")
  }
  
  # conditional mean of y|b 
  beta_y <- fix_eff[3:4]
  nu <- mclapply(r_bx, function(x) apply(x, 1, function(x) c(1, t)%*%beta_y + x[3] + t*x[4] + gamma*(mean_x+x[1]+t*x[2])), mc.cores = cores)
  nu2 <- mclapply(r_bx2, function(x) apply(x, 1, function(x) c(1, t)%*%beta_y + x[3] + t*x[4] + gamma*(mean_x+x[1]+t*x[2])), mc.cores = cores)
  mu_cond_y <- mclapply(nu, link_function, mc.cores = cores)
  mu_cond_y2 <- mclapply(nu2, link_function, mc.cores = cores)
  
  # averaging
  average <- unlist(mclapply(mu_cond_y, mean, na.rm = TRUE, mc.cores = cores))
  average2 <- unlist(mclapply(mu_cond_y2, mean, na.rm = TRUE, mc.cores = cores))
  mu <- mean(average2 - average)
  result <- mean(mu)
  var_tot <- var(average2 - average)
  
  return(c(result = result, variance = var_tot))
}

# OUTPUT:
# table with first row the association coefficient estimations at each time point, and second row the variances


#------------------- DERIVATION OF THE ASSOCIATION COEFFICIENT IN JSM (second formulation) WITH THE CLOSED FORMULA (NORMAL CASE) ------------------- 

# INPUT:
# model: JSM model
# time: time points for which I want the estimation (number/vector)

coeff_jsm_new <- function(model, time) {
  # resampling
  temp <- inla.hyperpar.sample(10^5, model)
  
  var_x <- lapply(time, function(x) 1/temp[, 1] + 1/temp[, 4] + (x^2 * 1/temp[, 5]) + (2*x*temp[, 6]*sqrt(1/temp[, 4])*sqrt(1/temp[, 5])))
  var_m <- lapply(time, function(x) 1/temp[, 4] + (x^2 * 1/temp[, 5]) + (2*x*temp[, 6]*sqrt(1/temp[, 4])*sqrt(1/temp[, 5])))
  gamma <- temp[, 10]
  
  # closed formula
  coeff_sampled <- mapply(function(y, z) gamma*(y/z), var_m, var_x)
  coeff <- cbind(time, t(apply(coeff_sampled, 2, function(x) quantile(x, probs = c(0.025, 0.5, 0.975), na.rm = TRUE))))
  
  return(coeff)
}

# OUTPUT:
# coeff: matrix of association coefficients (and its credible intervals) estimated for each time point



#-------------------- DERIVATION OF THE ASSOCIATION COEFFICIENT WITH THE INTEGRATION IN JSM (second formulation) ----------------------

# INPUT:
# model: JSM model
# fix_eff: fixed_effects of JSM (vector)
# R: number of Montecarlo samples (number)
# n: number of hyperparameters resamples (number)
# a: value of x that we are interesting on (number)
# t: for which time points I want the estimations (number/vector)
# cores: number of cores available (number)
# gamma: gamma parameter of JSM (number)
# link_function: name of the link function used in the model (character)

marginal_jsm_new <- function(model, fix_eff, R, n, a, t, cores, gamma, link_function){
  temp <- inla.hyperpar.sample(n, model)
  
  # marginal of the random effects b
  mean_b <- rep(0, 4)
  var_b <- list(NA)
  for(i in 1:n){
    cov1 <- matrix(c(1/temp[i, 4], temp[i, 6]*sqrt(1/temp[i, 4]*1/temp[i, 5]), temp[i, 6]*sqrt(1/temp[i, 4]*1/temp[i, 5]), 1/temp[i, 5]), nrow=2, ncol=2)
    cov2 <- matrix(c(1/temp[i, 7], temp[i, 9]*sqrt(1/temp[i, 7]*1/temp[i, 8]), temp[i, 9]*sqrt(1/temp[i, 7]*1/temp[i, 8]), 1/temp[i, 8]), nrow=2, ncol=2)
    var_b[[i]] <- as.matrix(bdiag(cov1, cov2)) 
  }
  
  # marginal of x
  beta_x <- fix_eff[1:2]
  mean_x <- c(1, t) %*% beta_x
  sd_x <- apply(temp, 1, function(x) sqrt(1/x[4]+(1/x[5])*(t^2)+2*t*(x[6]*sqrt(1/x[4])*sqrt(1/x[5]))+1/x[1]))
  
  # posterior distribution of the random effects b|x
  vector_bx <- lapply(var_b, function(x) x[1, ] + t*x[, 2])
  mean_bx <- var_bx <- r_bx <- mean_bx2 <- r_bx2 <- list(NA)
  for(i in 1:n){
    mean_bx[[i]] <- mean_b + vector_bx[[i]]*sd_x[i]^(-2)*as.vector(a-mean_x)
    var_bx[[i]] <- as.matrix(round(var_b[[i]] - (vector_bx[[i]] * sd_x[i]^(-2)) %*% t(vector_bx[[i]]), 6))
    r_bx[[i]] <- mvtnorm::rmvnorm(R, mean_bx[[i]], var_bx[[i]], method = "svd")
    mean_bx2[[i]] <- mean_b + vector_bx[[i]]*sd_x[i]^(-2)*as.vector((a+1)-mean_x)
    r_bx2[[i]] <- mvtnorm::rmvnorm(R, mean_bx2[[i]], var_bx[[i]], method = "svd")
  }
  
  # conditional mean of y|b 
  beta_y <- fix_eff[3:4]
  nu <- mclapply(r_bx, function(x) apply(x, 1, function(x) c(1, t)%*%beta_y + x[3] + t*x[4] + gamma*(mean_x+x[1]+t*x[2])), mc.cores = cores)
  nu2 <- mclapply(r_bx2, function(x) apply(x, 1, function(x) c(1, t)%*%beta_y + x[3] + t*x[4] + gamma*(mean_x+x[1]+t*x[2])), mc.cores = cores)
  mu_cond_y <- mclapply(nu, link_function, mc.cores = cores)
  mu_cond_y2 <- mclapply(nu2, link_function, mc.cores = cores)
  
  # averaging
  average <- unlist(mclapply(mu_cond_y, mean, na.rm = TRUE, mc.cores = cores))
  average2 <- unlist(mclapply(mu_cond_y2, mean, na.rm = TRUE, mc.cores = cores))
  mu <- mean(average2 - average)
  result <- mean(mu)
  var_tot <- var(average2 - average)
  
  return(c(result = result, variance = var_tot))
}

# OUTPUT:
# table with first row the association coefficient estimations at each time point, and second row the variances
