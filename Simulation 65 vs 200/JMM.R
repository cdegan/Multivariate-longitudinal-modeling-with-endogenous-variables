###################### FUNCTIONS TO SIMULATE A JOINT MIXED MODEL #####################################

#-------------- ESTIMATION OF THE VARIANCE-COVARIANCE MATRIX OF THE RANDOM EFFECTS ---------------

# INPUT:
# model: INLA model
# index: indexes of the precisions and correlations in summary.hyperpar (vector)
# n1: number of random effects (number)

obtainVarCov <- function(model, index, n1){
  # initialization
  SD <- diag(n1) 
  CI <- matrix(nrow = length(index), ncol = 2) 
  cor <- matrix(1, nrow=n1, ncol=n1)  
  
  # credible intervals
  tryCatch(CI[c(1:n1),] <- t(sapply(index[1:n1], function(i) unlist(inla.zmarginal(inla.tmarginal(function(x) (1/x), model$marginals.hyperpar[[i]]), silent=TRUE)[c(3,7)]))), 
           error = function(e) print("Some variance elements CI not well defined"))     # CI for the variances
  tryCatch(CI[(n1+1):length(index),] <- unlist(model$summary.hyperpar[index[-c(1:n1)],c(3,5)]), 
           error = function(e) print("Some correlation elements CI not well defined"))  # CI for the correlation values
  
  # standard deviations
  tryCatch(SD <- sapply(index[1:n1], function(i) inla.zmarginal(inla.tmarginal(function(x) sqrt(1/x), model$marginals.hyperpar[[i]]), silent=TRUE)$mean), 
           error = function(e) print("Some SD elements not well defined"))
  
  # correlations
  cor[lower.tri(cor)] <- model$summary.hyperpar[index[-c(1:n1)],1]
  cor[upper.tri(cor)] <- t(cor)[upper.tri(cor)]
  
  # variance-covariance matrix
  answer <- diag(SD) %*% cor %*% diag(SD)
  
  return(list(Cov = answer, CI = CI))
}

# OUTPUT: 
# Cov: variance-covariance matrix of the random effects (matrix)
# CI: credible intervals of each element of the variance-covariance matrix (matrix)



#---------------- FUNCTION TO ESTIMATE THE JMM --------------------------

# INPUT:
# data: dataset (data.frame)
# N: number of patients (number)
# sd_x, sd_y: sd values used as starting points for the hyperprior of the variance-covariance matrix of the random effects (vectors)
# family: distribution (e.g., gaussian, beta, gamma, etc.) of the [1] endogenous variable, [2] response variable (vector)
# a: values of the covariate in which we estimate the association (number)
# time: time points in which I estimate the coefficient (number/vector)
# cores: number of cores available (number)
# link: link function in the model of y (character)

jmm <- function(data, N, sd_x, sd_y, family, a, time, cores, link){
  failed <<- 0 
  N_obs <- nrow(data)
  
  fixed.effects <- list(Intercept_x = c(rep(1, N_obs), rep(NA, N_obs)),
                        t_x = c(data$t, rep(NA, N_obs)),
                        Intercept_y = c(rep(NA, N_obs), rep(1, N_obs)),
                        t_y = c(rep(NA, N_obs), data$t),
                        t = c(data$t, data$t))
  
  random.effects <- list(Random_Intercept = c(data$id, data$id+N),
                         Random_Slope = c(data$id+2*N, data$id+3*N))
  
  INLA_data <- c(fixed.effects, random.effects)
  INLA_data$Y <- list(c(data$x, rep(NA, N_obs)),
                      c(rep(NA, N_obs), data$y))
  
  INLA_formula <- Y ~ -1 + Intercept_x + t_x + Intercept_y + t_y +
    f(Random_Intercept, model = "iid4d", n = 4*N, hyper = list(theta1 = list(prior = "wishart4d", param = c(10, c(sd_x[1], sd_y[1], sd_x[2], sd_y[2]), 0)))) + 
    f(Random_Slope, t, copy = "Random_Intercept")
  
  tryCatch(JMM_INLA <- inla(INLA_formula, 
                            family = family,
                            data = INLA_data, 
                            control.compute = list(cpo = TRUE, dic = TRUE, waic = TRUE, config=TRUE),
                            control.predictor = list(compute=TRUE, link = c(rep(1, N_obs), rep(2, N_obs)))),
           error = function(e) {assign('failed', 1, envir=globalenv()); failed<<-1},
           warning = function(e) {assign('failed', 1, envir=globalenv()); failed<<-1})
  if (failed==1){
    print("Returning NA's")
    return(list(fixed_eff = rep(NA, 4),
                coeff_jmm = cbind(matrix(NA, nrow = length(time), ncol=3)),
                sds = NA, 
                hyperpar_marginals = NA,
                var_hyper = NA, 
                mlik = NA, 
                overall_dic = NA, y_dic = NA, 
                overall_waic = NA, y_waic = NA, 
                overall_cpo = NA, y_cpo = NA))
  }
  
  fix_eff <- JMM_INLA$summary.fixed[, 1]
  ass_coeff <- sapply(time, function(y) marginal_jmm(JMM_INLA, fix_eff, R, n_rep, a, y, cores, link))
  if (family[2] == "gaussian") {ass_coeff_old <- coeff_fun(JMM_INLA, time)} else {ass_coeff_old <- NA}
  return(list(fix_eff = fix_eff,
              coeff_jmm = ass_coeff, 
              coeff_jmm2 = ass_coeff_old, 
              precisions = JMM_INLA$summary.hyperpar,
              hyperpar_marginals = JMM_INLA$marginals.hyperpar,
              var_hyper = obtainVarCov(JMM_INLA, c(3:12), 4)$Cov,
              mlik = JMM_INLA$mlik[1], 
              overall_dic = JMM_INLA$dic$dic, 
              y_dic = sum(JMM_INLA$dic$local.dic[(N_obs+1):(2*N_obs)], na.rm=TRUE), 
              overall_waic = JMM_INLA$waic$waic, 
              y_waic = sum(JMM_INLA$waic$local.waic[(N_obs+1):(2*N_obs)], na.rm = TRUE),
              overall_cpo = -sum(log(JMM_INLA$cpo$cpo), na.rm=TRUE), 
              y_cpo = -sum(log(JMM_INLA$cpo$cpo[(N_obs+1):(2*N_obs)]), na.rm=TRUE)))
}

# OUTPUT:
# fix_eff: fixed effects of x and y (vector)
# coeff_jmm: estimation of association coefficient (output of function 'marginal_jmm' in file 'Ass_coeff_derivation.R') 
# coeff_jmm2: second way of estimation used only if the outcome is normally distributed (output of function 'coeff_fun' in file 'Ass_coeff_derivation.R')
# precisions: posterior summary of the hyperparameters of the model (matrix)
# hyperpar_marginals: posterior marginal distributions of the hyperparameters (list of matrices)
# var_hyper: variance-covariance matrix of the random effects (matrix)
# mlik: marginal likelihood (number)
# overall_dic, y_dic: DIC on the overall model and on the y part, respectively (numbers)
# overall_waic, y_waic: WAIC on the overall model and on the y part, respectively (numbers)
# overall_cpo, y_cpo: CPO on the overall model and on the y part, respectively (numbers)