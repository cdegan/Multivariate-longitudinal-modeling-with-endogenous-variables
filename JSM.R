###################### FUNCTIONS TO SIMULATE A JOINT SCALED MODEL #####################################

#-------------- ESTIMATION OF THE VARIANCE-COVARIANCE MATRIX OF THE RANDOM EFFECTS ---------------

# INPUT:
# model: INLA model
# index1: indexes of the precisions and correlations in summary.hyperpar of the random effects of x (vector)
# index2: indexes of the precisions and correlations in summary.hyperpar of the random effects of y (vector)
# n1: number of random effects of x (number)
# n2: number of random effects of y (number)

obtainVarCov_pairwise <- function(model, index1, index2, n1, n2){
  # variance-covariance matrix
  cov1 <- obtainVarCov(model, index1, n1)$Cov
  cov2 <- obtainVarCov(model, index2, n2)$Cov
  VarCov <- bdiag(cov1, cov2)
  
  # credible intervals
  tryCatch(CI <- rbind(obtainVarCov(model, index1, n1)$CI[1:n1, ], obtainVarCov(model, index2, n2)$CI[1:n1, ], # CI of variances
                       obtainVarCov(model, index1, n1)$CI[(n1+1):(length(index1)), ], obtainVarCov(model, index2, n2)$CI[(n2+1):(length(index2)), ]), # CI of correlations
           error = function(e) print("Some covariance elements CI not well defined"))
  
  return(list(Cov = VarCov, CI = CI))
}

# OUTPUT: 
# Cov: variance-covariance matrix of the random effects (matrix)
# CI: credible intervals of each element of the variance-covariance matrix (matrix)



#---------------- FUNCTION TO ESTIMATE THE JSM (first formulation) --------------------------

# INPUT:
# data: dataset (data.frame)
# N: number of patients (number)
# time: time points in which I estimate the coefficient (number/vector)
# sd_x, sd_y: sd values used as starting points for the hyperprior of the variance-covariance matrix of the random effects (vectors)
# family: distribution (e.g., gaussian, beta, gamma, etc.) of the [1] endogenous variable, [2] response variable (vector)
# cores: number of cores available (number)
# a: values of the covariate in which we estimate the association (number)
# link: link function in the model of y (character)

jsm <- function(data, N, time, sd_x, sd_y, family, cores, a, link){
  failed <<- 0
  fix_eff <- c(NA, 4)
  gamma_model <- NA
  ass_coeff <- NA
  N_obs <- nrow(data)
  
  fixed.effects <- list(Intercept_y = c(rep(NA, N_obs), rep(1, N_obs)),
                        t_y = c(rep(NA, N_obs), data$t),
                        t = c(data$t, data$t))
  
  random.effects <- list(Intercept_x = c(rep(1, N_obs), rep(NA, N_obs)),
                         t_x = c(rep(1, N_obs), rep(NA, N_obs)),
                         Intercept_x_scaled = c(rep(NA, N_obs), rep(1, N_obs)),
                         t_x_scaled = c(rep(NA, N_obs), rep(1, N_obs)), 
                         Random_intercept_x = c(data$id, rep(NA, N_obs)),
                         Random_slope_x = c(data$id+N, rep(NA, N_obs)),
                         Random_intercept_x_scaled = c(rep(NA, N_obs), data$id),
                         Random_slope_x_scaled = c(rep(NA, N_obs), data$id+N),
                         Random_intercept_y = c(rep(NA, N_obs), data$id),
                         Random_slope_y = c(rep(NA, N_obs), data$id+N))
  
  INLA_data <- c(fixed.effects, random.effects)
  INLA_data$Y <- list(c(data$x, rep(NA, N_obs)),
                      c(rep(NA, N_obs), data$y))
  
  INLA_formula = Y ~ -1 + 
    # fixed effects of x 
    f(Intercept_x) + f(t_x, t) +
    # scaled fixed effects of x 
    f(Intercept_x_scaled, copy = "Intercept_x", hyper = list(beta = list(fixed=FALSE, prior = "normal", param = c(0, 0.0001)))) +
    f(t_x_scaled, t, copy = "t_x", same.as = 'Intercept_x_scaled', hyper = list(beta = list(fixed=FALSE, prior = "normal", param = c(0, 0.0001)))) +
    # fixed effects of y
    Intercept_y + t_y + 
    # random effects for x
    f(Random_intercept_x, model = "iid2d", n = 2*N, hyper = list(theta1 = list(prior = "wishart2d", param =  c(3, c(sd_x[1], sd_x[2]), 0)))) + f(Random_slope_x, t, copy = "Random_intercept_x") +
    # random effects of x_scaled
    f(Random_intercept_x_scaled, copy = 'Random_intercept_x', same.as = 'Intercept_x_scaled', hyper = list(beta = list(fixed=FALSE, prior = "normal", param = c(0, 0.0001)))) +
    f(Random_slope_x_scaled, t, copy = 'Random_intercept_x', same.as = 'Intercept_x_scaled', hyper = list(beta = list(fixed=FALSE, prior = "normal", param = c(0, 0.0001)))) +
    # random effects of y
    f(Random_intercept_y, model = "iid2d", n = 2*N, hyper = list(theta1 = list(prior = "wishart2d", param =  c(3, c(sd_y[1], sd_y[2]), 0)))) + f(Random_slope_y, t, copy = "Random_intercept_y")
  
  tryCatch(JSM_INLA <- inla(INLA_formula, 
                            family = family,
                            data = INLA_data, 
                            control.compute = list(cpo = TRUE, dic = TRUE, waic = TRUE, config=TRUE), 
                            control.predictor = list(compute=TRUE, link = c(rep(1, N_obs), rep(2, N_obs))),
                            control.inla = list(control.vb = list(emergency = 30))),
           error = function(e) {assign('failed',1,envir=globalenv()); failed<<-1},
           warning = function(e) {assign('failed',1,envir=globalenv()); failed<<-1})
  if (failed==1){
    print("Returning NA's")
    return(list(fix_eff = NA,
                random = NA,
                gamma = NA,
                precisions = NA,
                hyperpar_marginals = NA, 
                coeff_jsm = cbind(matrix(NA, ncol = length(time), nrow=2)), 
                var_hyper = NA, 
                mlik = NA, 
                overall_dic = NA, y_dic = NA, 
                overall_waic = NA, y_waic = NA, 
                overall_cpo = NA, y_cpo = NA))
  }
  
  fix_eff <- c(JSM_INLA$summary.random$Intercept_x[, 2], JSM_INLA$summary.random$t_x[, 2], JSM_INLA$summary.fixed$mean)
  gamma_model <- JSM_INLA$summary.hyperpar[11, 1]
  ass_coeff <- sapply(time, function(y) marginal_jsm_old(JSM_INLA, fix_eff, R, n_rep, a, y, cores, gamma_model, link))
  if (family[2] == "gaussian") {ass_coeff_old <- coeff_jsm_old(JSM_INLA, time)} else {ass_coeff_old <- NA}
  return(list(fix_eff = fix_eff,
              random = rbind("Intercept_x_scaled" = JSM_INLA$summary.random$Intercept_x_scaled[, 2], "t_x_scaled" = JSM_INLA$summary.random$t_x_scaled[, 2]),
              gamma = gamma_model, 
              precisions = JSM_INLA$summary.hyperpar,
              hyperpar_marginals = JSM_INLA$marginals.hyperpar,
              coeff_jsm = ass_coeff, 
              coeff_jsm2 = ass_coeff_old, 
              var_hyper = obtainVarCov_pairwise(JSM_INLA, c(5:7), c(8:10), 2, 2)$Cov,
              mlik = JSM_INLA$mlik[1], 
              overall_dic = JSM_INLA$dic$dic,
              y_dic = sum(JSM_INLA$dic$local.dic[(N_obs+1):(2*N_obs)], na.rm=TRUE), 
              overall_waic = JSM_INLA$waic$waic, 
              y_waic = sum(JSM_INLA$waic$local.waic[(N_obs+1):(2*N_obs)], na.rm = TRUE),
              overall_cpo = -sum(log(JSM_INLA$cpo$cpo), na.rm=TRUE), 
              y_cpo = -sum(log(JSM_INLA$cpo$cpo[(N_obs+1):(2*N_obs)]), na.rm=TRUE)))
}

# OUTPUT:
# fix_eff: fixed effects of x and y (vector)
# random: random effects 
# gamma: posterior mode of the gamma hyperparameter
# precisions: posterior summary of the hyperparameters of the model (matrix)
# hyperpar_marginals: posterior marginal distributions of the hyperparameters (list of matrices)
# coeff_jsm: estimation of association coefficient (output of function 'marginal_jmm_old' in file 'Ass_coeff_derivation.R') 
# coeff_jsm2: second way of estimation used only if the outcome is normally distributed (output of function 'coeff_jsm_old' in file 'Ass_coeff_derivation.R')
# var_hyper: variance-covariance matrix of the random effects (matrix)
# mlik: marginal likelihood (number)
# overall_dic, y_dic: DIC on the overall model and on the y part, respectively (numbers)
# overall_waic, y_waic: WAIC on the overall model and on the y part, respectively (numbers)
# overall_cpo, y_cpo: CPO on the overall model and on the y part, respectively (numbers)



#---------------- FUNCTION TO ESTIMATE THE JSM (second formulation) --------------------------

# INPUT:
# data: dataset (data.frame)
# x: name of the endogenous variable in the dataset (character)
# y: name of the response variable in the dataset (character)
# id: name of the variable with subject-specific id in the dataset (character)
# t: name of the time variable in the dataset (character)
# family: distribution (e.g., gaussian, beta, gamma, etc.) of the [1] endogenous variable, [2] response variable (vector)
# sd_x, sd_y: sd values used as starting points for the hyperprior of the variance-covariance matrix of the random effects (vectors)
# a: values of the covariate in which we estimate the association (number)
# cores: number of cores available (number)
# time: time points in which I estimate the coefficient (number/vector)

jsm_new <- function(data, x, y, id, t, family, sd_x, sd_y, a, cores, time){
  failed <<- 0
  data.y <- data %>% dplyr::filter(!is.na(y))
  data.x <- data %>% dplyr::filter(!is.na(x))
  
  x.new <- x[!is.na(x)]
  y.new <- y[!is.na(y)]
  
  # number of observations for each model
  N.y <-  nrow(data.y) 
  N.x <-  nrow(data.x)
  
  # number of patients for each model
  n.y <- length(unique(pull(data.y, id)))
  n.x <- length(unique(pull(data.x, id)))
  
  # joint response variable
  Y.joint <- list(
    c(x.new, rep(NA, 2*N.y)),
    c(rep(NA, N.x), rep(c(0, NA), each = N.y)),
    c(rep(NA, N.x), rep(NA, N.y), y.new))
  
  # indices for the random effects
  # between x and y there are no differences because I have to have the same patients
  idx.int.x <- pull(data.x, id)
  idx.slo.x <- idx.int.x + n.x
  idx.int.y <- pull(data.y, id)
  idx.slo.y <- idx.int.y + n.y
  
  # linear predictors
  covariates <- list(
    # x and pseudo-submodel
    # fixed effects
    Int.x = c(rep(1, N.x + N.y), rep(NA, N.y)),
    t.x = c(pull(data.x, t), pull(data.y, t), rep(NA, N.y)),
    # random effects
    id.random.int.x = c(idx.int.x, idx.int.y, rep(NA, N.y)),
    id.random.slo.x = c(idx.slo.x, idx.slo.y, rep(NA, N.y)),
    # copied terms
    idx.copy = c(rep(NA, N.x), 1:N.y, rep(NA, N.y)),
    w.copy = c(rep(NA, N.x), rep(-1, N.y), rep(NA, N.y)),
    
    # y
    # fixed effects
    Int.y = c(rep(NA, N.x), rep(c(NA, 1), each=N.y)),
    t.y = c(rep(NA, N.x + N.y), pull(data.y, t)),
    pred.x = c(rep(NA, N.x + N.y), 1:N.y),
    # random effects
    id.random.int.y = c(rep(NA, N.x + N.y), idx.int.y),
    id.random.slo.y = c(rep(NA, N.x + N.y), idx.slo.y)
  )
  
  # joint data
  joint.data <- c(covariates)
  joint.data$Y <- Y.joint
  
  # Joint Model
  tryCatch(JSM_INLA <- inla(Y ~ -1 +
                              # x and pseudo-submodel
                              Int.x + t.x + 
                              f(id.random.int.x, Int.x, model = "iid2d", n = 2*n.x, hyper = list(theta1 = list(prior = "wishart2d", param =  c(3, c(sd_x[1], sd_x[2]), 0)))) + f(id.random.slo.x, t.x, copy = "id.random.int.x") + # linear predictor of x
                              f(idx.copy, w.copy, model = "iid", hyper = list(prec = list(initial = -6, fixed = TRUE))) + # u part in the pseudo-submodel definition
                              # y
                              Int.y + t.y +
                              f(id.random.int.y, model = "iid2d", n = 2*n.y, hyper = list(theta1 = list(prior = "wishart2d", param =  c(3, c(sd_y[1], sd_y[2]), 0)))) + f(id.random.slo.y, t.y, copy = "id.random.int.y") + # linear predictor of y
                              f(pred.x, copy = "idx.copy", hyper = list(beta = list(fixed=FALSE, param = c(0, 0.0001)))), # copy and scale the linear predictor of x
                            data = joint.data,  
                            family = family, 
                            control.predictor = list(compute=TRUE, link = c(rep(1, N.x), rep(2, N.y), rep(3, N.y))),
                            control.compute = list(cpo = TRUE, dic = TRUE, waic = TRUE, config = TRUE),
                            control.inla = list(control.vb = list(emergency = 30))),
           error = function(e) {assign('failed',1,envir=globalenv()); failed<<-1},
           warning = function(e) {assign('failed',1,envir=globalenv()); failed<<-1})
  if (failed==1){
    print("Returning NA's")
    return(list(fixed = NA,
                gamma = NA,
                precisions = NA,
                hyperpar_marginals = NA,
                coeff_jsm = cbind(matrix(NA, ncol = length(time), nrow=2)),
                var_hyper = NA,
                mlik = NA,
                overall_dic = NA, y_dic = NA,
                overall_waic = NA, y_waic = NA,
                overall_cpo = NA, y_cpo = NA))
  }
  
  fix_eff <- JSM_INLA$summary.fixed$mean
  gamma_model <- JSM_INLA$summary.hyperpar$mean[which(rownames(JSM_INLA$summary.hyperpar) == "Beta for pred.x")]
  ass_coeff <- sapply(time, function(y) marginal_jsm_new(JSM_INLA, fix_eff, R, n_rep, a, y, cores, gamma_model, link))
  return(list(fixed = fix_eff,
              gamma = gamma_model,
              precisions = JSM_INLA$summary.hyperpar,
              hyperpar_marginals = JSM_INLA$marginals.hyperpar,
              coeff_jsm = ass_coeff,
              var_hyper = obtainVarCov_pairwise(JSM_INLA, c(4:6), c(7:9), 2, 2)$Cov,
              mlik = JSM_INLA$mlik[1],
              overall_dic = JSM_INLA$dic$dic,
              overall_waic = JSM_INLA$waic$waic,
              overall_cpo = -sum(log(JSM_INLA$cpo$cpo), na.rm=TRUE)))
}

# OUTPUT:
# fix_eff: fixed effects of x and y (vector)
# gamma: posterior mode of the gamma hyperparameter
# precisions: posterior summary of the hyperparameters of the model (matrix)
# hyperpar_marginals: posterior marginal distributions of the hyperparameters (list of matrices)
# coeff_jsm: estimation of association coefficient (output of function 'marginal_jmm_old' in file 'Ass_coeff_derivation.R') 
# var_hyper: variance-covariance matrix of the random effects (matrix)
# mlik: marginal likelihood (number)
# overall_dic: DIC on the overall model (numbers)
# overall_waic: WAIC on the overall model (numbers)
# overall_cpo: CPO on the overall model (numbers)
