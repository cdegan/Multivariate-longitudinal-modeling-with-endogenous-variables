#################################### FUNCTIONS TO PREPARE THE SIMULATION OF THE MODELS ###########################################

#--------------- ESTIMATION OF THE STARTING POINTS FOR THE HYPERPRIOR OF THE VARIANCE-COVARIANCE MATRIX OF THE RANDOM EFFECTS --------------------

# INPUT:
# data: dataset (data.frame)
# N: number of patients (number)
# family: distribution (e.g., gaussian, beta, gamma, etc.) of the [1] endogenous variable, [2] response variable (vector)

lmm_start <- function(data, N, family){
  N_obs <- nrow(data)
  data$time_random <- data$id+N
  
  model_x <- inla(x ~ t + f(id, model = "iid2d", n=2*N) + f(time_random, t, copy = "id"), 
                  data = data %>% filter(!is.na(x)), family = family[1], 
                  control.predictor = list(compute=TRUE, link = rep(1, nrow(data %>% filter(!is.na(x))))))
  
  model_y <- inla(y ~ t + f(id, model = "iid2d", n=2*N) + f(time_random, t, copy = "id"), 
                  data = data %>% filter(!is.na(x)), family = family[2], 
                  control.predictor = list(compute=TRUE, link = rep(1, nrow(data %>% filter(!is.na(x))))))
  
  require(brinla)
  sd_x <- bri.hyperpar.summary(model_x)[2:3, 1]
  sd_y <- bri.hyperpar.summary(model_y)[2:3, 1]
  
  return(list(sd_x = sd_x, sd_y = sd_y))
}

# OUTPUT: 
# sd_x: standard deviations of random intercept and random slopes of x (vector)
# sd_y: standard deviations of random intercept and random slopes of y (vector)



#--------------- ESTIMATION OF THE JOINT MODELS --------------------

# INPUT:
# model_list: place where to save the model results of the simulation (list)
# dataset: output of the function "single.data" in the file Data_creation.R (list)
# M: number of simulations
# default: indicating if using the default hyperprior or not for the var-cov of the random effects (boolean)
# type: name of the model we want to estimate: "JMM" or "JSM"(character)

estimation_models <- function(model_list, dataset, M, default, type){
  # initialization progress bar
  pb <- progress_bar$new(
    format = "[:bar] :percent :eta",
    total = length(1:M)
  )
  
  # inizialization estimation
  N <- length(unique(dataset[[type]][[1]]$id))
  start_iter <- length(model_list[[type]]) + 1
  
  for (i in start_iter:M) {
    # defining the starting point of the hyperprior of the variance-covariance matrix of the random effects
    if (default) {
      vars <- list("sd_x" = c(1, 1), "sd_y" = c(1, 1))
    } else {vars <- lmm_start(dataset[[type]][[i]], N, fam)}
    
    # estimating the models
    if (type == "JMM") {
      print(paste("JMM number", i))
      model_list[[type]][[i]] <- jmm(dataset[[type]][[i]], N, vars$sd_x, vars$sd_y, fam, a, time, cores, link_function)
    } else {
      print(paste("JSM number", i))
      model_list[[type]][[i]] <- jsm(dataset[[type]][[i]], N, time, vars$sd_x, vars$sd_y, fam, cores, a, link_function)
    }
    
    # saving everything
    if (default) {
      save(model_list, file = paste0("models_", type, "_beta_", N, "default.RData"))
    } else {save(model_list, file = paste0("models_", type, "_beta_", N, ".RData"))}
    
    # Update the progress bar
    pb$tick()
  }
}

# OUTPUT:
# a list of estimated models
