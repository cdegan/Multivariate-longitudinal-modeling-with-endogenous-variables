######################### FUNCTIONS TO CREATE DATASETS FROM JMM AND JSM WITH ONE GAUSSIAN ENDOGENOUS VARIABLE AND ONE BETA RESPONSE VARIABLE #############################

# ------------------ DEFINITION OF THE DATA STRUCTURE -------------------- 

# INPUT:
# N: number of patient (number)
# n: number of repeated measurements per patient (number)
# p1: percentages of observed (not-missing) values in the x (number)
# p2: percentages of observed (not-missing) values in the y (number)

create_data <- function(N, n, p1, p2){
  t_list <- list()
  for(i in 1:N) t_list[[i]] <- sort(sample(seq(3, 26, by = 0.10), n))
  
  data <- data.frame(
    id = rep(1:N, each = n),
    t = unlist(t_list), 
    observed_x = rbern(N*n, p1),
    observed_y = rbern(N*n, p2),
    x = NA, 
    y = NA
  )
  
  return(data)
}

# OUTPUT:
# data: dataset structure



#-------------------- SIMULATION FROM A JMM --------------------- 

# INPUT:
# data: output of function 'create_data' 
# N: number of patients
# true_betas: values that we assume to be the real estimation coefficients (vector)
# D: variance-covariance matrix of the random effects involved in JMM, that we assume to be the real one (matrix)
# precision: precision of the beta distribution of y that we assume to be the real one (number)
# sd_x: standard deviation of the error terms that we assume to be the real one (number)

create_data_jmm <- function(data, N, true_betas, D, precision, sd_x){
  # error terms 
  errorX <- rnorm(N*n, sd = sd_x)
  
  # correlated 4 random effects: first 2 are the random effects of x, second 2 are the random effects of y
  u <- mvrnorm(N, mu = rep(0, 4), Sigma = D)
  
  # variable x
  data$x <- true_betas[1] + data$t * true_betas[2] + rep(u[, 1], table(data$id)) + 
    data$t * rep(u[, 2], table(data$id)) + errorX
  
  # variable y
  eta <- true_betas[3] + data$t * true_betas[4] + rep(u[, 3], table(data$id)) + 
    data$t * rep(u[, 4], table(data$id))                    # linear predictor
  mean <- plogis(eta)                                       # logit function of the linear predictor
  a <- mean*precision                                       # shape parameter number 1
  b <- -mean*precision + precision                          # shape parameter number 2
  for(i in 1:(N*n)) data$y[i] <- rbeta(1, a[i], b[i])       # y follows a beta distribution
  data$y <- ifelse(data$y == 0, data$y + 0.0001, data$y)
  data$y <- ifelse(data$y == 1, data$y - 0.0001, data$y)
  
  # randomly create an unbalanced data setting
  data$x <- ifelse(data$observed_x == 0, NA, data$x)
  data$y <- ifelse(data$observed_y == 0, NA, data$y)
  
  return(data)
}

# OUTPUT:
# data: completed dataset created based on JMM 



#-------------------- SIMULATION FROM A JSM --------------------- 

# INPUT:
# data: output of function 'create_data' 
# N: number of patients
# true_betas: values that we assume to be the real estimation coefficients (vector)
# D_x: variance-covariance matrix of the random effects of x, that we assume to be the real one (matrix)
# D_y: variance-covariance matrix of the random effects of y, that we assume to be the real one (matrix)
# gamma: gamma parameter of the JSM that we assume to be the real one (number)
# precision: precision of the beta distribution of y that we assume to be the real one (number)
# sd_x: standard deviation of the error terms that we assume to be the real one (number)

create_data_jsm <- function(data, N, true_betas, D_x, D_y, gamma, sd_x, precision){
  # error terms
  errorX = rnorm(N*n, sd = sd_x)
  
  # reformulation of the true beta of y
  true_betas[3:4] <- true_betas[3:4] - gamma*true_betas[1:2]
  
  # random effects
  u_x <- mvrnorm(N, mu = rep(0, 2), Sigma = D_x)
  u_y <- mvrnorm(N, mu = rep(0, 2), Sigma = D_y)
  
  # linear predictor m
  m <- true_betas[1] + data$t * true_betas[2] + rep(u_x[, 1], table(data$id)) + 
    data$t * rep(u_x[, 2], table(data$id)) 
  
  # variable x
  data$x <- m + errorX
  
  # variable y
  eta <- true_betas[3] + data$t * true_betas[4] + rep(u_y[, 1], table(data$id)) + 
    data$t * rep(u_y[, 2], table(data$id)) + gamma*m      # linear predictor
  mean <- plogis(eta)                                     # logit function of the linear predictor
  a <- mean*precision                                     # shape parameter number 1
  b <- -mean*precision + precision                        # shape parameter number 2
  for(i in 1:(N*n)) data$y[i] <- rbeta(1, a[i], b[i])     # y follows a beta distribution 
  data$y <- ifelse(data$y == 0, data$y + 0.0001, data$y)
  data$y <- ifelse(data$y == 1, data$y - 0.0001, data$y)
  
  # randomly create an unbalanced data setting
  data$x <- ifelse(data$observed_x == 0, NA, data$x)
  data$y <- ifelse(data$observed_y == 0, NA, data$y)
  
  return(data)
}

# OUTPUT:
# data: completed dataset created based on JMM 



# ------------------ CREATION OF THE DATA -------------------- 

# INPUT:
# N: number of patient (number)
# M: number of simulations (number)

single.data <- function(N, M){
  data <- list()
  for(i in 1:M) data[[i]] <- create_data(N, n, p1, p2)
  
  data_jmm <- lapply(data, function(x) create_data_jmm(x, N, true_betas, D, precision, sd_x))
  data_jsm <- lapply(data, function(x) create_data_jsm(x, N, true_betas, D_x, D_y, gamma, sd_x, precision))
  datasets <- list("JMM" = data_jmm, "JSM" = data_jsm)
  
  return(datasets)
}

# OUTPUT:
# datasets: list of JMM and JSM datasets, [[1]] JMM data and [[2]] JSM data
