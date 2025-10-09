# number of simulations
M <- 1000

# parameters for the association coefficient estimation
time <- 1:20
R <- 10^5
n_rep <- 1
a <- 9

# maximum number of repetitive observations for each patient
n <- 12

# probability to observed each variable
p1 <- 1 - 0.28    # probability to observed x 
p2 <- 1 - 0.45    # probability to observed y

# beta coefficients
true_betas <- c(12.108,    # intercept of x
                -0.166,    # slope of x 
                4.666,     # intercept y
                -0.278)    # slope y

# variances of the error terms of x
sd_x <- 0.22

# variance-covariance matrices of the random effects
D_x <- matrix(c(0.243, -0.019, -0.019, 0.004), 2, byrow=TRUE)   # var-cov matrix of the random effects of x
D_y <- matrix(c(6.004, -0.38, -0.38, 0.0416), 2, byrow=TRUE)    # var-cov matrix of the random effects of y

D <- rbind(cbind(D_x, matrix(c(0.65399, -0.056380, -0.032162, 0.005315), 2, byrow=TRUE)),
           cbind(t(matrix(c(0.65399, -0.056380, -0.032162, 0.005315), 2, byrow=TRUE)), D_y)) # joint var-cov matrix of the random effects

# scalar factor of the JSM
gamma <- 2.57 

# precision of the beta distribution Beta(a, b) 
precision <- 32.77

# association coefficients over time
coef_jmm <- as.numeric(c("0.245330", "0.257392", "0.272877", "0.262923", "0.245552", "0.212141", "0.189932", "0.158045", 
                         "0.140864", "0.129088", "0.119382", "0.114131", "0.111976", "0.10920", "0.107785", "0.107058", 
                         "0.102536", "0.101513", "0.099039", "0.095824"))
var_coeff_jmm <- as.numeric(c("0.005213", "0.005684", "0.006386", "0.007439", "0.007777", "0.007718", "0.005925", "0.004422", 
                              "0.003124", "0.002819", "0.002649", "0.002234", "0.002072", "0.00206", "0.002087", "0.002204", 
                              "0.002254", "0.002171", "0.002073", "0.002174"))
coef_jsm <- as.numeric(c("0.2147602", "0.2245170", "0.2342511", "0.243881", "0.2537917", "0.263590", "0.2740380", "0.2850331", 
                         "0.2938918", "0.3002987", "0.3052803", "0.3057826", "0.3033340", "0.2999659", "0.2911779", "0.2825748", 
                         "0.274030", "0.263235", "0.250314", "0.240403"))
var_coef_jsm <- as.numeric(c("0.0005669", "0.0005247", "0.0004372", "0.000448", "0.0004046", "0.000399", "0.0003945", "0.0004129", 
                             "0.0004192", "0.0004757", "0.0005412", "0.0006155", "0.0007095", "0.0008437", "0.0008766", "0.0009753", 
                             "0.001036", "0.001099", "0.001257", "0.001236"))