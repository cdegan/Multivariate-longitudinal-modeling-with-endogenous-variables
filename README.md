# Multivariate longitudinal modeling with continuous time-varying endogenous covariates - Code
Here one can find the R-code used within the context of the paper '*Multivariate longitudinal modeling with continuous time-varying endogenous covariates*'.

The folders are the following:
- **Simulation 65 vs 200** contains all the files that are needed to compile the simulation that compare the models estimation when assuming small or big sample size. This simulation is shown in Section 5 of the paper. The folder contains:
     -   *Simulated models* folder cointaining the files with the saved simulated estimated models, distringuish between JMM and JSM, and between 65 and 200 patients (*models_JMM_beta_65.RData*, *models_JMM_beta_200.RData*, *models_JSM_beta_65.RData*, *models_JSM_beta_65.RData*) (RData). The simulated data files are not available here due to their large size. If needed, you can request them by emailing c.degan@lumc.nl
     -   *Data_parameters.R* file containing all the fixed parameters used for the simulation (R file);
     -   *Data_creation.R* file containing the functions used to create the datasets for the simulation (R file);
     -   *models_preparation.R* file containing the functions used to prepare estimation of JMM and JSM (R file);
     -   *JMM.R* file containing the specific functions used to estimate the JMM (R file);
     -   *JSM.R* file containing the specific functions used to estimate the JMM (R file);
     -   *Ass_coeff_derivation.R* file containing the functions used to estimate the association coefficient in both JMM and JSM (R file);
     -   *Simulation_65vs200.Rmd* file containing the code used to implement the simulation and look at the result (R-markdown file).
       
- **Simulation default vs informative** contains all the files that are needed to compile the simulation that compare the models estimation with or without an informative hyperprior for the variance-covariance matrix of the random effects. This simulation is shown in the supplementary material. The folder contains:
     -   *Data* folder   
  
