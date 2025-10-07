# Multivariate longitudinal modeling with continuous time-varying endogenous covariates
Here one can find the R-code used within the context of the paper '*Multivariate longitudinal modeling with continuous time-varying endogenous covariates*'.

The folders are the following:
- **Simulation 65 vs 200** contains all the files used for the simulation shown in Section 5 of the paper. It contains:
     -   *Data* folder
     -   *Data_parameters.R* R file containing all the fixed parameters used for the simulation;
     -   *Data_creation.R* R file containing the functions used to create the datasets for the simulation;
     -   *models_preparation.R* R file containing the functions used to prepare estimation of JMM and JSM;
     -   *JMM.R* R file containing the specific functions used to estimate the JMM;
     -   *JSM.R* R file containing the specific functions used to estimate the JMM;
     -   *Ass_coeff_derivation.R* R file containing the functions used to estimate the association coefficient in both JMM and JSM;
     -   *Simulation_65vs200.Rmd* R-markdown file containing the code used to implement the simulation and look at the result.
       
- **Simulation default vs informative** contains all the files that are needed to compile the simulation that compare the models estimation with or without an informative hyperprior for the variance-covariance matrix of the random effects. This simulation is shown in the supplementary material. The folder contains:
     -   *Data* folder   
  
