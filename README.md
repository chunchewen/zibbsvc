# zibbsvc
A Bayesian Zero-Inflated Spatially Varying Coefficients Model for Overdispersed Binomial Data

Authors: Chun-Che Wen, Rajib Paul, Kelly Hunt, James O'Malley, Hong Li, Elizabeth Hill, Angela Malek, Brian Neelon

This repository includes R code to simulate the Zero-Inflated Beta-Binomial Spatially Varying Coefficient (ZIBB-SVC) model (Simulation-ZIBBSVC.R). We export the MCMC samples for each parameter as an R data file, called Simulation-ZIBBSVC.Rda, for use in generating plots. In the Make Plot.R script, we import the MCMC samples from Simulation-ZIBBSVC.Rda and generate the plots demonstrated in the manuscript and supplementary material. Note: The SC_SVI.csv file contains the Social Vulnerability Index (SVI) for each county in South Carolina (SC), and SC_adj.csv is the adjacency matrix for SC.

# Files
 - Simulation-ZIBBSVC.R: Generates simulation data using the MCMC algorithm, calculates WAIC, and creates trace plots.

 - Make Plot.R: Creates plots for the manuscript.

 - SC_SVI.csv: Percentile-based SVI data for each SC county.

 - SC_adj.csv: South Carolina adjacency matrix.

 - Simulation-ZIBBSVC.Rda: Contains the MCMC samples from the Simulation-ZIBBSVC.R script.

# Parameters
  - Alpha: Fixed effect in the binary component.
    
  - Beta: Fixed effect in the Beta-Binomial (BB) component.
    
  - Rho: Correlation parameter.
    
  - phi1 - phi8: Spatial random effects (8n x 1) vector.
    
  - Sigma_phi: 8 x 8 matrix of spatial random effect variance-covariance.


