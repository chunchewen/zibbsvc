# zibbsvc
A Bayesian Zero-Inflated Spatially Varying Coefficients Model for Overdispersed Binomial Data

Authors: Chun-Che Wen, Rajib Paul, Kelly Hunt, James O'Malley, Hong Li, Elizabeth Hill, Angela Malek, Brian Neelon

This repository includes three folders: Simulation, Synthetic, and SC Info.

The Simulation folder contains the R code to simulate the Zero-Inflated Beta-Binomial Spatially Varying Coefficient (ZIBB-SVC) model (Simulation-ZIBBSVC.R), as described in the Simulation section. We export the MCMC samples for each parameter as an R data file (Simulatio-MCMCSampler.Rda) for use in generating plots. In the Make Plot.R script, we import the MCMC samples from Simulation-ZIBBSVC.Rda and generate the plots presented in the manuscript and supplementary material.

The Synthetic folder contains a synthetic version of the South Carolina (SC) prenatal data to facilitate analysis, as the actual data cannot be shared due to policies from the SC Revenue and Fiscal Affairs Office, Health and Demographics Section, and the SC Department of Health and Environmental Control. This folder includes the synthetic dataset (Synthetic-CRF-Data.Rda), the script to run the ZIBB-SVC model on the synthetic data (Synthetic-ZIBBSVC.R), and the exported MCMC samples (Synthetic-MCMCSampler.Rda) for generating plots (Make Plot.R). However, the results, such as the maps and trend plots, may differ due to the use of synthetic data.

Note: The SC_SVI.csv file contains the Social Vulnerability Index (SVI) for each county in South Carolina, and SC_adj.csv is the 46 × 46 adjacency matrix for SC. Both files are included in the SC Info folder.

# Folder

## Simulation
 - Simulation-ZIBBSVC.R: Generates simulation data using the MCMC algorithm, calculates WAIC, and creates trace plots.

 - Simulation-MCMCSampler.Rda: Contains the MCMC samples from the Simulation-ZIBBSVC.R script.

 - Make Plot.R: Generates plots for the simulation section of the manuscript.
  
## Synthetic 
 - Synthetic-CRF-Data.Rda: Synthetic version of the South Carolina (SC) prenatal data.

 - Synthetic-ZIBBSVC.R Fit ZIBB-SVC model to the synthetic dataset

 - Synthetic-MCMCSampler.Rda: Contains the MCMC samples from the Synthetic-ZIBBSVC.R script.

 - Make Plot.R: Generates plots for the analysis section of the manuscript. However, the results, such as the maps and trend plots, may differ due to the use of synthetic data.

## SC Info

  - SC_SVI.csv: Percentile-based SVI data for each SC county.

 - SC_adj.csv: South Carolina adjacency matrix.


# Parameters
  - Alpha: Fixed effect in the binary component.
    
  - Beta: Fixed effect in the Beta-Binomial (BB) component.
    
  - Rho: Correlation parameter.
    
  - phi1 - phi8: Spatial random effects (8n x 1) vector.
    
  - Sigma_phi: 8 x 8 matrix of spatial random effect variance-covariance.


