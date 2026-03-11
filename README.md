# BayesITR
R data and code for the paper: "Robust Bayesian Learning for Individualized Treatment Rules Under Unmeasured Confounding".

## Data

The data analyzed in Section 6 "Application: DIVAT Data Analysis" of the paper are from The French Computerized and Validated Data in Transplantation (DIVAT), which is a large-scale database that collects medical records for kidney transplantation at several French hospitals. The data are publicly available. Full details of the data are available at https://www.divat.fr/. R data for the simulation study and the DIVAT data analysis of the paper are available at the following link: https://drive.google.com/file/d/1ctH0AXDrb596skbTCQy9F0mbGGJbJz0L/view?usp=sharing.

## Code 

The R scripts in the folders "Bayes_ITR_Simulation_Scenario_I" and "Bayes_ITR_Simulation_Scenario_II" are for Section 5 "Simulation Study", the R scripts in the folder "Bayes_ITR_DIVAT" are for Section 6 "Application: DIVAT Data Analysis", the R scripts in the folders "Bayes_ITR_Supplement_C", "Bayes_ITR_Supplement_E1", "Bayes_ITR_Supplement_E2", "Bayes_ITR_Supplement_E3", and "Bayes_ITR_Supplement_E4" are for the additional simulation studies reported in Supplementary Material Sections C and E1–E4.

Libraries and Version Numbers: R version 4.4.2, Rcpp 1.0.14, RcppArmadillo 14.2.2-1, Matrix 1.7-1, MCMCpack 1.7-1, LaplacesDemon 16.1.6, splines2 0.5.3, truncnorm 1.0-9, loo 2.8.0, twangContinuous 1.0.0, ggplot2 3.5.1, gridExtra 2.3.

### Instructions for Use

In the folder "Bayes_ITR_Simulation_Scenario_I": 

* The R script "Simu_Main.R" reproduces Simulation Study Scenario I in Section 5, including Figure 3 in the manuscript and Figures S1, S3, and S4 in the supplementary material.
