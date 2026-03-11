# BayesITR
R data and code for the paper: "Robust Bayesian Learning for Individualized Treatment Rules Under Unmeasured Confounding".

## Data

The data analyzed in Section 6 "Application: DIVAT Data Analysis" of the paper are from The French Computerized and Validated Data in Transplantation (DIVAT), which is a large-scale database that collects medical records for kidney transplantation at several French hospitals. The data are publicly available. However, one needs to fill out a request form for access. Full details of the data are available at https://www.divat.fr/. R data for the simulation study and the DIVAT data analysis of the paper are available at: https://drive.google.com/file/d/1ctH0AXDrb596skbTCQy9F0mbGGJbJz0L/view?usp=sharing.

## Code 

The R scripts in the folders "Bayes_ITR_Simulation_Scenario_I" and "Bayes_ITR_Simulation_Scenario_II" are for Section 5 "Simulation Study", the R scripts in the folder "Bayes_ITR_DIVAT" are for Section 6 "Application: DIVAT Data Analysis", the R scripts in the folders "Bayes_ITR_Supplement_C", "Bayes_ITR_Supplement_E1", "Bayes_ITR_Supplement_E2", "Bayes_ITR_Supplement_E3", and "Bayes_ITR_Supplement_E4" are for the additional simulation studies reported in Supplementary Material Sections C and E1–E4.

Libraries and Version Numbers: R version 4.4.2, Rcpp 1.0.14, RcppArmadillo 14.2.2-1, Matrix 1.7-1, MCMCpack 1.7-1, LaplacesDemon 16.1.6, splines2 0.5.3, truncnorm 1.0-9, loo 2.8.0, twangContinuous 1.0.0, ggplot2 3.5.1, gridExtra 2.3.

### Instructions for Use

**In the folder "Bayes_ITR_Simulation_Scenario_I":**

* The R script "Simu_Main.R" reproduces Simulation Study Scenario I in Section 5, including Figure 3 in the manuscript and Figures S1, S3, and S4 in the supplementary material.

* The R scripts "Simu_Data_Generate_I(a).R" and "Simu_Data_Generate_I(b).R" generate the simulated datasets for Scenarios I(a) and I(b), respectively. "MCMC_R_Functions.R" provides the R functions used for MCMC, "MCMC_Rcpp_Functions.cpp" provides the corresponding Rcpp functions, and "MCMC_Main.R" contains the main MCMC function.

* The R data files "MCMC_Results_I(a)_Mode1.RData" and "MCMC_Results_I(a)_Mode2.RData" save the MCMC posterior samples for the two modes under Scenario I(a), while "MCMC_Results_I(b)_Mode1.RData" and "MCMC_Results_I(b)_Mode2.RData" save the MCMC posterior samples for the two modes under Scenario I(b).

**In the folder "Bayes_ITR_Simulation_Scenario_II":**

* The R script “Simu_Main.R” reproduces Simulation Study Scenario II in Section 5, including Figure 4 in the manuscript and Figure S5 in the supplementary material.

* The R script "Simu_Data_Generate.R" generates the simulated datasets for Scenario II, and "Select_Partial_Data.R" selects partial data to create the sub-scenarios II(a), II(b), and II(c). "MCMC_R_Functions.R" provides the R functions used for MCMC, "MCMC_Rcpp_Functions.cpp" provides the corresponding Rcpp functions, and "MCMC_Main.R" contains the main MCMC function. "SGD_Functions.R" provides the R functions used for SGD, while "SGD_Main_Conservative.R" and "SGD_Main_Expected.R" contain the main SGD functions for the proposed conservative policy optimization and the alternative non-conservative policy optimization, respectively.

* The R data files "SGD_Results_II(a).RData", "SGD_Results_II(b).RData", and "SGD_Results_II(c).RData" save the SGD results for Scenarios II(a), II(b), and II(c), respectively.

**In the folder "Bayes_ITR_DIVAT":**

* The R script "DIVAT_Main.R" reproduces the DIVAT data analysis in Section 6, including Figure 5 and Table 1 in the manuscript and Figure S6 in the supplementary material.

* The R script "Data_Process.R" processes the DIVAT data for analysis. "MCMC_R_Functions.R" provides the R functions used for MCMC, "MCMC_Rcpp_Functions.cpp" provides the corresponding Rcpp functions, and "MCMC_Main.R" contains the main MCMC function.

* The R data file "DIVAT_Data.RData" saves the DIVAT dataset for analysis, while "MCMC_Results.RData" and "MCMC_Results_NoUC.RData" save the MCMC posterior samples for the proposed method and the alternative "NoUC" method, respectively.

**In the folder "Bayes_ITR_Supplement_C":**

* The R script "Supplement_C_Main.R" reproduces the additional simulation studies in Supplementary Material Section C, including Figure C1 and Figure S2 in the supplementary material.

* The R script "Simu_Data_Generate.R" generates the simulated datasets for Supplementary Material Section C. "MCMC_R_Functions.R" provides the R functions used for MCMC, "MCMC_Rcpp_Functions.cpp" provides the corresponding Rcpp functions, and "MCMC_Main.R" contains the main MCMC function.

* The R data files "MCMC_Results_B0=3.RData" through "MCMC_Results_B0=10.RData" save the MCMC posterior samples under different degrees of freedom, ranging from B0=3 to B0=10, and "ELPD_Values.RData" saves the corresponding elpd values for all scenarios.

**In the folder "Bayes_ITR_Supplement_E1":**

* The R script "Supplement_E1_Main.R" reproduces the additional simulation studies in Supplementary Material Section E1, including Figures E1 and E2 in the supplementary material.

* The R scripts "Simu_Data_Generate_I(a).R" and "Simu_Data_Generate_I(b).R" generate the simulated datasets for Scenarios I(a) and I(b), respectively, in the additional simulation studies. "MCMC_R_Functions.R" provides the R functions used for MCMC, "MCMC_Rcpp_Functions.cpp" provides the corresponding Rcpp functions, and "MCMC_Main.R" contains the main MCMC function.

* The R data files "MCMC_Results_I(a)_Mode1.RData" and "MCMC_Results_I(a)_Mode2.RData" save the MCMC posterior samples for the two modes under Scenario I(a), while "MCMC_Results_I(b)_Mode1.RData" and "MCMC_Results_I(b)_Mode2.RData" save the MCMC posterior samples for the two modes under Scenario I(b), in the additional simulation studies.

**In the folder "Bayes_ITR_Supplement_E2":**

* The R script "Supplement_E2_Main.R" reproduces the additional simulation studies in Supplementary Material Section E2, including Table E3 in the supplementary material.

* The R script "Simu_Data_Generate.R" generates the simulated datasets for Supplementary Material Section E2. "MCMC_R_Functions.R" provides the R functions used for MCMC, "MCMC_Rcpp_Functions.cpp" provides the corresponding Rcpp functions, and "MCMC_Main.R" contains the main MCMC function.

* The R data files "MCMC_Results_S1.RData" through "MCMC_Results_S6.RData" save the MCMC posterior samples for all six sub-scenarios in the additional simulation studies presented in Supplementary Material Section E2.

**In the folder "Bayes_ITR_Supplement_E3":**

* The R script "Supplement_E3_Main.R" reproduces the additional simulation studies in Supplementary Material Section E3, including Figure E4 in the supplementary material.

* The R script "Simu_Data_Generate.R" generates the simulated datasets for Scenario II, and "Select_Partial_Data.R" selects partial data to create the sub-scenarios II(a), II(b), and II(c) in the additional simulation studies. "MCMC_R_Functions.R" provides the R functions used for MCMC, "MCMC_Rcpp_Functions.cpp" provides the corresponding Rcpp functions, and "MCMC_Main.R" contains the main MCMC function. "SGD_Functions.R" provides the R functions used for SGD, and "SGD_Main_Conservative.R" contains the main SGD function for the proposed conservative policy optimization.

* The R data files "SGD_Results_Alpha.RData", "SGD_Results_Kappa.RData", and "SGD_Results_Sigma.RData" save the SGD results under different choices of the hyperparameters alpha, kappa, and sigma in the proposed policy optimization method, respectively.

**In the folder "Bayes_ITR_Supplement_E4":**

* The R script "Supplement_E4_Main.R" reproduces the additional simulation studies in Supplementary Material Section E4, including Figure E5 in the supplementary material.

* The R script "Simu_Data_Generate.R" generates the simulated datasets for Scenario II, and "Select_Partial_Data.R" selects partial data to create the sub-scenarios II(a), II(b), II(c), and II(d) in the additional simulation studies. "MCMC_R_Functions.R" provides the R functions used for MCMC, "MCMC_Rcpp_Functions.cpp" provides the corresponding Rcpp functions, and "MCMC_Main.R" contains the main MCMC function. "SGD_Functions.R" provides the R functions used for SGD, while "SGD_Main_Conservative.R" and "SGD_Main_Conservative_Mean.R" contain the main SGD functions for the proposed and alternative conservative policy optimization methods, respectively.

* The R data files "SGD_Results_II(a).RData", "SGD_Results_II(b).RData", "SGD_Results_II(c).RData", and "SGD_Results_II(d).RData" save the SGD results for sub-scenarios II(a), II(b), II(c), and II(d), respectively.
