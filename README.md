# Missing Data Sensitivity Analysis for Alcohol Research

This repository contains code and example data used to implement the MNAR sensitivity analyses described in the paper.

## Contents
- `Sensitivity Tutorial with Artifical Data.R` – executable R sript
- `Sensitivity Tutorial with Artifical Data.Rmd` – executable R markdown sript with detailed annotations
- `Blimp Studio Scripts` - folder containing scripts for Blimp GUI
- `.imp` files – Blimp model specifications
- `mnar.csv` – example dataset

## How to run
1. Install Blimp and the rblimp package (wwww.appliedmissingdata.com/blimp and CRAN)
2. Open the R or Rmd script to run sensitivity analyses in rblimp package
3. Open Blimp Studio GUI to run .imp files.
4. The R script reads the raw data .csv from this site. Save the .csv in the same directory as the .imp files if using Blimp Studio.

## Notes
- Results will vary slightly due to MCMC sampling
- These scripts are provided for reproducibility and illustration
