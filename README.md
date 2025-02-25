# MOFA2 Analysis Pipeline

## Overview
This script provides a complete workflow for analyzing multi-omics data using the **MOFA2** package in R. It includes data preprocessing, MOFA model training, and visualization of the results.

## Prerequisites
Ensure you have **R** installed and set up a suitable Python environment for MOFA2. The script installs necessary packages and dependencies.

## Steps in the Script

### 1. Install Required Packages
- Sets the CRAN repository.
- Installs `BiocManager` if not available.
- Installs MOFA2 and other required packages (`devtools`, `ggplot2`, etc.).

### 2. Load Required Libraries
- Loads essential libraries including `MOFA2`, `ggplot2`, `tidyverse`, `reticulate`, `randomForest`, `utils`, `survival`, and `survminer`.

### 3. Configure Python Environment
- Specifies a Python executable for MOFA2 operations.

### 4. Load Omics Data
- Reads multiple CSV files from `../Dataset/cll_data`.
- Stores them in a named list (`CLL_data`).

### 5. Load Metadata
- Reads metadata from `../Dataset/cll_metadata`.

### 6. Convert Data to Matrix Format
- Converts each omic dataset into a data matrix for MOFA2 compatibility.

### 7. Create and Prepare MOFA Model
- Initializes a **MOFA object** with the omics data.
- Configures model training options (number of factors, convergence mode, random seed).
- Runs MOFA2 analysis and saves the trained model.

### 8. Data Inspection
- Checks sample names consistency.
- Adds metadata to MOFA object.
- Inspects slot names, expectations, and data dimensions.

### 9. Data Visualization
- Plots correlation between factors.
- Plots variance explained by factors.
- Correlates factors with clinical covariates.
- Plots feature weights, heatmaps, and scatter plots for different views.

### 10. Factor Analysis
- Visualizes important factors using violin plots and heatmaps.
- Analyzes relationships between omics data and clinical parameters.

## Output
- The trained MOFA object is saved as `MOFA2_CLL.rds`.
- Various plots and visualizations are generated to interpret multi-omics relationships.

## Notes
- Ensure that the specified paths for data and metadata exist.
- Adjust the number of factors and training options based on the dataset.
- Use `reticulate` to properly configure Python dependencies if needed.

## Contact
For issues or questions, feel free to reach out!

