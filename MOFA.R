## ----------------------------------------------------------------
# Set CRAN repository
options(repos = c(CRAN = "https://cloud.r-project.org/"))

# Install BiocManager if not already installed
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Load BiocManager
library(BiocManager)

# Install MOFA2 package
BiocManager::install("MOFA2", force = TRUE)

# Install devtools package
install.packages("devtools")

# Install a specific version of ggplot2
devtools::install_version("ggplot2", version = "3.5.1")


## ----------------------------------------------------------------
# Load required libraries
library(devtools)
library(MOFA2)
library(data.table)
library(ggplot2)
library(tidyverse)
library(reticulate)
library(randomForest)
library(utils)
library(survival)
library(survminer)


## ----------------------------------------------------------------
# Configure Python environment
configured_python <- Sys.getenv("MOFA_PYTHON", unset = "")
if (nzchar(configured_python)) {
  message(sprintf("Using Python from MOFA_PYTHON: %s", configured_python))
  use_python(configured_python, required = TRUE)
} else {
  message("MOFA_PYTHON is not set; using reticulate default Python discovery.")
}

python_config <- tryCatch(
  py_config(),
  error = function(e) {
    stop(
      paste(
        "Unable to initialize Python via reticulate.",
        "Install Python and required packages, or set MOFA_PYTHON to a valid interpreter path.",
        sprintf("Original error: %s", conditionMessage(e))
      )
    )
  }
)

if (!py_module_available("mofapy2")) {
  stop(
    paste(
      "Python was found, but the MOFA backend package 'mofapy2' is unavailable.",
      "Install it in the active Python environment (for example: pip install mofapy2)",
      "or set MOFA_PYTHON to a Python interpreter where 'mofapy2' is installed."
    )
  )
}


## ----------------------------------------------------------------
# Resolve project root and dataset paths robustly
script_path <- tryCatch(normalizePath(sys.frame(1)$ofile), error = function(e) NA)
if (is.na(script_path)) {
  script_path <- normalizePath(getwd())
}
project_root <- if (file.info(script_path)$isdir) script_path else dirname(script_path)
data_path <- file.path(project_root, "Dataset", "cll_data")
meta_path <- file.path(project_root, "Dataset", "cll_metadata")

# Validate data and metadata directories
if (!dir.exists(data_path)) {
  stop(
    sprintf(
      "Expected data directory does not exist: %s\nCreate the directory and place view CSV files under Dataset/cll_data/{Drugs,Methylation,mRNA,Mutations}.",
      data_path
    )
  )
}

if (!dir.exists(meta_path)) {
  stop(
    sprintf(
      "Expected metadata directory does not exist: %s\nCreate the directory and place one metadata CSV under Dataset/cll_metadata/.",
      meta_path
    )
  )
}

# Validate expected view CSV files
expected_views <- c("Drugs", "Methylation", "mRNA", "Mutations")
expected_view_files <- file.path(data_path, expected_views, paste0(expected_views, ".csv"))
missing_view_files <- expected_view_files[!file.exists(expected_view_files)]
if (length(missing_view_files) > 0) {
  stop(
    sprintf(
      "Missing expected view CSV file(s):\n- %s\nExpected files:\n- %s",
      paste(missing_view_files, collapse = "\n- "),
      paste(expected_view_files, collapse = "\n- ")
    )
  )
}

# Validate metadata CSV selection is deterministic and unique
meta_files <- sort(list.files(path = meta_path, pattern = "\\.csv$", full.names = TRUE))
if (length(meta_files) != 1) {
  stop(
    sprintf(
      "Expected exactly 1 metadata CSV in %s, found %d.\nEnsure only one metadata file exists or rename extras outside this directory.",
      meta_path,
      length(meta_files)
    )
  )
}

# Configure MOFA model source
# Allowed values:
# - "train": train a new model and save outputs.
# - "load_local": load an existing model from local .rds or .hdf5.
# - "load_remote_demo": load the public CLL demo model from EBI (demo-only path).
model_source <- "train"
allowed_model_sources <- c("train", "load_local", "load_remote_demo")
if (!model_source %in% allowed_model_sources) {
  stop(
    sprintf(
      "Invalid model_source '%s'. Allowed values are: %s",
      model_source,
      paste(allowed_model_sources, collapse = ", ")
    )
  )
}
message(sprintf("MOFA model source active: %s", model_source))

# Output/input paths used by train/load_local modes
mofa_hdf5_path <- file.path(project_root, "Dataset", "MOFA2_CLL.hdf5")
mofa_rds_path <- file.path(project_root, "MOFA2_CLL.rds")

# Explicit opt-in required for remote demo mode
allow_remote_demo <- identical(tolower(Sys.getenv("MOFA_ALLOW_REMOTE_DEMO", unset = "false")), "true")

# Initialize a named list to store omics data (deterministic view mapping)
CLL_data <- setNames(vector("list", length(expected_views)), expected_views)

# Loop through expected views in deterministic order
for (view_name in expected_views) {
  view_file <- file.path(data_path, view_name, paste0(view_name, ".csv"))
  omic <- read.csv(view_file, check.names = FALSE)

  if (!"X" %in% colnames(omic)) {
    stop(
      sprintf(
        "View '%s' is missing the feature ID column named 'X': %s",
        view_name,
        view_file
      )
    )
  }

  rownames(omic) <- omic$X
  omic$X <- NULL
  CLL_data[[view_name]] <- omic
}

CLL_data


## ----------------------------------------------------------------
# Read metadata CSV file
metadata <- read.csv(meta_files[1])

if (!"sample" %in% colnames(metadata)) {
  stop(
    sprintf(
      "Metadata file must include a 'sample' column: %s",
      meta_files[1]
    )
  )
}

if (anyDuplicated(metadata$sample) > 0) {
  duplicated_samples <- unique(metadata$sample[duplicated(metadata$sample)])
  stop(
    sprintf(
      "Metadata contains duplicate sample IDs in column 'sample': %s",
      paste(duplicated_samples, collapse = ", ")
    )
  )
}

metadata

# Convert each omic data frame to a data matrix
for (i in seq_along(CLL_data)) {
  CLL_data[[i]] <- data.matrix(CLL_data[[i]])
}


## ----------------------------------------------------------------
# Create a MOFA object from omics data list
MOFAobject <- create_mofa(CLL_data)
MOFAobject


## ----------------------------------------------------------------
# Plot an overview of the input data
plot_data_overview(MOFAobject)


## ----------------------------------------------------------------
# Get default data options
data_opts <- get_default_data_options(MOFAobject)
data_opts

# Get default model options and set number of factors
model_opts <- get_default_model_options(MOFAobject)
model_opts$num_factors <- 15
model_opts

# Get default training options and set convergence mode and seed
train_opts <- get_default_training_options(MOFAobject)
train_opts$convergence_mode <- "slow"
train_opts$seed <- 42
train_opts


## ----------------------------------------------------------------
# Prepare the MOFA object with set options
if (model_source == "train") {
  MOFAobject <- prepare_mofa(
    MOFAobject,
    data_options = data_opts,
    model_options = model_opts,
    training_options = train_opts
  )

  MOFAobject <- run_mofa(MOFAobject, outfile = mofa_hdf5_path, use_basilisk = FALSE)
  saveRDS(MOFAobject, mofa_rds_path)
  message(sprintf("Trained MOFA model and saved outputs:\n- %s\n- %s", mofa_hdf5_path, mofa_rds_path))
} else if (model_source == "load_local") {
  if (file.exists(mofa_rds_path)) {
    MOFAobject <- readRDS(mofa_rds_path)
    message(sprintf("Loaded local MOFA model from RDS: %s", mofa_rds_path))
  } else if (file.exists(mofa_hdf5_path)) {
    MOFAobject <- load_model(mofa_hdf5_path)
    message(sprintf("Loaded local MOFA model from HDF5: %s", mofa_hdf5_path))
  } else {
    stop(
      sprintf(
        "model_source='load_local' but no local model file found.\nChecked:\n- %s\n- %s",
        mofa_rds_path,
        mofa_hdf5_path
      )
    )
  }
} else if (model_source == "load_remote_demo") {
  if (!allow_remote_demo) {
    stop(
      paste(
        "Remote demo loading is disabled by default.",
        "Set MOFA_ALLOW_REMOTE_DEMO=true to explicitly opt in to model_source='load_remote_demo'."
      )
    )
  }
  MOFAobject <- readRDS(url("http://ftp.ebi.ac.uk/pub/databases/mofa/cll_vignette/MOFA2_CLL.rds"))
  message("Loaded remote demo MOFA model from EBI FTP.")
}


## ----------------------------------------------------------------
# Inspect slot names of the MOFA object
slotNames(MOFAobject)

# Inspect names of the data
names(MOFAobject@data)

# Get dimensions of data in 'Drugs' view
dim(MOFAobject@data$Drugs$group1)

# Inspect names of expectations
names(MOFAobject@expectations)

# Get dimensions of latent factors (Z)
dim(MOFAobject@expectations$Z$group1)

# Get dimensions of weights (W)
dim(MOFAobject@expectations$W$mRNA)


## ----------------------------------------------------------------
# Check if sample names are consistent between MOFA and metadata
mofa_samples <- unlist(samples_names(MOFAobject))

missing_in_metadata <- setdiff(mofa_samples, metadata$sample)
if (length(missing_in_metadata) > 0) {
  stop(
    sprintf(
      "Metadata is missing %d sample(s) required by MOFA object: %s",
      length(missing_in_metadata),
      paste(missing_in_metadata, collapse = ", ")
    )
  )
}

extra_in_metadata <- setdiff(metadata$sample, mofa_samples)
if (length(extra_in_metadata) > 0) {
  stop(
    sprintf(
      "Metadata contains %d sample(s) not present in MOFA object: %s",
      length(extra_in_metadata),
      paste(extra_in_metadata, collapse = ", ")
    )
  )
}

# Reorder metadata to match MOFA sample order and verify strict alignment
metadata <- metadata[match(mofa_samples, metadata$sample), , drop = FALSE]

if (!identical(metadata$sample, mofa_samples)) {
  mismatch_idx <- which(metadata$sample != mofa_samples)[1]
  stop(
    sprintf(
      "Sample ordering mismatch after metadata alignment at position %d: metadata='%s' vs MOFA='%s'",
      mismatch_idx,
      metadata$sample[mismatch_idx],
      mofa_samples[mismatch_idx]
    )
  )
}

# Add metadata to MOFA object
samples_metadata(MOFAobject) <- metadata

# Plot correlation between factors
plot_factor_cor(MOFAobject)


## ----------------------------------------------------------------
# Plot variance explained by each factor
plot_variance_explained(MOFAobject, max_r2 = 15)


## ----------------------------------------------------------------
# Plot total variance explained (second plot output)
plot_variance_explained(MOFAobject, plot_total = TRUE)[[2]]


## ----------------------------------------------------------------
# Correlate factors with covariates
correlate_factors_with_covariates(MOFAobject,
                                  covariates = c("Gender", "died", "age"),
                                  plot = "log_pval"
)


## ----------------------------------------------------------------
# Plot factor 1 colored by 'Factor1'
plot_factor(MOFAobject,
            factors = 1,
            color_by = "Factor1"
)


## ----------------------------------------------------------------
# Plot top weights in the 'Mutations' view for factor 1
plot_weights(MOFAobject,
             view = "Mutations",
             factor = 1,
             nfeatures = 10,
             scale = TRUE
)


## ----------------------------------------------------------------
# Plot top weights in the 'Mutations' view for factor 1
plot_top_weights(MOFAobject,
                 view = "Mutations",
                 factor = 1,
                 nfeatures = 10,
                 scale = TRUE
)


## ----------------------------------------------------------------
# Plot factor 1 colored by IGHV, with violin plots
plot_factor(MOFAobject,
            factors = 1,
            color_by = "IGHV",
            add_violin = TRUE,
            dodge = TRUE
)


## ----------------------------------------------------------------
# Plot factor 1 colored by Gender, with violin plots
plot_factor(MOFAobject,
            factors = 1,
            color_by = "Gender",
            dodge = TRUE,
            add_violin = TRUE
)


## ----------------------------------------------------------------
# Plot weights in the 'mRNA' view for factor 1
plot_weights(MOFAobject,
             view = "mRNA",
             factor = 1,
             nfeatures = 10
)


## ----------------------------------------------------------------
# Plot scatter of data points colored by IGHV for top positive feature
plot_data_scatter(MOFAobject,
                  view = "mRNA",
                  factor = 1,
                  features = 4,
                  sign = "positive",
                  color_by = "IGHV"
) + labs(y = "RNA expression")


## ----------------------------------------------------------------
# Plot scatter of data points colored by IGHV for top negative feature
plot_data_scatter(MOFAobject,
                  view = "mRNA",
                  factor = 1,
                  features = 4,
                  sign = "negative",
                  color_by = "IGHV"
) + labs(y = "RNA expression")


## ----------------------------------------------------------------
# Plot heatmap for top 25 features, scale by row
plot_data_heatmap(MOFAobject,
                  view = "mRNA",
                  factor = 1,
                  features = 25,
                  cluster_rows = FALSE, cluster_cols = FALSE,
                  show_rownames = TRUE, show_colnames = FALSE,
                  scale = "row"
)


## ----------------------------------------------------------------
# Plot heatmap for top 25 features with denoising, scale by row
plot_data_heatmap(MOFAobject,
                  view = "mRNA",
                  factor = 1,
                  features = 25,
                  denoise = TRUE,
                  cluster_rows = FALSE, cluster_cols = FALSE,
                  show_rownames = TRUE, show_colnames = FALSE,
                  scale = "row"
)


## ----------------------------------------------------------------
# Plot weights in the 'Mutations' view for factor 3
plot_weights(MOFAobject,
             view = "Mutations",
             factor = 3,
             nfeatures = 10,
             abs = FALSE
)


## ----------------------------------------------------------------
# Plot factor 3 colored by trisomy12, with violin plots
plot_factor(MOFAobject,
            factors = 3,
            color_by = "trisomy12",
            dodge = TRUE,
            add_violin = TRUE
)


## ----------------------------------------------------------------
# Plot scatter of data points colored by trisomy12 for top positive feature
plot_data_scatter(MOFAobject,
                  view = "Drugs",
                  factor = 3,
                  features = 4,
                  sign = "positive",
                  color_by = "trisomy12"
) + labs(y = "Drug response (cell viability)")


## ----------------------------------------------------------------
# Plot heatmap for top 25 features with denoising, scale by row
plot_data_heatmap(MOFAobject,
                  view = "mRNA",
                  factor = 3,
                  features = 25,
                  denoise = TRUE,
                  cluster_rows = TRUE, cluster_cols = FALSE,
                  show_rownames = TRUE, show_colnames = FALSE,
                  scale = "row"
)


## ----------------------------------------------------------------
# Plot factors 1 and 3 with color and shape
p <- plot_factors(MOFAobject,
                  factors = c(1, 3),
                  color_by = "IGHV",
                  shape_by = "trisomy12",
                  dot_size = 2.5,
                  show_missing = TRUE
)

p <- p +
  geom_hline(yintercept = -1, linetype = "dashed") +
  geom_vline(xintercept = (-0.5), linetype = "dashed")

print(p)


## ----------------------------------------------------------------
# Load the random forest library
suppressPackageStartupMessages(library(randomForest))

# Extract factors 1 and 2
df <- as.data.frame(get_factors(MOFAobject, factors = c(1, 2))[[1]])


## ----------------------------------------------------------------
# Add IGHV data as a factor
df$IGHV <- as.factor(MOFAobject@samples_metadata$IGHV)
# Train random forest model on samples with non-missing IGHV
model.ighv <- randomForest(IGHV ~ ., data = df[!is.na(df$IGHV), ], ntree = 10)
# Remove IGHV
df$IGHV <- NULL


## ----------------------------------------------------------------
# Predict IGHV using the model
MOFAobject@samples_metadata$IGHV.pred <- stats::predict(model.ighv, df)


## ----------------------------------------------------------------
# Add trisomy12 data as a factor
df$trisomy12 <- as.factor(MOFAobject@samples_metadata$trisomy12)
# Train random forest model on samples with non-missing trisomy12
model.trisomy12 <- randomForest(trisomy12 ~ ., data = df[!is.na(df$trisomy12), ], ntree = 10)
# Remove trisomy12
df$trisomy12 <- NULL

# Predict trisomy12 using the model
MOFAobject@samples_metadata$trisomy12.pred <- stats::predict(model.trisomy12, df)
# Add new column based on missing IGHV status 
MOFAobject@samples_metadata$IGHV.pred_logical <- c("True", "Predicted")[as.numeric(is.na(MOFAobject@samples_metadata$IGHV)) + 1]

# plot predicted factors
p <- plot_factors(MOFAobject,
                  factors = c(1, 3),
                  color_by = "IGHV.pred",
                  shape_by = "IGHV.pred_logical",
                  dot_size = 2.5,
                  show_missing = TRUE
)
p <- p +
  geom_hline(yintercept = -1, linetype = "dashed") +
  geom_vline(xintercept = (-0.5), linetype = "dashed")

print(p)


## ----------------------------------------------------------------
# Extract weights
weights <- get_weights(MOFAobject,
                       views = "all",
                       factors = "all",
                       as.data.frame = TRUE
)
head(weights)


## ----------------------------------------------------------------
# Extract factors
factors <- get_factors(MOFAobject,
                       factors = "all",
                       as.data.frame = TRUE
)
head(factors)


## ----------------------------------------------------------------
# Extract data
data <- get_data(MOFAobject,
                 views = "all",
                 as.data.frame = TRUE
)
head(data)


## ----------------------------------------------------------------
# Impute missing values in MOFA object
MOFAobject <- impute(MOFAobject)

# Show a subset of the imputed values
MOFAobject@data$mRNA[[1]][1:5, 190:195]


## ----------------------------------------------------------------
# Prepare Surv object for survival analysis
SurvObject <- Surv(MOFAobject@samples_metadata$TTT, MOFAobject@samples_metadata$treatedAfter)
# Extract latent factors
Z <- get_factors(MOFAobject)[[1]]
# Fit a Cox proportional hazards model
fit <- coxph(SurvObject ~ Z)
fit


## ----------------------------------------------------------------
# Summarize the Cox model
s <- summary(fit)
coef <- s[["coefficients"]]

# Create data frame for plotting
df <- data.frame(
  factor = factor(rownames(coef), levels = rev(rownames(coef))),
  p = coef[, "Pr(>|z|)"],
  coef = coef[, "exp(coef)"],
  lower = s[["conf.int"]][, "lower .95"],
  higher = s[["conf.int"]][, "upper .95"]
)

# Create plot of hazard ratios
ggplot(df, aes(x = factor, y = coef, ymin = lower, ymax = higher)) +
  geom_pointrange(col = "#619CFF") +
  coord_flip() +
  scale_x_discrete() +
  labs(y = "Hazard Ratio", x = "") +
  geom_hline(aes(yintercept = 1), linetype = "dotted") +
  theme_bw()


## ----------------------------------------------------------------
# Create a dataframe for survival analysis
df <- data.frame(
  time = SurvObject[, 1],
  event = SurvObject[, 2],
  Z1 = Z[, 1]
)
# Get a cutpoint for survival analysis based on factor 1
cut <- surv_cutpoint(df, variables = "Z1")
# Create a binary factor based on the cutpoint
df$FactorCluster <- df$Z1 > cut$cutpoint$cutpoint
# Fit survival model
fit <- survfit(Surv(time, event) ~ FactorCluster, df)
# Plot survival curves
ggsurvplot(fit,
           data = df,
           conf.int = TRUE, pval = TRUE,
           fun = function(y) y * 100,
           legend = "top", legend.labs = c(paste("low LF 1"), paste("high LF 1")),
           xlab = "Time to treatment", ylab = "Survival probability (%)", title = "Factor 1"
)$plot


## ----------------------------------------------------------------
# Print session info
sessionInfo()
