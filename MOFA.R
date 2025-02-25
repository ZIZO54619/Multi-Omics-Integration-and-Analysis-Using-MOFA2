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
py_config()

# Use specific Python executable
use_python("../python.exe")


## ----------------------------------------------------------------
# Define data path
path <- "../Dataset/cll_data"

# List files in the data path
files <- list.files(path)
files

# Initialize an empty list to store omics data
CLL_data <- list()

# Loop through each directory in the main path
for (file_name in files) {
    sub_path <- file.path(path, "/", file_name)
    # Get list of CSV files within subfolders
    sub_files <- list.files(path = sub_path, pattern = "\\.csv$", full.names = TRUE)
    # Loop through the CSV files
    for (i in seq_along(sub_files)) {
      omic <- read.csv(sub_files[i])
      rownames(omic) <- omic$X
      omic$X <- NULL
      CLL_data[[length(CLL_data) + 1]] <- omic
    }
  }

# Set names of omics data list
names(CLL_data) <- c("Drugs", "Methylation", "mRNA", "Mutations")
CLL_data


## ----------------------------------------------------------------
# Define metadata path
meta_path <- "../Dataset/cll_metadata"

# List metadata files
meta_files <- list.files(path = meta_path, pattern = "\\.csv$", full.names = TRUE)
meta_files

# Read metadata CSV file
metadata <- read.csv(meta_files[1])
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
MOFAobject <- prepare_mofa(MOFAobject,
                           data_options = data_opts,
                           model_options = model_opts,
                           training_options = train_opts
)
# Run MOFA analysis
MOFAobject <- run_mofa(MOFAobject, outfile = "../Dataset/MOFA2_CLL.hdf5", use_basilisk = FALSE)

# Save MOFA object
saveRDS(MOFAobject, "MOFA2_CLL.rds")

# Load MOFA object from an external URL
MOFAobject <- readRDS(url("http://ftp.ebi.ac.uk/pub/databases/mofa/cll_vignette/MOFA2_CLL.rds"))


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
stopifnot(all(sort(metadata$sample) == sort(unlist(samples_names(MOFAobject)))))

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

