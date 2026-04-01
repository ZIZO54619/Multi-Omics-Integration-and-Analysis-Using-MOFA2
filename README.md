# MOFA2 Analysis Pipeline (CLL Multi-Omics)

This repository contains an R workflow (`MOFA.R`) for training or loading a MOFA2 model on CLL multi-omics data, then running downstream visualization and analysis.

## Repository layout

```text
.
├── MOFA.R
├── MOFA.md
├── Dataset/
│   ├── cll_data/
│   │   ├── Drugs/Drugs.csv
│   │   ├── Methylation/Methylation.csv
│   │   ├── mRNA/mRNA.csv
│   │   └── Mutations/Mutations.csv
│   └── cll_metadata/
│       └── <single_metadata_file>.csv
└── Plots/
```

## 1) Environment setup (run once)

`MOFA.R` does **not** install dependencies automatically. Install R and Python dependencies first.

### R packages

```r
install.packages(c(
  "data.table", "ggplot2", "tidyverse", "reticulate",
  "randomForest", "survival", "survminer"
))
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("MOFA2")
```

### Python package (for MOFA backend)

Install `mofapy2` in the Python environment that reticulate will use:

```bash
pip install mofapy2
```

## 2) Python selection

You can optionally pin the Python interpreter via `MOFA_PYTHON`.

```bash
# Linux/macOS
export MOFA_PYTHON=/opt/miniconda3/envs/mofa/bin/python
Rscript MOFA.R
```

```powershell
# Windows PowerShell
$env:MOFA_PYTHON="C:\Users\you\miniconda3\envs\mofa\python.exe"
Rscript .\MOFA.R
```

If `MOFA_PYTHON` is not set, reticulate default discovery is used.

## 3) Model source modes

In `MOFA.R`, choose one mode via `model_source`:

- `"train"`: train a new model and save outputs.
- `"load_local"`: load `MOFA2_CLL.rds` (or fallback HDF5) from local paths.
- `"load_remote_demo"`: load remote demo model (disabled unless `MOFA_ALLOW_REMOTE_DEMO=true`).

## 4) Runtime behavior and fail-fast checks

`MOFA.R` performs explicit checks before analysis:

- required R packages available
- Python initializes successfully via reticulate
- Python module `mofapy2` available
- expected dataset directories/files exist
- metadata file count is deterministic (exactly one CSV)
- each omics view contains feature-id column `X`
- metadata has `sample`, no duplicates, and aligns exactly with MOFA sample order

## 5) Run

From repo root:

```bash
Rscript MOFA.R
```

## Notes on scientific settings

The script intentionally keeps analytical choices unchanged (unless explicitly approved), including:

- `model_opts$num_factors <- 15`
- `train_opts$convergence_mode <- "slow"`
- covariates `Gender`, `died`, `age`
- random forest logic for missing IGHV/trisomy12
- survival endpoint/cutpoint flow (`TTT`, `treatedAfter`)
