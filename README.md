# Lambda Rector
## An R package to correct relatvie abundance  of sequeneincg reads into 16S rRNA gene copy-number based on an internal Lambda Phage standard.

**NOTE: this package is only tunned for four filtersets, Biomass, pH, Riboflavine, and DO**
# Installation

## 1. Isntall and library devtools pakcage on your machine
```R
install.packages("devtools")
library(devtools)
```

## 2. Download and isntall biolecter
```R
devtools::install_github("farhadm1990/lambda_rector")
library(lambda_rector)
```

## 3. Rinning `lambda_rector` function
```R
test = lambda_rector(
  ps,
  lamba_id = "Lambda",
  out_path = "./",
  negative_cont = NULL,
  negative_filt = TRUE,
  rare_depth = 10000,
  taxa_level = "Kingdom"
)

# This will return a list of differnt phyloseq objects and saves the output plots
```




