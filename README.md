# Lambda Rector
## An R package to correct relatvie abundance  of sequeneincg reads into 16S rRNA gene copy-number based on an internal Lambda Phage standard.

**This is the supporting package for paper DOIXXXXX**
# Installation

## 1. Isntall and library devtools pakcage on your machine
```R
install.packages("devtools")
library(devtools)
```

## 2. Download and isntall `lambda_rector`
```R
devtools::install_github("farhadm1990/lambda_rector")
library(lambda.rector)
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
  taxa_level = "Order"
)

# This will return a list of differnt phyloseq objects and saves the output plots
```
## 3. Output examples

![plot1](https://github.com/farhadm1990/lambda_rector/blob/main/pix/plot_without_bad_samples.jpeg)

### Fig2. An example of filtering output of by the package on suspicious samples.

![plot2](https://github.com/farhadm1990/lambda_rector/blob/main/pix/Order_relative.jpeg)

### Fig3. Relative abundance of 16S rRNA gene sequencing reads at order level in different Mock and Lambda concentrations.

![plot3](https://github.com/farhadm1990/lambda_rector/blob/main/pix/Order_copy_number.jpeg)

### Fig3. 16S rRNA gene Copy-number corrected reads at order level in different Mock and Lambda concentrations.