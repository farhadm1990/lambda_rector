# Lambda Rector
## An R package to correct relatvie abundance  of sequeneincg reads into 16S rRNA gene copy-number based on an internal Lambda Phage standard.

**This is the supporting package for paper DOIXXXXX**

# Installation

## 1. Install and library devtools pakcage on your machine
```R
if(!require("devtools")){install.packages("devtools")}
library(devtools)
```

## 2. Download and install `lambda_rector`
```R
devtools::install_github("farhadm1990/lambda_rector")
library(lambda.rector)
```
## 3. Creating a phyloseq object based on test dataset
```R
count = read.table("./lambda.rector/tests/count_test.tsv")
metadata = read.table("./lambda.rector/tests/metadata_test.tsv")
taxa = read.table("./lambda.rector/tests/taxonomy.tsv", header = F) %>% column_to_rownames("V1") %>% 
tidyr::separate( col = "V2", sep = ";", into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")) %>% apply( 2, function(x){gsub("[a-zA-Z]+__", "", x)}) # parssing the taxa column and tidying the names

#IMPORTANT: giving an arbiterary name to the lambda standard to be passed on to the function later. 

cbind(taxa[,1] %>% data.frame() %>% rename("."="king"), count) %>% group_by(king) %>% summarise_all(sum)

taxa[taxa[,1]=="Unassigned",] <- "Lambda"

ps = phyloseq(otu_table(count, taxa_are_row = TRUE), tax_table(as(taxa, "matrix")), sample_data(metadata))

```

## 4. Rinning `lambda_rector` function
```R
test = lambda_rector(
  ps,
  lamba_id = "Lambda",
  out_path = "./",
  negative_cont = NULL,
  negative_filt = TRUE,
  rare_depth = 10000,
  taxa_level = "Species",
  std_threshold = 1.4
)

# This will return a list of differnt phyloseq objects and saves the output plots
```
## 3. Output examples

![plot1](https://github.com/farhadm1990/lambda_rector/blob/main/pix/plot_without_bad_samples.jpeg)

### Fig1. An example of filtering output by the package on suspicious samples.

![plot2](https://github.com/farhadm1990/lambda_rector/blob/main/pix/Order_relative.jpeg)

### Fig2. Relative abundance of 16S rRNA gene sequencing reads at order level in different Mock and Lambda concentrations.

![plot3](https://github.com/farhadm1990/lambda_rector/blob/main/pix/Order_copy_number.jpeg)

### Fig3. 16S rRNA gene Copy-number corrected reads at order level in different Mock and Lambda concentrations.