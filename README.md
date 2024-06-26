<div style="display: flex; justify-content: center; align-items: center; height: 100vh;">
    <img width="300" height="300" src="https://raw.githubusercontent.com/farhadm1990/lambda_rector/main/pix/logo.jpg" alt="Logo">
</div>



# Lambda Rector
### An R package to correct relatvie abundance  of sequeneincg reads into 16S rRNA gene copy-number based on an internal Lambda Phage standard.

**This is the supporting package for paper DOIXXXXX**

# Installation

## 1. Install and library devtools pakcage on your machine
```R
if(
    !require("devtools")
    ){
        install.packages("devtools")
        }
library(
    devtools
    )
```



## 2. Download and install `lambda_rector`
```R
devtools::install_github(
    "farhadm1990/lambda_rector",
    dependencies=TRUE
    )
library(
    lambda.rector
    )
```

### Alternatively you can clone this repository

```bash
git clone https://github.com/farhadm1990/lambda_rector.git


devtools::install_local("./lambda_rector")


```


## 3. Creating a phyloseq object based on test dataset
```R
# In the test subdirecotry of lambda_rector you can find test dataset.
count = read.table(
    "./lambda.rector/tests/count_test.tsv"
    )
metadata = read.table(
    "./lambda.rector/tests/metadata_test.tsv"
    )

taxa = read.table(
    "./lambda.rector/tests/taxonomy.tsv", header = F
    ) %>% 
    column_to_rownames(
        "V1"
        ) %>% 
tidyr::separate( 
    col = "V2", sep = ";", 
    into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
    ) %>% apply( 
        2, function(x){
            gsub("[a-zA-Z]+__", "", x)
        }
        ) # parssing the taxa column and tidying the names

#IMPORTANT: giving an arbiterary name to the lambda standard to be passed on to the function later. 

cbind(
    taxa[,1] %>% data.frame() %>% rename("."="king"), count
    ) %>% 
    group_by(king) %>% 
    summarise_all(sum)

taxa[taxa[,1]=="Unassigned",] <- "Lambda"

ps = phyloseq(
    otu_table(
        count, taxa_are_row = TRUE), 
        tax_table(as(taxa, "matrix")), 
        sample_data(metadata)
        )

```

## 4. Rinning `lambda_rector` function
```R
test_ps = lambda_rector(
                        ps, 
                          lambda_id = "Lambda",
                          singletone_threshold = 1, 
                          out_path = "./", 
                          negative_cont = NULL,
                          negative_filt = FALSE, 
                          rare_depth = 10000, 
                          taxa_level = "Kingdom", 
                          std_threshold = 1.48,
                          rel_abund_threshold = c(0.01,0.99)
                          )


# This will return a list of differnt phyloseq objects and saves the output plots



# Extracting the copy-corrected talbe
cbind(sample_data(test_ps$copy_corrected_ps), test_ps$copy_corrected_matrix) %>% 
data.frame() %>% 
rownames_to_column("barcodes") %>% 
select(-loaded_copy_lambda, -samp_id, -volum_mock, -volum_lambda) %>% 
group_by(barcodes, lambda_ng_ul, mock_ng_ul) %>% 
summarise_all(mean) 

```
## 3. Output examples

![plot1](https://github.com/farhadm1990/lambda_rector/blob/main/pix/plot_with_bad_samples.jpeg)
![plot1](https://github.com/farhadm1990/lambda_rector/blob/main/pix/plot_without_bad_samples.jpeg)
### Fig1. An example of filtering output by the package on suspicious samples.

![plot2](https://github.com/farhadm1990/lambda_rector/blob/main/pix/Genus_relative.jpeg)

### Fig2. Relative abundance of 16S rRNA gene sequencing reads at Genus level in different Mock and Lambda concentrations.

![plot3](https://github.com/farhadm1990/lambda_rector/blob/main/pix/Genus_copy_number_noLambda.jpeg)

### Fig3. 16S rRNA gene Copy-number corrected reads at Genus level in different Mock and Lambda concentrations.