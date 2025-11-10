# Barracoda-2.0
This repository contains Barracoda-2.0 a tool for analysing DNA barcode sequencing data. 

### Program dependencies
install perl 

install python 

install bowtie2 

install R 

install GNU parallel

### Required R Packages

The following R libraries must be installed:

install.packages(c("squash", "xlsx", "tidyverse", "ggplot2", 
                   "openxlsx", "reshape2", "dplyr", "data.table", 
                   "scales", "tools"))


For edgeR and limma, install from Bioconductor:

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("edgeR", "limma"))

### Dowload Barracoda-2.0 
git clone git@github.com:SRHgroup/Barracoda-2.0.git

### Modify paths 
Modify paths within the barracoda-2.0.sh script

### Test Barracoda  
./barracoda-2.0.sh -h
