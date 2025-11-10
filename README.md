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

### Prepare data files 

### Run Barracoda (-n for Nanopore) 
./barracoda-2.0.sh -f data/Nanopore/PAS39055_pass_barcode24_small.fastq.gz -m data/Nanopore/sample_id_table.xlsx -a data/Nanopore/barcode_annotations.xlsx -A data/Nanopore/sample_id_tags.fasta -B GAAGTTCCAGCCAGCGTCACAGTTT -C 6 -D data/Nanopore/a_epitope_tag.fasta -E GGTCAGCATCATTTCC -F data/Nanopore/b_epitope_tag.fasta -G 6 -H CAATCTTGAGCGTGACTTAAG -n

### Run Barracoda (Illumina)
./barracoda-2.0.sh -f data/Illumina/test_10k.fastq -m data/Illumina/sample-idenfication-table-small.xlsx -a data/Illumina/Barcode_annotations_small.xlsx -A data/Illumina/sample.fasta -B GAAGTTCCAGCCAGCGTCACAGTTT -C 6 -D data/Illumina/a.fasta -E GGTCAGCATCATTTCC -F data/Illumina/b.fasta -G 6 -H GTTATCGGCTCGTTCACACTCGA


