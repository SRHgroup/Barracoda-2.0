# Barracoda-2.0

This repository contains Barracoda-2.0 a tool for analysing DNA barcode sequencing data.

# Quick Start

### Install Anaconda:

Install Anaconda:
If Anaconda is not installed on your system, then install it using the link below.

https://repo.anaconda.com/archive/

### Clone the repository:

Then clone or download the repository: 

```
git clone https://github.com/SRHgroup/Barracoda-2.0.git
cd Barracoda-2.0
```

### Install program dependensies

Program dependensies can be found in environment.yml.
Using conda alows you have all program dependensies in a conda enviroment called barracoda_env

```
conda env create -f environment.yml
```

Activate conda environment

```
conda activate barracoda_env
```

### Modify barracoda paths
Modify the "barracoda_dir" path within the barracoda-2.0.sh script

### Test the installation:

Change permission to barracoda-2.0.sh to make it exicutable. 

```
chmod 775 barracoda-2.0.sh
```

### Test Barracoda-2.0

```
./barracoda-2.0.sh -h 
```

# Example Usage

Barracoda can run using either Illumina/IonTorrent or Nanopore input data

### Run Barracoda-2.0 using Nanopore data 

```
conda activate barracoda_env
./barracoda-2.0.sh -m data/Nanopore_data/sample_identification.xlsx -a data/Nanopore_data/barcode_annotation.xlsx -A data/Nanopore_data/sample-identification-tag.fasta -B GAAGTTCCAGCCAGCGTCACAGTTT -C 6 -D data/Nanopore_data/EpitopeTagA.fasta -E GGTCAGCATCATTTCC -F data/Nanopore_data/EpitopeTagB.fasta -G 6 -H CAATCTTGAGCGTGACTTAAG -f data/Nanopore_data/Nanopore_data -o test_results/Nanopore_results -n
```

Add the plate setup using `-p data/Nanopore_data/barcode-plate-setup.xlsx` (optional)

### Run Barracoda-2.0 using IonTorrent data 
```
conda activate barracoda_env
./barracoda-2.0.sh -m data/IonTorrent_data/sample_identification.xlsx -a data/IonTorrent_data/barcode_annotation.xlsx -A data/IonTorrent_data/sample-identification-tag.fasta -B GAAGTTCCAGCCAGCGTCACAGTTT -C 6 -D data/IonTorrent_data/EpitopeTagA.fasta -E GGTCAGCATCATTTCC -F data/IonTorrent_data/EpitopeTagB.fasta -G 6 -H GTTATCGGCTCGTTCACACTCGA -f data/IonTorrent_data/IonTorrent_data.fastq -o test_results/IonTorrent_results
```

Add the plate setup using `-p data/IonTorrent_data/barcode-plate-setup.xlsx` (optional)


# Plot Results

The script `plot_results.R` can be used as an example of how to visualize the results.
Example plots produced using this script can be found in the `plots/` folder.





