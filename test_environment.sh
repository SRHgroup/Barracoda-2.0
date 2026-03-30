#!/bin/bash

echo "=== Testing environment tools ==="

# -------------------------
# Python
# -------------------------
echo "Python test:"
python - <<'PYTHON'
import sys

print("Python is working!")

# List of modules to check
modules = ["os", "gzip", "regex", "argparse", "zipfile", "shutil", "Bio.SeqIO"]

for mod in modules:
    try:
        __import__(mod.split('.')[0])
        print(f"Module '{mod}' is installed.")
    except ImportError:
        print(f"Module '{mod}' is NOT installed!")
PYTHON

echo

# -------------------------
# R
# -------------------------
echo "R test:"
Rscript - <<'RSCRIPT'
packages <- c("squash", "xlsx", "tidyverse", "ggplot2", "openxlsx",
              "reshape2", "dplyr", "data.table", "scales", "tools",
              "edgeR", "limma")

cat("R is working!\n")
for (pkg in packages) {
  if (suppressWarnings(require(pkg, character.only = TRUE))) {
    cat("Package", pkg, "is installed.\n")
  } else {
    cat("Package", pkg, "is NOT installed!\n")
  }
}
RSCRIPT

echo

# -------------------------
# Perl
# -------------------------
echo "Perl test:"
perl -e 'print "Perl is working!\n"'

# -------------------------
# Bowtie2
# -------------------------
echo "Bowtie2 test:"
bowtie2 --version | head -n1

# -------------------------
# GNU parallel
# -------------------------
echo "GNU parallel test:"
parallel --version | head -n1

# -------------------------
# Full paths of dependencies
# -------------------------
echo
echo "=== Full paths of all dependencies ==="
echo "Python:    $(which python)"
echo "R:         $(which R)"
echo "Perl:      $(which perl)"
echo "Bowtie2:   $(which bowtie2)"
echo "Parallel:  $(which parallel)"
