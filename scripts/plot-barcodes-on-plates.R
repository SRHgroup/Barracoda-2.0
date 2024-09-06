library(openxlsx)
library(reshape2)
library(ggplot2)
options(stringsAsFactors = FALSE, scipen = 1000)

# source("/home/tuba/marquard/git/Barracoda/webservice/src/scripts/functions.R")


#------- COMMAND LINE ARGUMENTS 

cArgs <- commandArgs(TRUE)


base_dir    <- cArgs[1]
plateFile   <- cArgs[2]
annotations <- cArgs[3]



#-------   GET PLATE SETUP
plates <- readPlateLayoutFromExcel(plateFile)


#-------   GET BARCODE ANNOTATIONS
anno <- read.annotation(annotations)
all_barcodes <- unlist(lapply(anno, rownames))


#-------- Load and plot the read count matrices

file_in <- "unfiltered-readcounts--with-clonality-red.All.txt"

out_dir <- file.path(base_dir, "plate-setup")
dir.create(out_dir)

file_out <- sub("\\..*$", "", file_in)

x <- read.delim(file.path(base_dir, file_in), as.is = TRUE, row.names = 1)

# Remove annotation columns
x <- x[,!colnames(x) %in% colnames(anno)]

# Remove the X that R previously added to the akeys
colnames(x) <- sub("^X", "", colnames(x))

plotDataOnPlates(x, paste0(out_dir, "/", file_out, "--plate-view"), plates, skip = 0, log = TRUE, invertSet = !rownames(x) %in% all_barcodes)   


#-------- Load and plot the log_fold_change and pvalues and plot them

for (dir in list.dirs(base_dir, recursive = FALSE)) {

	if(!file.exists(file.path(dir, "all-pvalues.txt"))) { next }

	out_dir <- file.path(dir, "plate-setup")
	dir.create(out_dir)

	ps <- read.delim(file.path(dir, "all-pvalues.txt"), row.names = 1)
	plotDataOnPlates(-log10(ps), paste0(out_dir, "/pvals--plate-view"), plates, skip = 0, ValueName = "-log10(p)", addSummaries = FALSE)   # 

	logfc <- read.delim(file.path(dir, "all-logfoldchange.txt"), row.names = 1)
	plotDataOnPlates(logfc, paste0(out_dir, "/logfc--plate-view"), plates, skip = 0, ValueName = "log fold change", symmetric = TRUE, addSummaries = FALSE)   # 

}

