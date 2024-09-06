options(stringsAsFactors = FALSE, scipen = 1000)

#------- COMMAND LINE ARGUMENTS 


if(!interactive()) {
	cArgs <- commandArgs(TRUE)
} else {
	cArgs <- c("out")
}

base_dir     <- cArgs[1]


for (dir in list.dirs(base_dir, recursive = FALSE)) {

	# list of files with fold change data for each sample in experiment 1
	fc_files <- list.files(file.path(dir, "fold_change"), full.names = TRUE)

	if(length(fc_files) == 0) next
	
	# read each fold change table 
	dats <- lapply(fc_files, read.delim)

	# rename according to sample name (found in file name)
	names(dats) <- sub("\\..*$", "", basename(fc_files)) 

	# collect p values
	ps <- sapply(dats, function(i) setNames(i$p, i$barcode))
	write.table(ps,    file = file.path(dir, "all-pvalues.txt"),       sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE)

	# collect p values
	logfc <- sapply(dats, function(i) setNames(i$log_fold_change, i$barcode))
	write.table(logfc, file = file.path(dir, "all-logfoldchange.txt"), sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE)


}

