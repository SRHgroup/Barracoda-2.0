#
###### This script generates the following files:
# all-readcounts.All.txt
# all-readcounts.xlsx
# all-readcounts--with-clonality-red.All.txt
# all-readcounts--with-clonality-red.xlsx
# barplot--total-read-count.pdf
# barplot--total-read-count.png
# experiment_<nn>   (folders for each experiment), each containing:
#   barplot--total-read-count.pdf
#   barplot--total-read-count.png
#   barplot--mean-uniq-count--with-signif.pdf
#   barplot--mean-uniq-count--with-signif.png
#   barplot--fold-change--with-signif.pdf
#   barplot--fold-change--with-signif.png
#   mean-readcounts--with-clonality-red.All.txt
#   mean-readcounts--with-clonality-red.xlsx
#   fold_change.xlsx
#   fold_change  (a folder containing multiple txt files)
#

################################
###        PACKAGES          ###
################################

library(limma)
library(edgeR)
library(reshape2)
library(ggplot2)
library(scales)
library(openxlsx)

################################
###  COMMAND LINE ARGUMENTS  ###
################################

if(!interactive()) {
	cArgs <- commandArgs(TRUE)
}

file_in    <- cArgs[1]
sample_map <- cArgs[2]
dir_out    <- cArgs[3]

pepA.fa    <- cArgs[4]
pepB.fa    <- cArgs[5]

NA_length    <- cArgs[6]
NB_length    <- cArgs[7]

annotations <- cArgs[8]

sum_of_counts <- cArgs[9]

# set sum_of_counts to FALSE as default
if (is.na(sum_of_counts)==TRUE) {sum_of_counts=FALSE}


cat("file_in:", file_in, "\n")
cat("sample_map", sample_map, "\n")
cat("dir_out", dir_out, "\n")
cat("pepA.fa ", pepA.fa , "\n")
cat("pepB.fa", pepB.fa, "\n")

cat("NA_length", NA_length, "\n")
cat("NB_length", NB_length, "\n")
cat("annotations", annotations, "\n")

###############################
###   READ AND PREPARE DATA  ###
################################

## Read the mapped reads
library(data.table)
x <- fread(file_in)



## SAMPLE NAMES AND MULTIPLE EXPERIMENTS

sid <- read.samplemap(sample_map)

sid$akey <- sid$akey_num
sid$akey_num <- as.numeric(sub("^.*_(\\d+)$", "\\1", sid$akey_num)) 


# A vector used to map between sample names and barcode ids
# Values = sample names ("input"...), Names = sample barcode ids (A-Key...)
# sid_map <- setNames(sid[,2], sid[,1]) 

# A vector used to map between experiment ids and barcode ids
# Values = experiment
# sid_series <- setNames(sid[,3], sid[,1]) 

# A list of samples used in each experiment
# one list item per experiment
# each item is a vector of sample names, names of vector are a akey numbers (like sid_map)
exp_map <- lapply(split(sid, sid$experiment), function(x) setNames(x[,2], x[,1]))

# Keep unchanged akey
x$akey <- x$sample_id


## Convert to numeric (remove this section later!)
# Convert Sample ID to pure numeric
x$sample_id <- as.numeric(sub("^.*_(\\d+)$", "\\1", x$sample_id)) 

# Add sample name to data (based on akey number)
# x$sample <- sid_map[x$sample_id]


## PlOT THE NUMBER OF READS per sample key (all keys on the chip, not just those in the current experiment)
tab <- aggregate(x$akey, by = list(x$akey), FUN = length)

p1 <- plotReadsPerSample(tab, sid)
printGgplot(p1, "total-reads-per-key", dir_out, height = calcBarplotHeight(tab[,1]))


## A AND B OLIGOS

## do that in the main perl script instead !!! ###
# # Perhaps use these to check that all barcodes from anno are represented in the fasta files??
# pepA <- getBarcodeNames(pepA.fa, ".*_(A[0-9]+)\\s?$")
# pepB <- getBarcodeNames(pepB.fa, ".*_(B[0-9]+)\\s?$")

# Simplify oligo names so they simply become "A1", "B20", etc..
## Peptide A
x$pepA <- sub("^.*(A\\d+)$", "\\1", x$pepA)

## Peptide B
x$pepB <- sub("^.*(B\\d+)$", "\\1", x$pepB)

## AB barcode
x$Barcode <- ABpepFun(x$pepA, x$pepB)




## BARCODE ANNOTATIONS

anno_list <- read.annotation(annotations)

if(length(anno_list) == 1) {

  # if just one anno file but >1 experiment, they must all use the same anno file
  if(length(exp_map) > 1) anno_list <- rep(anno_list, length(exp_map))

  # if anno came from text file (no sheet name), adopt exp name from sample key map
  if(is.null(names(anno_list))) names(anno_list) <- names(exp_map)

  }


# Sort annoations tables and experiment map the same way. Not sure it matters.
experiments <- sort(names(exp_map))

# Check whether annotation tables exist for all experiments in the sample key map.
found_annotation <- experiments %in% names(anno_list)

if(!all(found_annotation))
   stop("[ERROR] Annotation table not found for some experiment names: ",
        paste(experiments[!found_annotation], collapse = ", "))

# Sort the same way
anno_list <- anno_list[experiments]
exp_map <- exp_map[experiments]



################################
###       ANALYSE DATA       ###
################################


###      N6 sequences        ###
## Combine the two N6 sequences to make the N12
#x$N12 <- makeN12(x, NA_length, NB_length)

x[ , N12 := ifelse(nchar(N6A) == NA_length & nchar(N6B) == NB_length, paste0(N6A, N6B), NA)]


###     TURN INTO MATRIX  (all barcodes and samples)              ###
mat_all <- readCountMatrix(x, dim.names = c("Barcode", "sample"))

###    CLONALITY REDUCTION (all barcodes)     ###
matUniq_all <- clonReducedMatrix(x, dim.names = c("Barcode", "sample"))

# SAVE ALL READ COUNTS (all)
WriteMultiSheetExcel(
	x = lapply(ListEachColumn(mat_all), AnnotateTables, data.frame(row.names = rownames(mat_all)), "barcode"),
	dir = dir_out, file = "unfiltered-readcounts", andTxt = "first")

# SAVE ALL READ COUNTS, WITH CLONALITY REDUCTION
WriteMultiSheetExcel(
	x = lapply(ListEachColumn(matUniq_all), AnnotateTables, data.frame(row.names = rownames(matUniq_all)), "barcode"),
	dir = dir_out, file = "unfiltered-readcounts--with-clonality-red", andTxt = "first")



### FOR EACH EXPERIMENT, SEPARATELY
#

for (i in names(exp_map)) {
	print(paste('experiment:',i))
	dir_out_tmp <- file.path(dir_out, paste0("experiment_", i))
	dir.create(dir_out_tmp)

	dir_out_plot_logfc <- file.path(dir_out_tmp, "graphs_logfc")
	dir.create(dir_out_plot_logfc)

	dir_out_plot_counts <- file.path(dir_out_tmp, "graphs_counts")
	dir.create(dir_out_plot_counts)
	
	# Annotations and samples for this experiment
	my_anno <- anno_list[[i]]
	my_samples <- exp_map[[i]]
	
	# Split by HLA
	myHLA <- findMHCheader(colnames(my_anno), i)
	anno_split <- c(list(all = my_anno), split(my_anno, my_anno[,myHLA]))

        # Check if there is more than one occurrence of each MHC type in myHLA.
        # There must be more than one occurrence of the MHC type in order to run the MHC specific analysis.

        MHC_table <- table(my_anno[,myHLA])
        less_than_one <- names(MHC_table[MHC_table<=1])
        included_MHCs <- c('all', names(MHC_table[MHC_table>1]))

        if (length(less_than_one) > 0) {
          message('#WARNING: Skipping the MHC specific analysis for ', less_than_one, '. Since this MHC was found only once.')
          anno_split <- anno_split[names(anno_split) %in% included_MHCs]
        }

	# List to collect the HLA-specific analyses
	list_of_tables_mod_combined <- list()

		
	for (h in names(anno_split)) {
	  print(h) # print HLA
	  # if(h=="all")
	  #   next

	  anno <- anno_split[[h]]
	  
	  # Assign again, as I change it during this loop..!
	  my_samples <- exp_map[[i]]
#	print("my_samples")  
#	print(my_samples)
  	# Which AB oligos used in this experiment?	
  	peps <- rownames(anno)
  	
#	print("peps")
#	print(peps)
#	print("x")
#	print(x)	
  	###     TURN INTO MATRIX  (only specified barcodes and samples)   ###
  	mat <- readCountMatrix(x, peps, my_samples, c("Barcode", "sample"))
#  	print("mat")
#	print(mat)
  	###    CLONALITY REDUCTION   ###
  	matUniq <- clonReducedMatrix(x, peps, my_samples, c("Barcode", "sample"))
  	
        samples <- sampleList(my_samples)


# 	# make sum of counts for duplicated experiments
# 	if (sum_of_counts==TRUE) {
# 	 print('sum_of_counts==TRUE')
# 	 print(samples) 
# 	 #i=0
#          for (s in samples) {
# 	   print(s)
#           # i=i+1
#            matUniq[,s[1]] <- as.matrix(rowSums(matUniq[,s]))
#            print(matUniq)
#            deselect <- s[s!= s[1]]
# 	   print(deselect)
#            matUniq <- matUniq[,!colnames(matUniq)==deselect]
#            print(matUniq)
#            #samples[i] <- s[1]
#          }
# 	i=0
# 	for (s in samples) {
# 	   print(s)	
# 	   i=i+1
#            mat[,s[1]] <- as.matrix(rowSums(mat[,s]))
#            print(mat)
#            deselect <- s[s!= s[1]]
#            print(deselect)
#            mat <- mat[,!colnames(mat)==deselect]
#            print(mat)
#            samples[i] <- s[1]
#          }
# 	print(class(my_samples))
#         print(my_samples)
#         print(samples)
# 	my_samples <- my_samples[names(my_samples) %in% samples]
# 	print("TRUE LOOP DONE")
# 	}
# 	print("my_samples new")
# 	print(my_samples)
# 	print("Samples")
# 	print(samples)
# 	
# 	print("mat")
# 	print(mat)
# 	print("matUniq")
# 	print(matUniq)
	
  	if( h == "all" ) {
  
  	  # SAVE ALL READ COUNTS  (only your samples and barcodes)
    	WriteMultiSheetExcel(
    	  x = lapply(ListEachColumn(mat), AnnotateTables, anno, "barcode"),
    	  dir = dir_out_tmp, file = "all-readcounts", andTxt = "first")
    	
    	# SAVE ALL READ COUNTS, WITH CLONALITY REDUCTION
    	WriteMultiSheetExcel(
    	  x = lapply(ListEachColumn(matUniq), AnnotateTables, anno, "barcode"),
    	  dir = dir_out_tmp, file = "all-readcounts--with-clonality-red", andTxt = "first")
    
    	# ###     PLOT SAMPLE TOTALS   ###
    	df2 <- prepareTotals(mat, matUniq, my_samples)
    	
    	# PLot totals
    	p2 <- plotTotals_Exp(df2) + ggtitle(paste("Experiment", i))
    	printGgplot(p2, "annotated-reads-per-sample", dir_out_tmp, height = calcBarplotHeight(my_samples) )
  
    }  	
  	
  	### Check for samples with zero read counts ------------------------------
  	akeys_with_more_than_zero_reads <- removeSamplesWithZeroReads(names(my_samples), matUniq)
  	my_samples <- my_samples[akeys_with_more_than_zero_reads]
  	
  	
  	###        REPLICATES        ###
  	matMean     <- ReplicateMeans(mat, my_samples)
  	matUniqMean <- ReplicateMeans(matUniq, my_samples)
  	
  	if( h == "all" ) {
  	  # Save for later for plotting together with HLA-specific pvalues
  	  matUniqMean_all <- matUniqMean
  	  
  	  # Save so we know later which samples have enough counts when we consider all HLAs
  	  my_samples_all <- my_samples
  	  
    	### SAVE MEAN READ COUNTS, WITH CLONALITY REDUCTION
    	WriteMultiSheetExcel(
    	  x = lapply(ListEachColumn(matMean), AnnotateTables, anno, "barcode"),
    	  dir = dir_out_tmp, file = "mean-readcounts", andTxt = "first")
    	WriteMultiSheetExcel(
    	  x = lapply(ListEachColumn(matUniqMean), AnnotateTables, anno, "barcode"),
    	  dir = dir_out_tmp, file = "mean-readcounts--with-clonality-red", andTxt = "first")

  	}  	
  	
  	###  FIND ENRICHED EPITOPES  ###-----------------------------------------------------------
  	samples <- sampleList(my_samples)

  	if(length(samples) == 0) {
  	  message("WARNING: Could not determine enriched barcodes for experiment '",
  	          i, "' MHC '", h,
  	          "' because all samples had zero counts.")
  	  next
  	}
  	
  	input <- names(my_samples)[my_samples == "input"]

  	if(length(input) == 0) {
  	  message("WARNING: Could not determine enriched barcodes for experiment '",
  	          i, "' MHC '", h, "' because all input samples had zero counts.")
  	  next
  	}
  
  	fc <- calcFoldChange(samples, input, matUniq)
	
  	if(ncol(fc$logFC) == 0) {
  	  message("WARNING: No samples left with enough reads. Skipping analysis of experiment '", i, "' for MHC '", h, "'")
  	  next
  	}

        if (ncol(fc$logFC)!=length(names(samples))) {
          # Remove samples not found in the fc dataframe from samples and my_samples 
          sample_differences=setdiff(names(samples),colnames(fc$logFC))
          samples <- samples[!names(samples) %in% sample_differences]
          my_samples <-  my_samples[! my_samples %in% sample_differences]
          message("WARNING: the following sample ", sample_differences, " have been skipped for MHC ", h,".")
        }
  	
  	fdr <- 0.001
  	enr <- findEnrichedBarcodes(fc$logFC, fc$FDR, fdr, matUniq)
  
  	
  	###    NORMALIZE READ COUNTS BY EDGER CALCULATED NORM.FACTORS (CORRECTING FOR LIBRARY SIZE) ###
    matUniqNorm <- normReadCounts(samples, my_samples, matUniq, fc, enr)
  	
    p3 <- plotNormReadCounts_xy(matUniqNorm, fdr)
    printGgplot(p3,
                paste0("scatterplot--exp", i, "--", h),
                dir_out_plot_logfc,
                height = 1+3*ceiling(length(samples)/2),
                width = 8)

    
    if( h == "all" ) {
      
    	###     PLOT READ COUNTS     ###
    	counts_df <- prepareCounts(matUniqMean, unique(my_samples), peps, enr, anno, i)
    	colors <- setNames(c("black", CbFCols()[1]), levels(counts_df$Enriched))
    	counts_p1 <- plotCounts(counts_df, colors)
    	printGgplot(counts_p1, paste0("barplot--exp", i, "--overall"), dir_out_plot_counts, width = 12)
    
    	counts_p2 <- plotCounts_HLA(counts_df, colors)
    	printGgplot(counts_p2, paste0("barplot--exp", i, "--overall--col-by-MHC"), dir_out_plot_counts, width = 12)
    
    }  	
  
  	###    SAVE TABLE AS EXCEL   ###
  	list_of_tables <- prepareTables(fc, my_samples, matUniq, c("Barcode", "sample"))
  	
  	list_of_tables_mod <- lapply(list_of_tables, function(x) {
  	  
  	  # add minus log 10 transformed p values
  	  x$`-log10(p)` <- - log10(x$p)
  
  	  # add a masked pvalue which is 1 when logfc is negative
  	  x$`masked_p (p = 1 if logFC < 0)` <- ifelse(x$log_fold_change < 0, 1, x$p)
  	  # and its minus log10 transofrm:
  	  x$`-log10(masked_p)` <- ifelse(x$log_fold_change < 0, 0, x$`-log10(p)`)
  	  
  	  # add edger normed counts:
  	  sample_i <- x$sample[1]
  	  norm_i <- subset(matUniqNorm, sample == sample_i)
  	  x$`count.normalised (edgeR)` <- norm_i[match(x$Barcode, norm_i$Barcode), ]$count
  	  x$`input.normalised (edgeR)` <- norm_i[match(x$Barcode, norm_i$Barcode), ]$input
  
  	  return(x)	
  	  
  	})
  	
  	list_of_tables_mod_combined[[h]] <- list_of_tables_mod
	
	}
	
	# Combine across HLAs
	#my_samples <- names(sampleList(exp_map[[i]]))
	my_samples <- names(sampleList(my_samples_all))  ## we stored this earlier in the h="all" HLA iteration loop

	# Pull out the tables where statistics was calculated on all barcodes, irrespective of HLA
	list_of_tables_overall     <- lapply(my_samples, function(s) list_of_tables_mod_combined[["all"]][[s]] ) 
	
	# Combine the tables where statistics was calculated within each group of HLA specific barcodes
	list_of_tables_hlaspecific <-
	  lapply(my_samples, function(s) {
	    by_hla <- lapply(names(anno_split), function(h) {
	      list_of_tables_mod_combined[[h]][[s]]
	    })
	    names(by_hla) <- names(anno_split)
	    by_hla[["all"]] <- NULL
	    join_hla <- do.call(rbind, by_hla)
	    ind <- match(rownames(my_anno), join_hla$Barcode)
	    join_hla <- join_hla[ind, ]
      
	    # fill out in case of NAs(occurs if an HLA did not have enough reads to calculate log fold change)
	    join_hla$Barcode <- rownames(my_anno)
	    join_hla$sample <- s
	    return(join_hla)
	  })

	# Preserve names of samples, so excel sheets will get the names
	names(list_of_tables_overall)     <- my_samples
	names(list_of_tables_hlaspecific) <- my_samples
	
	### SAVE FOLD CHANGE AND P VALUES
	WriteMultiSheetExcel(x = lapply(list_of_tables_overall, AnnotateTables, my_anno, "barcode", by.x = 1),
		dir = dir_out_tmp, file = "fold_change", andTxt = "all")

	### SAVE FOLD CHANGE AND P VALUES
	WriteMultiSheetExcel(x = lapply(list_of_tables_hlaspecific, AnnotateTables, my_anno, "barcode", by.x = 1),
	                     dir = dir_out_tmp, file = "fold_change_by_MHC", andTxt = "all")
	
	

	###     PLOT READ COUNTS AGAIN, BUT WITH HLA-SPECIFIC P-VALUES    ###
	enr_all <- findEnrichedBarcodes(
	  logfc = sapply(list_of_tables_hlaspecific, function(x) x$log_fold_change),
	  p = sapply(list_of_tables_hlaspecific, function(x) x$p),
	  fdr, matUniqMean_all)
	
	counts_df_all <- prepareCounts(matUniqMean_all, unique(my_samples), rownames(my_anno), enr_all, my_anno, i)
	
	colors <- setNames(c("black", CbFCols()[1]), levels(counts_df_all$Enriched))
	counts_p1 <- plotCounts(counts_df_all, colors)
	printGgplot(counts_p1, paste0("barplot--exp", i, "--p-by-MHC"), dir_out_plot_counts, width = 12)
	
	counts_p2 <- plotCounts_HLA(counts_df_all, colors)
	printGgplot(counts_p2, paste0("barplot--exp", i, "--p-by-MHC--col-by-MHC"), dir_out_plot_counts, width = 12)
	
	
		
}











