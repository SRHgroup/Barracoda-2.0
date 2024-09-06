################################
###  COMMAND LINE ARGUMENTS  ###
################################

library(openxlsx)
library(tidyverse)

################################
###  COMMAND LINE ARGUMENTS  ###
################################

if(!interactive()) {cArgs <- commandArgs(TRUE)}
sample_map    <- cArgs[1]
sample_map_out    <- cArgs[2]
annotations <- cArgs[3]


################################
###    SAMPLE ID TABLE       ###
################################

# Read sample_map file
sid <- read.samplemap(sample_map)

# Save as .txt file
write.table(sid, file = sample_map_out, row.names = F, col.names = F, sep = "\t", quote = F)

# Make exp_ma
sid$akey <- sid$akey_num
sid$akey_num <- as.numeric(sub("^.*_(\\d+)$", "\\1", sid$akey_num)) 
exp_map <- lapply(split(sid, sid$experiment), function(x) setNames(x[,2], x[,1]))

# Read sample id table and checks for duplicated samples or A-keys (within the same experiment)
any.duplications(sid, type='sample', file=sample_map_out)
any.duplications(sid, type='akey', file=sample_map_out)

################################
###    ANNOTATION FILE       ###
################################

# Read annotations (this also checks for duplicsted barcodes)
anno_list <- read.annotation(annotations)
if(length(anno_list) == 1) {
    if(length(exp_map) > 1) anno_list <- rep(anno_list, length(exp_map))
    if(is.null(names(anno_list))) names(anno_list) <- names(exp_map)
}


################################
###    EXPERIMENTS           ###
################################

# Experiments in sample id table
experiments <- sort(names(exp_map))

# Check if the number of experiments in sample id table and annotation file is the same
if(length(experiments) != length(names(anno_list))){
  stop("[ERROR] Uneven number of experiments: ",
       paste("Sample id table file has", length(experiments), "experiments and the annotion file has", length(names(anno_list)), "experiments"))
}

# Check whether all experiments in sample id table is present in the annotation file (as sheet)
found_annotation <- experiments %in% names(anno_list)
if(!all(found_annotation)){
  stop("[ERROR] Annotation table not found for some experiment names: ",
       paste(experiments[!found_annotation], collapse = ", "))
}
  
# Check whether all experiments in the annotation file (sheets) is present in the sample id table
found_sampe_id_table <- names(anno_list) %in% experiments
if(!all(found_annotation)){
  stop("[ERROR] Experiment not found in sample id table: ",
       paste(experiments[!found_annotation], collapse = ", "))
}


################################
###    CHECK HLA             ###
################################

# Annotation file - check HLA header is correct
for (i in names(exp_map)) {
  my_anno <- anno_list[[i]]
  myHLA <- findMHCheader(colnames(my_anno), i)
}

