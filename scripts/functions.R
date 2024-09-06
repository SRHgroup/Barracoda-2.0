#--- this files is run ahead of all the r scripts of barracoda.
#--- consider making an r package instead (one day...)

library(tidyverse)

################################
###        FUNCTIONS         ###
################################

## read annotation file whether its excel or plain text

read.annotation <- function(f) {
  require(openxlsx)
  
  #--- find out whether its a txt file or an excel sheet
  filetype <- system(paste("file -b", f), intern = TRUE)
  
  isText  <- grepl("text", filetype, ignore.case = TRUE)
  
  isExcel <- grepl("Microsoft", filetype, ignore.case = TRUE) ||
    grepl("Excel", filetype, ignore.case = TRUE)
  
  if(!isText & !isExcel) stop("[ERROR] Annotation file does not seem to be either a text file or an excel sheet.")
  
  if(isText) {
    
    anno_tmp <- read.delim(file_w_last_newline(f), as.is = TRUE, row.names = NULL)
    
    if(any.duplications(anno_tmp[,1])) stop()
    
    anno <- read.delim(f, as.is = TRUE, row.names = 1)
    
    return(list(anno))
  }
  
  if(isExcel) {
    
    sheets <- getSheetNames(f)
    names(sheets) <- sheets
    
    anno <- lapply(sheets, function(i) {

      anno_tmp <- read.xlsx(f, sheet = i, rowNames = FALSE)
      
      if(any.duplications(anno_tmp[,1], filename = i)) stop()
      
      read.xlsx(f, sheet = i, rowNames = TRUE)
    
    })
    
  }
  return(anno)  
}

# Handle text files that are missing the trailing newline on last line
file_w_last_newline <- function(f) {

  pipe(paste("sed -e '$a\\'", f))
  
  }



#any.duplicated.barcodes <- function(barcodes, filename = "") {
#  
#  dups <- duplicated(barcodes)
#  
#  if(any(dups)) {
#    
#    message("ERROR: Barcodes cannot be repeated in annotation file ",
#            filename,
#            ". Found one or more duplicated barcodes (",
#            paste(unique(barcodes[dups]), collapse = ", "),
#            ")")
#    
#    return(TRUE)
#  
#  } else {
#      return(FALSE)
#  }
#
#}

### NEW add check for buplicated barcodes, samples and A-keys
any.duplications <- function(data, type='Barcodes', filename = "") {
  if (type=='Barcodes') {
    df <- data
  } 
  if (type=='sample') {
    df <- data %>% filter(sample!="input") %>% select(sample,experiment) %>% group_by(experiment)     
  } 
  if (type=='akey') {
    df <- data %>% filter(sample!="input") %>% select(akey,experiment) %>% group_by(experiment)     
  } 
  
  # find duplicates
  dups <- df %>% duplicated()
  
  if(any(dups)) {
    if (type=='Barcodes') {
      message("ERROR: Barcodes cannot be repeated in barcode annotation ", filename,
              ". Found one or more duplicated barcodes (", paste(unique(df[dups]), collapse = ", "),")")
    }
    if (type=='sample') {
      message("ERROR: Samples cannot be repeated in the sample idenfication table file ", filename,
              ". Found one or more duplicated entries (", paste(df[dups,]$sample, collapse = ", "),")")
    }
    if (type=='akey') {
      message("ERROR: A-keys cannot be repeated in the sample idenfication table file ", filename,
              ". Found one or more duplicated A-key (", paste(df[dups,]$akey, collapse = ", "),")")
    }
    return(TRUE)
  } else {
    return(FALSE)
  }
}


## read sample map file whether its excel or plain text
# only reads the first sheet of an excel file
# if only 2 columns, add experiment column and put the value 1 in each row
# give column names

read.samplemap <- function(f) {
  require(openxlsx)
  
  #--- find out whether its a txt file or an excel sheet
  filetype <- system(paste("file -b", f), intern = TRUE)
  
  isText  <- grepl("text", filetype, ignore.case = TRUE)
  
  isExcel <- grepl("Microsoft", filetype, ignore.case = TRUE) || grepl("Excel", filetype, ignore.case = TRUE)
  
  if(!isText & !isExcel) stop("[ERROR] Sample-key-map file does not seem to be either a text file or an excel sheet.")
  
  if(isText) sid <- read.delim(file_w_last_newline(f), header = FALSE, as.is = TRUE)

  if(isExcel) sid <- read.xlsx(f, colNames = FALSE)
  
  # If experiment-column is missing, all must be from same experiment
  if(ncol(sid) == 2) sid <- data.frame(sid, V3 = rep(1, nrow(sid)))
  
  colnames(sid) <- c("akey_num", "sample", "experiment")
  
  return(sid)
}


### Possible barcode names
#
# Get the names of the barcode A and B tag sequences, from the fasta files.
# Only keep the number at the end, eg.
#    sequence header:  >Pep_2OS_1_Oligo_A3 38 62
#    sequence id:      >Pep_2OS_1_Oligo_A3
#    keep this number: 3
getBarcodeNames <- function(file, pattern) {
  fasta_headers <- system(paste("grep '>'", file), intern = TRUE) # Get the fasta header lines
  fasta_ids <- strsplit(fasta_headers, " ")  # Remove the trailing numbers and spaces
  pep_names <- sapply(fasta_ids , function(i) {
    sub(pattern, "\\1", i[1], perl = TRUE) # Extract the name based on the pattern
  })
  return(pep_names)
}

## The function used to make a combined name for a pair of A- and B- barcodes.
ABpepFun <- function(a,b) {
  paste0(a,b)
}

## Combine A and B barcodes in all possible ways
makeABpeps <- function(A,B) {
  as.vector(t(outer(A, B, FUN = ABpepFun)))
}

makeN12 <- function(x, lengthA, lengthB) {
  perfect_N6 <- nchar(x$N6A) == lengthA & nchar(x$N6B) == lengthB
  N12 <- rep(NA, length(perfect_N6))
  N12[perfect_N6] <- apply(x[perfect_N6, c("N6A", "N6B")], 1, function(x) paste(x, collapse = ""))
  return(N12)
}

## Make a matrix of read counts per sample-barcode pair
# mat <- table( barcode = factor(x$barcode, levels = peps), sample = x$sample_id)
readCountMatrix <- function(x, peps, sid_map, dim.names) {
  if(missing(peps) | missing(sid_map)) {
    peps    <- sort(unique(x$Barcode))
    samples <- sort(unique(x$sample_id))
  } else {
    samples  <- names(sid_map)
  }
  
  index.list <- list(
    factor(x$Barcode, levels = peps),
    factor(x$sample_id, levels = samples)
  )
  names(index.list) <- dim.names
  m <- tapply(	X = vector(length = nrow(x)),
               INDEX = index.list,
               FUN = length)
  m[is.na(m)] <- 0
  return(m)
}

clonReducedMatrix <- function(x, peps, sid_map, dim.names) {
  if(missing(peps) | missing(sid_map)) {
    peps    <- sort(unique(x$Barcode))
    samples <- sort(unique(x$sample_id))
  } else {
    samples  <- names(sid_map)
  }
  
  ## Split by sample, to make a list
  xl <- split(x, factor(x$sample_id, levels = samples))
  
  ## Clonality reduction
  m <- sapply(xl, function(sample) {
    sapply(split(sample, factor(sample$Barcode, levels = peps)), function(pep) {
      length(unique(pep$N12[!is.na(pep$N12)]))
    })
  })
  names(dimnames(m)) <- dim.names
  return(m)	
}

ReplicateMeans <- function(x, xmap) {
  # x is a matrix of rows = barcodes, cols = samples
  # xmap is a vector with
  #     values = meaningful sample names, which may be duplicated
  #     names of the values = colnames of x
  
  # first remove any columns of x not mentioned in xmap
  x <- x[ , names(xmap), drop = FALSE]
    
  if (any(duplicated(xmap))) {
    xMean <- sapply(split(names(xmap), factor(xmap, levels = unique(xmap))), function(s) {
      if (length(s) > 1) { return(rowMeans(x[,as.character(s)])) }
      else { return(x[ ,as.character(s)])	}
    })
    names(dimnames(xMean)) <- names(dimnames(x))
  } else {
    xMean <- x
    colnames(xMean) <- xmap[colnames(x)]
  }
  names(dimnames(xMean)) <- names(dimnames(x))
  return(xMean)
}

## why does this return a list rather than a vector? -AE
sampleList <- function(sid_map) {
  s <- split(names(sid_map), factor(sid_map, levels = unique(sid_map)))
  s$input <- NULL
  return(s)
}


removeSamplesWithZeroReads <- function(samples, mat) {
  ## Check for samples with zero read counts ------------------------------
  ##  And remove them. Works on both vectors and lists of vectors.
  
  if(is.list(samples)) {
    samples_new <-lapply(samples, removeSamplesWithZeroReads, mat)
    return(samples_new[!sapply(samples_new, is.null)] )
  }
  
    libsizes <- colSums(mat[ , samples, drop = FALSE])
    
    if(all(libsizes == 0)) {
      message("WARNING: the following akeys have been skipped because they had zero reads: ",
              paste(samples, collapse = ", "))
      return(NULL)
    }
    
    if(any(libsizes == 0)) {
      message("WARNING: the following akeys have been skipped because they had zero reads: ",
              paste(samples[libsizes == 0], collapse = ", "))
      return(samples[libsizes > 0])
      
    }
    
    return(samples)
    
}

calcFoldChange <- function(sample.list, input.samples, mUniq){
  # source("http://bioconductor.org/biocLite.R")
  # biocLite("edgeR")
  require(edgeR)
	#print("cal log fold change input")
	#print(sample.list)
	#print(input.samples)
	#print(mUniq)
  ## Run edgeR analysis ---------------------------------------------------
  edge <- lapply(names(sample.list), function(i) {
    s <- sample.list[[i]]

    y  <- DGEList(counts = mUniq[, c(s, input.samples)], 
                  group = c(rep(2, length(s)), rep(1, length(input.samples)))) 
    #print("y")
    #print(y)
    # What to do if there are samples with zero reads?
    # currently causes the function (and program) to crash, so we should notice these and ignore them
    # calc the library size first, and give as argument to calcNormFactors
    # how it is cal inside calcNormFactors if not provided:  lib.size <- colSums(x)
    # do this ourselves, and then make sure to skip any that have lib.size = 0
    yn  <- calcNormFactors(y, method = "TMM")
    if (any(is.na(yn$samples$norm.factors))) {
      ## Too few counts in this sample, could not calculate norm.factors
      message("WARNING: Not enough reads in sample ", i, " to determine Enriched barcodes.")
      return(data.frame(logFC = rep(NA, nrow(mUniq)), FDR = rep(NA, nrow(mUniq))))
    }
    ## AE: dispersion = 0.1 is based on fitting of undetected peptides in "mixed" data
    etn <- exactTest(yn, dispersion = 0.1)  
    ttn <- as.data.frame(topTags(etn, n = nrow(mUniq), sort.by = "none"))
    
    # save norm.factors for later
    norm.factors <- yn$samples
    norm.factors$analysis <- i
    norm.factors$akey <- row.names(norm.factors)
    row.names(norm.factors) <- NULL
    
    return( list(logFC = ttn[ ,"logFC"],
                 FDR   = ttn[ ,"FDR"],
                 norm = norm.factors)
    )
  })
  

  results <- lapply(c(logFC = "logFC", FDR = "FDR"), function(i) {
    j <- do.call(cbind, lapply(edge, function(x) x[[i]]))
    colnames(j) <- names(sample.list)
    j <- j[,!apply(j, 2, function(x) all(is.na(x))), drop = FALSE]
    rownames(j) <- rownames(mUniq)
    names(dimnames(j)) <- names(dimnames(mUniq)) 
    return(j)
  })
  
  results$norm.factors <- do.call(rbind, lapply(edge, function(x) x[["norm"]]))
  
  
  return(results)
}

findEnrichedBarcodes <- function(logfc, p, fdr, mUniq) {
  
  enriched <- p < fdr & logfc > 0
  enriched <- apply(enriched, 2, as.numeric)
  rownames(enriched) <- rownames(mUniq)
  names(dimnames(enriched)) <- names(dimnames(mUniq)) 
  
  return(enriched)
}


normReadCounts <- function(sample.list, akey_names, mUniq, edgeR_output, enriched_barcodes) {
  
  result_list <- lapply(names(sample.list), function(i) {
    s <- sample.list[[i]]

    # get norm factors for this analysis    
    norm.facs <- subset(edgeR_output$norm.factors, analysis == i)
    
    # normalize read counts
    mNorm <- sapply(seq_along(norm.facs$akey), function(j) {
      akey <- norm.facs$akey[j]
      normf <- norm.facs$norm.factors[j]
      counts_norm <- mUniq[,akey] / ( sum(mUniq[,akey]) * normf )      
      return(counts_norm)
    })
    colnames(mNorm) <- norm.facs$akey
    
    # calc mean of replicates
    m <- ReplicateMeans(mNorm, akey_names[colnames(mNorm)])  
    colnames(m)[colnames(m) == i] <- "count"
    
    res <- data.frame(Barcode = rownames(m), sample = i, m)
    rownames(res) <- NULL
    return(res)
  })
  
  result_df <- do.call(rbind, result_list)
  
  logfc <- melt(edgeR_output$logFC, value.name = "logFC")
  fdr <- melt(edgeR_output$FDR, value.name = "FDR")
  enriched <- melt(enriched_barcodes, value.name = "Enriched")
  
  result_df <- merge(result_df, logfc)
  result_df <- merge(result_df, fdr)
  result_df <- merge(result_df, enriched)
  
  return(result_df)
}

plotNormReadCounts_xy <- function(m, fdr) {
  require(squash)
  require(ggplot2)
  m$Enriched <- factor(c("Not","Enriched")[1+m$Enriched], levels = c("Not","Enriched"))
  
  minFc <- min(-3, pretty(m[,"logFC"])[1] )
  maxFc <- max( 3, rev(pretty(m[,"logFC"]))[1]      )

  col_labs <- c("< -2","< -1","-1 to 1","> 1","> 2")
  m$logFC_factor <- cut(m[,"logFC"],
     breaks = c(minFc,-2,-1,1,2,maxFc),
     labels = col_labs)
  m$logFC_factor <- factor(m$logFC_factor, levels = rev(col_labs))
  
  p <- ggplot(m,
         aes(x = input,
             y = count,
             color = logFC_factor)) +
    geom_point(shape = 16, size = 1, alpha = 0.8, stroke = 0.5) +
    scale_color_manual("Log fold change", values = rev(setNames(blueorange(5), col_labs)), drop = FALSE) +
    scale_x_log10() + scale_y_log10() + ylab("Sample") +
    facet_wrap(~sample, scales = "free", ncol = 2) +
    # theme_classic() +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
	  text = element_text(size=6), legend.text=element_text(size=5), legend.title = element_text(size=5), legend.key.size = unit(0.3, 'cm')
	  ) +
  
    # theme(panel.spacing=grid::unit(2, "lines")) +
    geom_abline(slope = 1)
  
  m_signif <- subset(m, Enriched == "Enriched")
  if(nrow(m_signif) > 0) {
    m_signif$Enriched <- paste("FDR <", fdr, "\nand logFC > 0")
    p <- p + geom_point(data = m_signif,
                        aes(x = input,
                            y = count,
                            shape = Enriched),
                         color = 1, size = 1, stroke = 0.5) +
      scale_shape_manual("", values = c(1,1))
  }
  
  p <- p + ggtitle("Fraction of reads per barcode (normalized by edgeR)")
  
  }



MakeBins <- function(
  x,           # The AB barcodes for which you need to assign bins.
  ABpeps,      # The unique AB-combinations in use. A1B1, A1B2, etc... (eg. factor levels)
  ### ABpeps is currently ignored!!! -Aron
  sep = "-",
  splitB = 2   # Split B's into 2 bins So if there are 24 B's, make the bins from 1-12 and 13-24
) {
  x <- as.character(x)
  x.A <- as.integer(sub("^A(\\d+)B(\\d+)$", "\\1", x, perl = TRUE))
  x.B <- as.integer(sub("^A(\\d+)B(\\d+)$", "\\2", x, perl = TRUE))
  allA <- sort(unique(x.A))
  allB <- sort(unique(x.B))
  bins <- expand.grid(B = allB, A = allA)[, 2:1]
  rownames(bins) <- paste0("A", bins$A, "B", bins$B)
  bins$Bgroup <- cut(bins$B, breaks = splitB, labels = FALSE)
  
  Bsplit <- split(allB, cut(allB, splitB))
  from <- sapply(Bsplit, min)
  to <- sapply(Bsplit, max)
  
  bins$label <- paste0('A', bins$A, 'B', from[bins$Bgroup], sep, to[bins$Bgroup])
  
  out <- bins[x, "label"]
  factor(out, levels = unique(out))
}


CbFCols <- ColorblindFriendlyColors <- function(n = 7) {
  repeats <- ceiling(n / 7)
  colors <- c(orange = rgb( 230, 159,   0,  max=255),
              sky_blue = rgb( 86,  180, 233,  max=255),
              bluish_green = rgb( 0,   158, 155,  max=255),
              yellow = rgb( 240, 228,  66,  max=255),
              blue = rgb( 0,   114, 178,  max=255),
              vermillion = rgb( 213, 94,    0,  max=255),
              reddish_purple = rgb( 204, 121, 167,  max=255))
  return( rep(colors, repeats)[1:n] )
}

completeFactors <- function(x, levels) {
  # Add the empty factors, for plotting nicely
  addRows <- x[NULL,]
  for (i in levels) {
    # If there is no data for this barcode, and a row with NA
    if(nrow(subset(x, Barcode == i)) == 0) {
      addRows <- rbind(addRows, data.frame(i, x[1,2], rep(NA, ncol(x)-2)), deparse.level = 2)
    }
  }
  colnames(addRows) <- colnames(x)
  return( rbind(x, addRows) )
}



### FUNCTIONS ###
MarkReplicates <- function(labels) {
  i <- 1
  group <- rep(i, length(labels))
  repeat {
    if (any(duplicated(labels[group == i]))) {
      group[group == i][duplicated(labels[group == i])] <- i+1
      i <- i+1
    } else {
      break
    }
  }
  names(group) <- names(labels)
  return(group)
}




prepareCounts <- function(m, sample_levels, pep_levels, p, a, exp_name) {
  require(reshape2)
  
  # order by input counts
  barcode_index <- rownames(m)[order(m[,"input"])]
  
  tmp <- melt(m, value.name = "Count", as.is = TRUE)
  df <- merge(subset(tmp, sample == "input", select = c("Barcode", "Count")),
              subset(tmp, sample != "input"), by = "Barcode", suffixes = c("_input", ""))
  
  df_norm <- unsplit(lapply(split(df, f = df$sample), function(x) {
    x$Count_input_norm <- sum(x$Count) * x$Count_input/sum(x$Count_input)
    x$Frac_input <- x$Count_input/sum(x$Count_input)
    x$Frac_sample <- x$Count/sum(x$Count)
    return(x)
  }), f = df$sample)
    
  df <- df_norm
    
  sample_levels <- sample_levels[sample_levels != "input"]
  df$sample <- factor(df$sample, levels = sample_levels)
  
  df$Barcode <- factor(as.character(df$Barcode), levels = barcode_index)
  
  # df$Count <- ifelse(df$Count > 0 & df$Count <= 1, 1.1, df$Count)
  
  # Add enrichment
  df_p <- melt(p, value.name = "Enriched", as.is = TRUE)
  df <- merge(df, df_p, all = TRUE)
  
  # # Add HLA info
  hla_index <- findMHCheader(colnames(a), exp_name)
  
  hla <- a[ , hla_index, drop = FALSE]
  colnames(hla) <- "MHC_restriction"
  df <- merge(df, hla, by.x = "Barcode", by.y = "row.names", all = TRUE)
  
  
  # Simplify HLA names
  df$MHC_restriction <- sub("HLA", "", df$MHC_restriction)
  df$MHC_restriction <- sub("-", "", df$MHC_restriction)
  df$MHC_restriction <- sub(":", "", df$MHC_restriction)
  df$MHC_restriction <- sub("\\*", "", df$MHC_restriction)
  
  df$Enriched <- ifelse(is.na(df$Enriched), 0, df$Enriched)
  df$Enriched <-  c("p >= 0.001", "p < 0.001")[match(df$Enriched, c(0, 1))]
  df$Enriched <- factor(df$Enriched, levels = c("p >= 0.001", "p < 0.001"))
  
  return(df)
}


findMHCheader <- function(headers, filename = "") {
  
  # patterns to find allele header
  mhc_patterns <- c("MHC", "allele", "HLA", "H2", "H-2")
  
  matching_headers <- unique(
    unlist(
      sapply(mhc_patterns,
             grep,
             headers,
             value = TRUE,
             ignore.case = TRUE
             )
        )
    )
  
  if(length(matching_headers)  > 1) {
    
    hla_index <- matching_headers[1]
    
    message("WARNING: found more than one possible MHC column in annotation file ",
            filename,
            "(",
            paste(matching_headers, collapse = ", "),
            "). Using the first one (",
            hla_index,
            ")."
            )
  
  } else if(length(matching_headers) == 0) {
    
    stop("ERROR: Could not find HLA info in annotation file ",
            filename,
            ". All annotation files must include a column with HLA information. ", 
            "Make sure the column header includes one of the following words so that I can find it: ",
            paste(mhc_patterns, collapse = ", "),
            ".")
    
  } else if(length(matching_headers) == 1) {
    hla_index <- matching_headers[1] 
  }
  
  return(hla_index)
  
}


plotCounts <- function(x, cols) {
  require(ggplot2)
  x$fill_legend <- "Input"
  p <- ggplot(x,
              aes(x = Barcode, y = Count))
  p <- p + geom_bar(aes(y = Count_input, fill = fill_legend), stat = "identity")
  p <- p + geom_point(aes(shape = Enriched), size = 0.5, stroke = 0.5)
  # p <- p + theme_bw()
  p <- p + theme(axis.ticks = element_blank(),
                 axis.text.x = element_blank(),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 text = element_text(size=6), legend.text=element_text(size=5), legend.title = element_text(size=5), legend.key.size = unit(0.3, 'cm')
  ) 
  maxCount <- max(x[,c("Count", "Count_input")],na.rm=TRUE)
  p <- p + scale_y_log10(limits=c(1,10^(log10(maxCount)*1.05)), expand = c(0, 0))
  # p <- p + scale_y_log10(limits=c(1,10^1.05+max(x[,c("Count", "Count_input")],na.rm=TRUE)), expand = c(0, 0))
  # p <- p + scale_y_log10(expand = c(0, 0))
  p <- p + scale_color_manual(values = cols)
  p <- p + scale_shape_manual(values = setNames(c(1,8), nm = levels(x$Enriched)))
  p <- p + scale_fill_manual(name = "Baseline", values = "grey80")
  p <- p + ggtitle("Read count (with clonality reduction)")
  p <- p + facet_grid(sample ~ ., scales = "free")
  return(p)
}

plotCounts_HLA <- function(x, cols) {
  require(ggplot2)
  x$fill_legend <- "Input"
  x$Barcode <- factor(x$Barcode, levels = unique(x$Barcode[order(x$MHC_restriction, x$Count_input)]))
  p <- ggplot(x,
              aes(x = Barcode, y = Count))
  p <- p + geom_bar(aes(x = Barcode, y = Count_input, fill = MHC_restriction), stat = "identity")
  p <- p + geom_point(aes(shape = Enriched), size=0.5, stroke = 0.5)
  # p <- p + theme_bw()
  p <- p + theme(axis.ticks = element_blank(),
                 axis.text.x = element_blank(),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 text = element_text(size=6), legend.text=element_text(size=5), legend.title = element_text(size=5), legend.key.size = unit(0.3, 'cm')
  ) 
  
  maxCount <- max(x[,c("Count", "Count_input")],na.rm=TRUE)
  p <- p + scale_y_log10(limits=c(1,10^(log10(maxCount)*1.05)), expand = c(0, 0))
  #p <- p + scale_y_log10(limits=c(1,10^1.08+max(x[,c("Count", "Count_input")],na.rm=TRUE)), expand = c(0, 0))
  # p <- p + scale_y_log10(expand = c(0, 0))
  p <- p + scale_fill_manual(name = "Baseline", values = unname(ColorblindFriendlyColors(length(unique(x$MHC_restriction)))))
  p <- p + scale_color_manual(values = cols)
  p <- p + scale_shape_manual(values = setNames(c(1,8), nm = levels(x$Enriched)))
  # p <- p + scale_fill_manual("", values = "grey80")
  p <- p + ggtitle("Read count (with clonality reduction)")
  p <- p + facet_grid(sample ~ ., scales = "free")
  return(p)
}


# plotCounts_HLA <- function(x, cols) {
#   require(ggplot2)
#   x$fill_legend <- "Input"
#   p <- ggplot(x,
#               aes(x = Barcode, y = Count))
#   p <- p + geom_bar(aes(y = Count_input, fill = fill_legend), stat = "identity")
#   p <- p + geom_point(aes(color = Enriched))
#   p <- p + theme(axis.ticks = element_blank(),
#                  axis.text.x = element_blank(),
#                  panel.grid.major = element_blank(),
#                  panel.grid.minor = element_blank()) 
#   p <- p + scale_y_log10(limits=c(1,10^1.05+max(x[,c("Count", "Count_input")],na.rm=TRUE)), expand = c(0, 0))
#   # p <- p + scale_y_log10(expand = c(0, 0))
#   p <- p + scale_color_manual(values = cols)
#   p <- p + scale_fill_manual("", values = "grey80")
#   p <- p + ggtitle("Read count (with clonality reduction)")
#   p <- p + facet_grid(sample ~ HLA, scales = "free")
#   return(p)
# }

# plotCounts_HLA <- function(x, cols) {
#   require(ggplot2)
# 
#   hlas <- sub("HLA", "", x$HLA)
#   hla_levels <- unique(hlas)
#   hla_loci <- sub(".*([ABC]).*", "\\1",  hla_levels)
#   hla_number <- sub(".*?([0-9]+).*", "\\1",  hla_levels)
#   hla_levels <- hla_levels[order(hla_loci, hla_number)]
#   
#   p <- ggplot(x,
#               aes(x = Barcode, y = Count))
#   p <- p + geom_bar(stat = "identity")
#   p <- p + theme(axis.ticks = element_blank(),
#                  axis.text.x = element_blank(),
#                  panel.grid.major = element_blank(),
#                  panel.grid.minor = element_blank()) 
#   p <- p + scale_y_log10(limits=c(1,10^1.05+max(x[,"Count"],na.rm=TRUE)), expand = c(0, 0))
#   p <- p + ggtitle("Read count (with clonality reduction)")
#   p <- p + facet_grid(sample ~ HLA, scales = "free", space = "free_x")
#   return(p)
# }


### PREPARE FOR PLOTTING ###
prepareFoldChange <- function(m, sample_levels, pep_levels, p) {
  require(reshape2)
  
  df <- melt(m, value.name = "logFC", as.is = TRUE)
  df$sample <- factor(df$sample, levels = sample_levels)
  df$Barcode <- factor(as.character(df$Barcode), levels = pep_levels)
  
  df <- completeFactors(df, pep_levels)
  df$logFC <- ifelse(is.na(df$logFC), 0, df$logFC)
  
  df$barcode_bin <- MakeBins(df$Barcode, pep_levels, splitB = 2)
  
  # Add enrichment
  df_p <- melt(p, value.name = "Enriched", as.is = TRUE)
  df <- merge(df, df_p, all = TRUE)
  
  df$Enriched <- ifelse(is.na(df$Enriched), 0, df$Enriched)
  df$Enriched <-  c("p >= 0.01", "p < 0.01")[match(df$Enriched, c(0, 1))]
  df$Enriched <- factor(df$Enriched, levels = c("p >= 0.01", "p < 0.01"))
  
  return(df)
}

### MAKE GGPLOT ###
plotFoldChange <- function(x, cols) {
  require(ggplot2)
  
  p <- ggplot(x, aes(x = Barcode, y = logFC, fill = Enriched))
  p <- p + geom_bar(stat = "identity", position = "dodge")
  p <- p + theme(axis.ticks = element_blank(),
                 axis.text.x = element_blank(),
                 strip.text.x = element_text(size = 7, angle = 90),
                 legend.direction = "vertical",
                 text = element_text(size=6), legend.text=element_text(size=5), legend.title = element_text(size=5), legend.key.size = unit(0.3, 'cm')
  ) 
  p <- p + scale_fill_manual(values = cols) # 
  p <- p + ggtitle("Log fold change")
  p <- p + facet_grid(sample ~ barcode_bin, scales = "free")
  return(p)
}


### PREPARE FOR PLOTTING ###
prepareTotals <- function(m, mUniq, sid_map) {
  df <- rbind(
    data.frame(
      Sample = MarkIfReplicate(sid_map[colnames(m)]),
      Reads = colSums(m),
      clonal = "Before"
    ),
    data.frame(
      Sample = MarkIfReplicate(sid_map[colnames(mUniq)]),
      Reads = colSums(mUniq),
      clonal = "After"
    )
  )
  
  df$clonal <- factor(df$clonal, levels = c("After", "Before"))
  # df$Sample <- factor(df$Sample, levels = unique(sid_map))
  # df$Replicate <- factor(df$Replicate, levels = 1:max(df$Replicate))
  # df$Replicate <- as.character(df$Replicate)
  return(df)
}

MarkIfReplicate <- function(x) {
  unsplit( lapply( split(x, f = x), function(y) {
    if(length(y) > 1) paste(y, seq_along(y), sep = ".") else y
  }), f = x)
}

### MAKE GGPLOT TO PLOT READS PER SAMPLE - FOR ALL THE KEYS ###
plotReadsPerSample <- function(tab, sid) {
  require(scales)
  require(ggplot2)
  
  df <- merge(tab, sid, by.x = "Group.1", by.y = "akey", all = TRUE)
  colnames(df) <- c("Key", "Number.of.reads", "akey_num", "sample", "Experiment")
  
  # df <- data.frame(Key = names(tab), 'Number.of.reads' = as.numeric(tab), stringsAsFactors = FALSE)
  df$Key[df$Key == ""] <- "No key"
  df$Key[df$Key == "tie"] <- "Tie"
  
  # key_numbers <- as.numeric(sub(".*_(\\d+)$", "\\1", df$Key))
  df$Key <- sub(".*_(\\d+)$", "\\1", df$Key)
  df$Key <- factor(df$Key, levels = unique(df$Key[order(as.numeric(df$Key))]))
  
  # ensure experiment is a factor
  df$Experiment <- factor(df$Experiment)
  
  nudge <- 0.01 * max(df$Number.of.reads)
  ncols <- length(unique(df$Experiment))
  col <- ColorblindFriendlyColors(ncols)
  nkeys <- length(unique(df$Key))
  
  g <- ggplot(data = df,
              aes(x = Key,
                  y = Number.of.reads,
                  label = comma(Number.of.reads),
                  fill = Experiment)) +
    geom_bar(stat = "identity", position = "dodge") +
    scale_x_discrete(limits = rev(levels(df$Key))) +
    scale_y_continuous(labels = comma) +
    scale_fill_manual(values = unname(col), na.value = "grey") +
    coord_flip() +
    geom_text(vjust = 0.5, hjust = 1, size = 1.5, nudge_y = -nudge) +
    geom_text(aes(x = Key, label = sample, y = 0),
              vjust = 0.5, hjust = 0, size = 1.5, nudge_y = nudge) +
    theme(plot.margin = grid::unit(c(1,1,1,1), "cm"),
	  text = element_text(size=6), legend.text=element_text(size=5), legend.title = element_text(size=5), legend.key.size = unit(0.3, 'cm')
	  ) +
    guides(fill = guide_legend(nrow = min(nkeys, ncols)))
  
}

### MAKE GGPLOT ###
plotTotals_Exp <- function(x) {
  require(scales)
  require(ggplot2)
  require(reshape2)
  
  x$Sample <- as.character(x$Sample)
  mdf <- melt(x, id.var = c("Sample","clonal"))
  x <- dcast(mdf, Sample+clonal~variable, drop = FALSE, fill = NA)
  
  # x$Replicate <- factor(x$Replicate, levels = rev(sort(as.numeric(unique(x$Replicate)))))
  
  p2 <- ggplot(x, aes(x = Sample, y = Reads, fill = clonal)) +
    # geom_bar(position = "dodge", stat = "identity", fill = CbFCols(1), colour = "black") +
    geom_bar(stat = "identity", position = "dodge") +
    # facet_grid(~clonal, as.table = FALSE) +
    scale_y_continuous(labels = comma) +
    scale_x_discrete(limits = rev(levels(factor(x$Sample)))) +
    scale_fill_discrete("Clonality reduction", breaks = c("Before", "After")) +
    xlab("Sample") +
    theme(axis.ticks = element_blank(),
          strip.text.x = element_text(size=5),
	  text = element_text(size=6), legend.text=element_text(size=5), legend.title = element_text(size=5), legend.key.size = unit(0.3, 'cm')
	  ) +
    coord_flip()
  return(p2)
}

calcBarplotHeight <- function(n, base = 4, div = 6) {
  if(length(n) > 1) {
    n <- length(n)
  }
  floor(n/div) + base
}


### PREPARE TABLES ###
prepareTables <- function(fc, sid_map, mUniq, dim.names) {
  # Table of log fold change and p values
  table_pvals <- merge(
    melt(fc$logFC, value.name = "log_fold_change", as.is = TRUE),
    melt(fc$FDR, value.name = "p", as.is = TRUE)
  )
  
  # Table of read counts per barcode per sample
  table_counts <- melt(mUniq, value.name = "count", as.is = TRUE)
  
  # Specify the replicate number (1, 2, 3, ... etc)
  table_counts$rep <- MarkReplicates(sid_map)[table_counts$sample]
  
  # Change sample id to common sample name
  table_counts$sample <- sid_map[table_counts$sample]
  
  # Split by sample
  list_of_tables <- split(table_counts, table_counts$sample)
  
  # Reshape from long melted format to wide format (1 column per replicate, instead of 1 row)
  list_of_tables <- lapply(list_of_tables, function(i) {
    reshape(i, idvar = dim.names, timevar = "rep", direction = "wide")
  })
  
  # Table with just read counts from inputs
  table_input <- list_of_tables[["input"]]
  colnames(table_input) <- gsub("^count", "input", colnames(table_input))
  table_input$sample <- NULL
  
  # Remove the input from the data
  list_of_tables$input <- NULL
  
  # For each sample, merge with input counts, log fold change and p value
  list_of_tables <- lapply(list_of_tables, function(i) {
    m <- merge(i, table_input, sort = FALSE)  # input
    merge(m, subset(table_pvals, sample == i[1,"sample"]), sort = FALSE, all = TRUE) # log Fold Change and p values
  })
  
  correct.order <- order(match(names(list_of_tables), sid_map))
  list_of_tables <- list_of_tables[correct.order]
  
  return(list_of_tables)
}


# # Write x to an excel file. And optionally also to a txt file
# WriteExcel <- function(x, txtFile, append = FALSE, ...) {
# 	require(xlsx) #load the package
# 	write.xlsx(x = x, row.names = FALSE, append = append,  ...)
# 	if(!is.na(txtFile)) {
# 		headers <- !file.exists(txtFile)
# 		write.table(x, file = txtFile, sep = "\t", row.names = FALSE, col.names = headers, quote = FALSE, append = !headers)
# 	}	
# }

removeIfFound <- function(f) {
  if(file.exists(f)) {
    file.remove(f)
  }
}

# WriteMultiSheetExcel <- function(x, dir, file, andTxt = c("first", "all", "append", "none")[1], ...) {
# 	file <- file.path(dir, file)
# 	if(andTxt == "append") {
# 		txtFile <- rep(paste0(file, ".txt"), length(x))
# 		removeIfFound(txtFile[1])  # Make sure it doesn't exist already
# 	} else if(andTxt == "first") {
# 		txtFile <- paste0(file, ".", names(x), ".txt")
# 		txtFile[-1] <- NA
# 		removeIfFound(txtFile[1])  # Make sure it doesn't exist already
# 	} else if(andTxt == "all") {
# 		dir.create(file)
# 		txtFile <- paste0(names(x), ".txt")
# 		txtFile <- file.path(file, txtFile)
# 	} else {
# 		txtFile <- rep(NA, length(x))
# 	}
# 	names(txtFile) <- names(x)

# 	file <- paste0(file, ".xlsx")
# 	removeIfFound(file)

# 	for (i in names(x)) {
# 		WriteExcel(x = x[[i]], file = file, append = TRUE, sheetName = i, txtFile = txtFile[i], ...)
# 		# write table
# 	}
# }

WriteMultiSheetExcel <- function(x, dir, file, andTxt = c("first", "all", "append", "none")[1], ...) {
  require(openxlsx)
  file <- file.path(dir, file)
  
  # Determine how to print the txt files
  #
  # First:  Only make a txt file from the first element of the list x
  # All:    Make a new txt file for each element of the list x
  # Append: Append all elements of the list x to the same txt file
  # None:   Do not make any txt files
  #
  if(andTxt == "append") {
    txtFile <- rep(paste0(file, ".txt"), length(x))
    removeIfFound(txtFile[1])  # Make sure it doesn't exist already
  } else if(andTxt == "first") {
    txtFile <- paste0(file, ".", names(x), ".txt")
    txtFile[-1] <- NA
    removeIfFound(txtFile[1])  # Make sure it doesn't exist already
  } else if(andTxt == "all") {
    dir.create(file)
    txtFile <- paste0(names(x), ".txt")
    txtFile <- file.path(file, txtFile)
  } else {
    txtFile <- rep(NA, length(x))
  }
  names(txtFile) <- names(x)
  
  file <- paste0(file, ".xlsx")
  removeIfFound(file)

  # Write as txt as well
  for (i in names(x)) {
    if(!is.na(txtFile[i])) {
      headers <- !file.exists(txtFile[i])
      write.table(x[[i]], file = txtFile[i], sep = "\t", row.names = FALSE, col.names = headers, quote = FALSE, append = !headers)
    }
  }
  
  # don't write excel workbook if no. of samples  is too high (write.xlsx will crash with too many sheets)
  if(length(x) <= 1000) {
    write.xlsx(x, file = file, row.names = FALSE)
  }
  

}

# Take a dataframe and make a list where each element is a dataframe of ncol = 1 for each column.
# And the first element of the list is the whole table.
ListEachColumn <- function(x) {
  l <- c(list(x), lapply(colnames(x), function(i) {
    x[,i,drop=FALSE]
  }))
  names(l) <- c("All", colnames(x))
  # l <- lapply(l, function(i) {
  # rownames(i) <- rownames(x)
  # i
  # })
  return(l)
}

# Merge x with a dataframe of annotations. Default is to merge by rownames.
AnnotateTables <- function(x, a, rowHeader, by.x = "row.names", by.anno = "row.names", ...) {
  tab <- merge(x, a, sort = FALSE, by.x = by.x, by.y = by.anno)
  colnames(tab)[1] <- rowHeader
  return(tab)
}

printGgplot <- function(p, file, dir, height = 8, width = 8) {
  # avoid too large pngs, which will crash (arbitrary cutoff)
  if(height/width < 100) {
    ggsave(file=file.path(dir, paste0(file, ".png")), plot=p, width=width, height=height, units="cm") # height = 8, width = 8
    #png(file.path(dir, paste0(file, ".png")),height = round(600*height/width), width = 600)
    #print(p)
    #dev.off()
  }
  ggsave(file=file.path(dir, paste0(file, ".pdf")), plot=p, width=width, height=height, units="cm")
  #pdf(file.path(dir, paste0(file, ".pdf")), height = height, width = width)
  #print(p)
  #dev.off()
}





#----- Functions to plot data on lab-plate layout

readPlateLayoutFromExcel <- function(plateFile) {
	require(openxlsx)
	require(reshape2)

	#------- Get the plate layouts

	plates <- read.xlsx(plateFile, colNames = FALSE)

	if(nrow(plates) %% 17 != 0 ) stop("[ERROR] Number of rows must be divisible by 17 (16 rows plus 1 header per plate)")
	if(ncol(plates)       != 25) stop("[ERROR] Number of columns must be exactly 25 (24 columns plus 1 header per plate)")

	#------- Split into plates (one matrix per plate)

	groups <- ceiling(seq_len(nrow(plates))/17) # 17, because a plate has 16 rows, plus 1 header
	plates <- split(plates, groups)
	names(plates) <- sapply(plates, function(x) x[1,1]) # name is in the top left corner of each plate

	plates_new <- lapply(plates, function(x) {
		m <- as.matrix(x[-1,-1])                   # just well contents
		dimnames(m) <- list(x[-1,1], x[1, -1])     # row and column names
		allNA <- function(x) all(is.na(x))         # identify rows with all NA (empty values)
		m[!apply(m, 1, allNA),!apply(m, 2, allNA)] # remove those rows
	})


	#------- Reshape each matrix and combine in one table

	plates_melt <- lapply(plates_new, melt, varnames = c("row", "column")) # convert from xy matrix to columns x and y
	plates_melt <- lapply(names(plates_melt), function(x) data.frame(plates_melt[[x]], plate = x)) # add plate column
	all_plates <- do.call(rbind, plates_melt) # combine all plates in one table

	all_plates$row    <- factor(all_plates$row, levels = sort(levels(all_plates$row), decreasing = TRUE))
	all_plates$column <- factor(all_plates$column, levels = unique(all_plates$column))
	all_plates$plate  <- factor(all_plates$plate, levels = names(plates))

	all_plates$value[all_plates$value == ""] <- NA

	if(any(duplicated(na.omit(all_plates$value)))) stop("[ERROR] One or more barcodes were found in more than one well in the plate layout excel sheet. Each barcode should be unique!")

	return(all_plates)
}



plotDataOnPlates <- function(x, fileName, plates, invertSet, log = FALSE, skip = 20, symmetric = FALSE, ValueName = "Reads", addSummaries = TRUE) {
	x <- as.data.frame(x)

	barcodes <- plates[ ,"value"]

	# check that all barcodes are in the plates
	problemRows <- rowSums(x)>0 & ! rownames(x) %in% barcodes
	absent <- rownames(x[problemRows, , drop = FALSE])
	counts <-  rowSums(x[problemRows, , drop = FALSE])
	if(any(problemRows )) message("WARNING: Some barcodes found in data were not present in the supplied barcode plate layout: ",
	                              paste(
	                                paste0(
	                                  apply(
	                                    data.frame(absent, counts)[rev(order(counts)),],
	                                    1,
	                                    paste0,
	                                    collapse = " ("
	                                   ),
	                                  " reads)"
	                                  ),
	                                collapse = ", "
	                                )
	                              )

	if(addSummaries) {
	  dat <- cbind(x, sum_of_all = rowSums(x), present_in_any = rowSums(x != 0))
	} else {
	  dat <- cbind(x)
	}
	
	if(log & any(dat < 0)) stop("[ERROR] All values must be non-negative when log = TRUE.") 
	
	if(log) symmetric <- FALSE   # if log scale, scales shouldn't be symmetri around zero 

	twoColor <- !missing(invertSet)

	if(twoColor) symmetric <- TRUE # if two colors are used (with invertSet) the scale should be symmetric around zero, to ensure that values = 0 will appear grey

	if(twoColor & any(dat < 0)) stop("[ERROR] All values must be non-negative when invertSet is used.") 

	if(twoColor) names(invertSet) <- rownames(x)

	plot_plates <- lapply(colnames(dat), function(i) {

		# Combine plate layout with barcode read counts
		my_plates <- plates[,c("row", "column", "plate")]
		my_plates$value <- dat[barcodes, as.character(i)]

		if(twoColor)	{
			my_plates$invert <- invertSet[barcodes]
			my_plates$invert[is.na(my_plates$invert)] <- FALSE
		}

		# split by plate
		tmp <- split.data.frame(my_plates, my_plates$plate)

		# skip plates where all wells are NA (NA means not present in the dat matrix at all)
		dropPlates <- sapply(tmp, function(j) all(is.na(j$value)))
		
		# combine the remaining plates in one matrix
		do.call(rbind, tmp[!dropPlates])
	})
	
	names(plot_plates) <- colnames(dat)

	#-- Plot the read counts in the plate layout
	pdf(paste0(fileName, ".pdf"))
	
	print_plates <- list()

	for (i in 1:length(plot_plates)) {

		x <- plot_plates[[i]]

		#----- Generate matrices for excel
		if(twoColor) x$invert <- NULL
		
		each_plate <- split(x, x$plate)
		each_plate <- each_plate[sapply(each_plate, nrow) > 0] 

		each_matrix <- lapply(names(each_plate), function(i) {
			x <- each_plate[[i]]
			my_plate <- dcast(row ~ column, data = x, value.var = "value", fill = 0, drop = FALSE)
			my_plate <- my_plate[order(as.character(my_plate[,1])), ]
			colnames(my_plate)[1] <- i
			rbind(colnames(my_plate), as.matrix(my_plate), rep("", ncol(my_plate)))
		})
		each_matrix[[1]]

		all_matrices <- as.data.frame(do.call(rbind, each_matrix))
		all_matrices[,-1] <- apply(all_matrices[,-1], 2, as.numeric)

		#----- Store matrices in a list
		print_plates[[i]] <- all_matrices

		# ---- Skip plotting if an akey has values in less than 20 wells
		if( sum(x[,"value"] != 0, na.rm = TRUE) < skip) { next }

		x <- plot_plates[[i]]
		
		if(log) x$value <- log10(x$value + 1)
		if(twoColor) x$value[x$invert] <- - x$value[x$invert]

		message("\nAKEY", names(plot_plates)[i])
		
		#---- Generate plot
		p <- ggplot(data = x, aes(x = column, y = row, fill = value)) +
			facet_wrap(~plate, ncol = 2) + geom_tile() +
			ggtitle(names(plot_plates)[i]) +
			theme(axis.text=element_text(size=5),
			text = element_text(size=6), legend.text=element_text(size=5)
			)

		if(twoColor) {
			myBreaks <- myLabels <- pretty(abs(x$value))
			maxValue <- max(abs(x$value), na.rm = TRUE)
		} else {
			myBreaks <- myLabels <- pretty(x$value)
			maxValue <- max(abs(x$value), na.rm = TRUE)
		}

		if(log) {

			maxValue <- max(abs(x$value), na.rm = TRUE)
			maxValue <- max(1, maxValue)

			myBreaks <- 1:(1+ceiling(maxValue))
			myLabels  <- 10^(abs(myBreaks))
			myBreaks <- log10(myLabels - 1)

			myBreaks <- c(0, myBreaks)
			myLabels <- c(0, myLabels)
		}

		if(twoColor) {
			myBreaks <- myBreaks[myBreaks!=0]
			myBreaks  <- c(rev(-myBreaks), 0, myBreaks)
			
			myLabels <- myLabels[myLabels!=0]
			myLabels <-  c(rev(myLabels), 0, myLabels)

			myLabels <- as.character(myLabels)
			myLabels[1] <- "Not in panel"
			myLabels[length(myLabels)] <- "In panel"

			legend_name <- paste0("in\npanel\n\n\n", ValueName, "\n\n\nnot in\npanel")
			# cols <- c(scales::muted("red"), "grey85", scales::muted("blue"))
			# key_limits <- c(-ceiling(maxValue), ceiling(maxValue))

		} else {
			
			legend_name <- ValueName

		}

		if(symmetric) {

			cols <- c(scales::muted("red"), "grey85", scales::muted("blue"))
			key_limits <- c(-ceiling(maxValue), ceiling(maxValue))
		
		} else {

			cols <- c("grey85", scales::muted("blue"))
			key_limits <- c(0, ceiling(maxValue))
		}

			message("breaks: ", paste(myBreaks, collapse = ", "))
			message("labels: ", paste(myLabels, collapse = ", "))
			message("maxValue: ", maxValue)
			message("limits: ", paste(key_limits, collapse = ", "))

		p <- p + scale_fill_gradientn(
				name = legend_name,
				colours = cols,
				na.value = "grey85",
                breaks = myBreaks,
                labels = myLabels,
                limits = key_limits
                )
		p <- p + guides(fill = guide_colorbar( barheight = 20, title.position = "left", title.hjust = 0.5, draw.ulim = FALSE, draw.llim = FALSE))	
		p <- p + theme(panel.grid.minor = element_blank(), 
			       panel.grid.major = element_blank(), 
			       panel.background = element_blank(),
		               text = element_text(size=6), legend.text=element_text(size=5), legend.title = element_text(size=5), legend.key.size = unit(0.3, 'cm')
		)
	

		#----- Make plot		
		print(p)

	}

	dev.off()

	#----- Print all matrices to excel sheets
	names(print_plates) <- names(plot_plates)
	write.xlsx(print_plates, file = paste0(fileName, ".xlsx"), colNames = FALSE)

	invisible(plot_plates)

}

