#!/usr/bin/env Rscript
# Generates a tiny, deterministic smoke-test fixture for Barracoda:
#   2 a-keys (one input + one sample) x 2 antigens (same HLA), with synthesised
#   reads built in the exact Barracoda read layout so they map cleanly.
#
# Read layout (verified against real IonTorrent reads):
#   sampleTag(10) + primerA(-B,25) + N6A(6) + epitopeA(25) + anneal(-E,16)
#   + epitopeB(25) + N6B(6) + revcomp(primerB, -H)
#
# Run from the repo root:  Rscript test_data/smoke/make_smoke_data.R
suppressMessages(library(openxlsx))
set.seed(42)  # deterministic N6 barcodes -> reproducible fixture

src   <- "test_data"                       # source example data (gitignored, from server)
outdir<- "tests/smoke"                      # committed fixture lives here
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

## --- constant sequences (the canonical IonTorrent primers, from test_data/howtorun) ---
primerA <- "GAAGTTCCAGCCAGCGTCACAGTTT"   # -B
anneal  <- "GGTCAGCATCATTTCC"            # -E
primerH <- "GTTATCGGCTCGTTCACACTCGA"     # -H (appears reverse-complemented in the read)
revcomp <- function(s) {
  chartr("ACGT", "TGCA", paste(rev(strsplit(s, "")[[1]]), collapse = ""))
}
primerB_in_read <- revcomp(primerH)

## --- read a fasta into a named character vector ---
read_fasta <- function(f) {
  ln <- readLines(f); ln <- ln[nzchar(trimws(ln))]
  hdr <- grepl("^>", ln)
  setNames(ln[!hdr], sub("^>\\s*", "", ln[hdr]))
}
tags <- read_fasta(file.path(src, "sample_id_tags.fasta"))
epiA <- read_fasta(file.path(src, "a.fasta"))
epiB <- read_fasta(file.path(src, "b.fasta"))

## --- choices: experiment RCC_17, 1 input + 1 sample a-key, 2 antigens (HLA A0201) ---
akeys <- c(input = "A-Key_2OS_F1_70", sample = "A-Key_2OS_F1_66")
akey_sample_name <- c("A-Key_2OS_F1_70" = "input", "A-Key_2OS_F1_66" = "RCC 17 PE+")
# antigens as (A-oligo, B-oligo, peptide, HLA); both A0201 so the MHC path runs.
# Use DISTINCT A and B oligos so both epitope fastas hold >=2 sequences:
# dissect-barcodes--fast.pl mis-handles an epitope fasta with a single sequence
# (its 1-entry fasta reader flattens the structure, breaking the ${pep}{id} lookup).
antigens <- data.frame(
  A   = c("A1", "A2"),
  B   = c("B61", "B62"),
  Barcode = c("A1B61", "A2B62"),
  Peptide = c(192, 193),
  HLA = c("A0201", "A0201"),
  RCC.sample = c("RCC 17", "RCC 17"),
  stringsAsFactors = FALSE
)

oligoA_name <- function(a) grep(paste0("Oligo_", a, "$"), names(epiA), value = TRUE)
oligoB_name <- function(b) grep(paste0("Oligo_", b, "$"), names(epiB), value = TRUE)

## --- sanity: required sequences exist ---
stopifnot(all(akeys %in% names(tags)))
for (a in unique(antigens$A)) stopifnot(length(oligoA_name(a)) == 1)
for (b in unique(antigens$B)) stopifnot(length(oligoB_name(b)) == 1)

## --- synthesise reads ---
rand_n6 <- function() paste(sample(c("A","C","G","T"), 6, replace = TRUE), collapse = "")
N_PER_CELL <- 40  # reads per (a-key x antigen)
fq <- character(0); rid <- 0
for (ak in akeys) {
  for (i in seq_len(nrow(antigens))) {
    oA <- epiA[[oligoA_name(antigens$A[i])]]
    oB <- epiB[[oligoB_name(antigens$B[i])]]
    for (k in seq_len(N_PER_CELL)) {
      rid <- rid + 1
      # epitope A is forward in the read; epitope B and primer B are reverse-complemented
      seq <- paste0(tags[[ak]], primerA, rand_n6(), oA, anneal, revcomp(oB), rand_n6(), primerB_in_read)
      fq <- c(fq, paste0("@SMOKE:", sprintf("%05d", rid)), seq, "+",
              paste(rep("I", nchar(seq)), collapse = ""))
    }
  }
}
writeLines(fq, file.path(outdir, "reads.fastq"))

## --- tiny tag / epitope fastas (only what the 2 a-keys / 2 antigens need) ---
writeLines(as.vector(rbind(paste0(">", names(tags[akeys])), tags[akeys])),
           file.path(outdir, "tags.fasta"))
aset <- unique(sapply(unique(antigens$A), oligoA_name))
bset <- unique(sapply(unique(antigens$B), oligoB_name))
writeLines(as.vector(rbind(paste0(">", aset), epiA[aset])), file.path(outdir, "epitopeA.fasta"))
writeLines(as.vector(rbind(paste0(">", bset), epiB[bset])), file.path(outdir, "epitopeB.fasta"))

## --- sampleID (no header, 3 cols) and annotation (sheet per experiment) ---
sid <- data.frame(akey = akeys, sample = akey_sample_name[akeys], experiment = "RCC_17")
write.xlsx(sid, file.path(outdir, "sampleID.xlsx"), colNames = FALSE)
write.xlsx(list(RCC_17 = antigens[, c("Barcode","Peptide","HLA","RCC.sample")]),
           file.path(outdir, "annotations.xlsx"))

cat("Wrote smoke fixture to", outdir, "\n")
cat("  reads.fastq :", rid, "reads (",N_PER_CELL,"per a-key x antigen )\n")
cat("  a-keys      :", paste(akeys, collapse=", "), "\n")
cat("  antigens    :", paste(antigens$Barcode, collapse=", "), "(HLA A0201)\n")
