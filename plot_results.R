library(tidyverse)
library(ggseqlogo)
library(openxlsx)
library(writexl)
library(readxl)
library(stringr)

########################################
# Functions
########################################

rbindAllSheets <- function(file) {
  # This function rbinds a single dataframe with the content of multiple sheets from the same Excel file
  # (assuming that all the sheets have the same columns)
  wb <- loadWorkbook(file)
  sheets <- sheets(wb)
  do.call(rbind, lapply(sheets, function(sheet) {
    readWorkbook(file, sheet) })
  ) 
}

#########################################################################
## Functions 
#########################################################################

rbindAllSheets <- function(file) {
  # This function rbinds a single dataframe with the content of multiple sheets from the same Excel file
  # (assuming that all the sheets have the same columns)
  wb <- loadWorkbook(file)
  sheets <- sheets(wb)
  
  df <- do.call(rbind, lapply(sheets, function(sheet) {
    df_sheet <- readWorkbook(file, sheet) 
    df_sheet 
  })
  ) 
  return(df)
}

# Working path ------------------------------------------------------------

data_dir <- 'Barracoda-2.0/results/Nanopore_results/barracoda_2026-03-25.1/experiment_exp1'
setwd(data_dir)

# Load data ------------------------------------------------------------

# set fold change threshold (default=2)
FC_threshold=2

# Use the rbindAllSheets() function to bind all sheets in the fold_change.xlsx file
#logFC_df <- rbindAllSheets(paste0(data_dir,"/fold_change.xlsx")) %>% as_tibble()

# Or used the fold_change_by_MHC.xlsx file for MHC corrected values
logFC_df <- rbindAllSheets(paste0(data_dir,"/fold_change_by_MHC.xlsx")) %>% as_tibble() %>% drop_na(log_fold_change)

# plot data ------------------------------------------------------------

# make dataframe for plotting
plot_df <- logFC_df %>% 
            mutate(
              colors = ifelse(log_fold_change > FC_threshold,"enriched", "not enriched"),
              colors = factor(colors,levels = c("enriched", "not enriched"))
            )

color_list <- c('enriched'='#b2182b', 'not enriched'='gray')

logFC_plot <- ggplot(data=plot_df, aes(x=barcode, y=log_fold_change, color=colors)) + 
  geom_point(alpha=0.8) + 
  #scale_y_continuous(limits = c(0, NA)) +
  #scale_x_discrete(expand = expansion(add = 0.1)) + 
  geom_hline(yintercept = FC_threshold, linetype = "dashed", color='#636363') +
  scale_color_manual(values=color_list) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size=14),
        panel.spacing.x = unit(0,"line"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 3),
        #legend.position="bottom", legend.direction = "horizontal",legend.box = "vertical",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.background.x = element_rect(color = NA,  fill=NA)
  ) + labs(title='', x="pMHC barcodes", y="Log2FC", color='') 

logFC_plot
#ggsave(plot=logFC_plot,filename='Barracoda-2.0/plots/log2FC_color_enrichment.png', height=6.0, width=13.0) 

# split per HLA
plot_per_HLA <- logFC_plot + facet_grid(~HLA, scales = 'free', space = 'free') + 
  theme(strip.text.x = element_text(size=9, angle=90))
plot_per_HLA

ggsave(plot=plot_per_HLA,filename='Barracoda-2.0/plots/log2FC_per_HLA_color_enrichment.png', height=4.5, width=14.0) 

# split per sample 
plot_per_sample <- logFC_plot + facet_wrap(~sample, scales = 'free_x')
plot_per_sample

#ggsave(plot=plot_per_sample,filename='Barracoda-2.0/plots/log2FC_per_sample_color_enrichment.png', height=10.0, width=13.0) 

# color each sample ---------------------------------------------------

# sample names
samples <- logFC_df %>% dplyr::pull(sample) %>% unique()

library(RColorBrewer)
base_colors <- brewer.pal(12, "Paired")
sample_colors <- colorRampPalette(base_colors)(length(samples))
color_list <- c(setNames(sample_colors, samples),'not enriched'='gray')

plot_df <- logFC_df %>% 
  mutate(colors = ifelse(log_fold_change > FC_threshold, sample, "not enriched"))

logFC_plot <- ggplot(data=plot_df, aes(x=barcode, y=log_fold_change, color=colors)) + 
  geom_point(alpha=0.8) + 
  #scale_y_continuous(limits = c(0, NA)) +
  #scale_x_discrete(expand = expansion(add = 0.1)) + 
  geom_hline(yintercept = FC_threshold, linetype = "dashed", color='#636363') +
  scale_color_manual(values=color_list) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size=14),
        panel.spacing.x = unit(0,"line"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 3),
        #legend.position="bottom", legend.direction = "horizontal",legend.box = "vertical",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.background.x = element_rect(color = NA,  fill=NA)
  ) + labs(title='', x="pMHC barcodes", y="Log2FC", color='') 

logFC_plot
#ggsave(plot=logFC_plot,filename='Barracoda-2.0/plots/log2FC_color_samples.png', height=6.0, width=13.0) 

# split per HLA 
plot_per_HLA <- logFC_plot + facet_grid(~HLA, scales = 'free', space = 'free') + 
  theme(strip.text.x = element_text(size=9, angle=90), legend.position="bottom") + labs(color='sample')
plot_per_HLA

ggsave(plot=plot_per_HLA,filename='Barracoda-2.0/plots/log2FC_per_HLA_color_samples.png', height=5.5, width=13.0) 
