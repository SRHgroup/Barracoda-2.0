#
# This script creates the plots "barplot--read-lengths.[png|pdf]"
#

################################
###        PACKAGES          ###
################################

library(ggplot2)

################################
###  COMMAND LINE ARGUMENTS  ###
################################

file_in     <- commandArgs(TRUE)[1:3]  # files containing read lengths (all reads, after "filter 2 of 3", after "filter A and B")
dir_out     <- commandArgs(TRUE)[4]    # output directory
exp_length  <- commandArgs(TRUE)[5]    # expected length of barcode DNA

cat(file_in)
################################
###   READ AND PREPARE DATA  ###
################################

names(file_in) <- c("All", "Filter 2 of 3", "A and B")

x <- do.call(rbind, lapply(names(file_in), function(f) {
	x <- read.delim(file_in[f], as.is = TRUE)
	colnames(x) <- c("length", "reads")
	x$File <- f
	return(x)
}))

x$File <- factor(x$File, levels = names(file_in))
x$col <- factor(ifelse(x$length == exp_length, "expected length", ""))


################################
###     PLOT READ LENGTHS    ###
################################

ColorblindFriendlyColors <- function(n = 7) {
  return(c(   orange = rgb( 230, 159,   0,  max=255),
            sky_blue = rgb( 86,  180, 233,  max=255),
        bluish_green = rgb( 0,   158, 155,  max=255),
              yellow = rgb( 240, 228,  66,  max=255),
                blue = rgb( 0,   114, 178,  max=255),
          vermillion = rgb( 213, 94,    0,  max=255),
      reddish_purple = rgb( 204, 121, 167,  max=255))[1:n])
}

cls <- setNames(ColorblindFriendlyColors()[c(3,7)], levels(x$col))

p <- ggplot(x, aes(x = length, y = reads, fill = col))
p <- p + geom_bar(stat = "identity", position = "dodge")
p <- p + scale_fill_manual("Expected length", breaks=c("expected length"), labels = "", values = cls) # breaks=c("expected length") 
p <- p + theme(legend.justification=c(0,0), legend.position=c(0,0), legend.direction = "horizontal", 
	       text = element_text(size=6),
               legend.key.size = unit(0.3, 'cm'), legend.title = element_text(size=5), legend.text=element_text(size=5)
)
p <- p + ggtitle("Lengths of sequencing reads") + xlab("Read length") + ylab("Number of reads")
p <- p + facet_grid(File ~ .)

#ggsave(file=file.path(dir_out, "read-lengths.png"), plot=p, width = 600, height=round(600*7/5))
ggsave(file=file.path(dir_out, "read-lengths.png"), plot=p, width=8, height=8, units="cm")

#png(file.path(dir_out, "read-lengths.png"), height = 600, width = round(600*7/5))
#print(p)
#dev.off()

#pdf(file.path(dir_out, "read-lengths.pdf"), height = 5)
#print(p)
#dev.off()



