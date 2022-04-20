#plot allele frequencies
library(tidyverse)
library(ggplot2)
library(here)
library(gplots)
#Get populations to loop analysis over
POP_SPEC <- stringr::str_split(read_file(file.path(here::here(), "vcf_analyses/admixture/pop-species.4.list")), pattern = ",")[[1]]
#Loop bringing in data by population, cleaning up for fixed frequencies, and plotting allele frequency distribution
#Import nucleotide diversity data
pi <- read_delim(file.path(here::here(), "vcf_analyses/popstats/nucleotide_diversity/FILE.windowed.pi"), delim = "\t", col_names = c("CHROM", "BIN_START", "BIN_END", "N_VARIANT", "Pi", "loc"), skip = 1)
#correct for 0 values excluded by vcftools
#get mean site pi values, removing "-nan"s first
pi_window_scores <- aggregate(pi[, 5], list(pi$loc), mean)
colnames(pi_window_scores)[1] <- "loc"
write.csv(pi_window_scores, file=file.path(here::here(), "vcf_analyses/popstats/nucleotide_diversity/FILE_pi_site_mean.csv"))
#Plot data; plots won't be accurate until 0 values are appended into data.frame, the number of which will depend on genome size and number of windows without variants
pdf(file=file.path(here::here(), "results/nuclear_stats/FILE_window_pi.pdf"), width = 6, height = 3)
p <- ggplot(pi, aes(x=loc, y=Pi, fill=loc)) + geom_violin() + geom_point(data=pi_window_scores, shape=24, size=2) + ylim(0,0.01)
p
dev.off()