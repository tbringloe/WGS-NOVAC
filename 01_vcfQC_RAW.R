# install dependencies
#install.packages("here", "tidyverse", "rprojroot", "ggplot2")
# load libraries
library(tidyverse)
library(here)
library(rprojroot)
library(ggplot2)

# Plot histogram of quality of site calls using a phrep scoring system
var_qual <- read_delim(file.path(here::here(), "vcf_analyses/vcf_stats/raw/FILE_raw.lqual"), delim = "\t",
                       col_names = c("chr", "pos", "qual"), skip = 1)
qual.gg <- ggplot(var_qual, aes(qual)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
qual.gg <- qual.gg + theme_light() + xlim(0, 999)
ggsave(filename = file.path(here::here(), "results/vcf_QC/raw/RAW_qual_FILE.pdf"),
       plot = qual.gg, 
       width = 8,
       height = 4,
       units = c("in"), bg = "white")

# Plot histogram of mean site depth, with an upper plot depth limit of 100
var_depth <- read_delim(file.path(here::here(), "vcf_analyses/vcf_stats/raw/FILE_raw.ldepth.mean"), delim = "\t",
                        col_names = c("chr", "pos", "mean_depth", "var_depth"), skip = 1)
site.depth.gg <- ggplot(var_depth, aes(mean_depth)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
site.depth.gg <- site.depth.gg + theme_light() + xlim(0, 100)
ggsave(filename = file.path(here::here(), "results/vcf_QC/raw/RAW_siteDepth_FILE.pdf"),
       plot = site.depth.gg, 
       width = 8,
       height = 4,
       units = c("in"), bg = "white")


# Plot histogram of mean missingness per site
var_miss <- read_delim(file.path(here::here(), "vcf_analyses/vcf_stats/raw/FILE_raw.lmiss"), delim = "\t",
                       col_names = c("chr", "pos", "nchr", "nfiltered", "nmiss", "fmiss"), skip = 1)
site.miss.gg <- ggplot(var_miss, aes(fmiss)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
site.miss.gg <- site.miss.gg + theme_light()
ggsave(filename = file.path(here::here(), "results/vcf_QC/raw/RAW_siteMiss_FILE.pdf"),
       plot = site.miss.gg, 
       width = 8,
       height = 4,
       units = c("in"), bg = "white")

# Plot histogram of minor allele frequency
var_freq <- read_delim(file.path(here::here(), "vcf_analyses/vcf_stats/raw/FILE_raw.frq"), delim = "\t",
                       col_names = c("chr", "pos", "nalleles", "nchr", "a1", "a2"), skip = 1)
var_freq$maf <- var_freq %>% select(a1, a2) %>% apply(1, function(z) min(z))
maf.gg <- ggplot(var_freq, aes(maf)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
maf.gg <- maf.gg + theme_light()
ggsave(filename = file.path(here::here(), "results/vcf_QC/raw/RAW_maf_FILE.pdf"),
       plot = maf.gg, 
       width = 8,
       height = 4,
       units = c("in"), bg = "white")

# Plot histogram of mean depth per individual
ind_depth <- read_delim(file.path(here::here(), "vcf_analyses/vcf_stats/raw/FILE_raw.idepth"), delim = "\t",
                        col_names = c("ind", "nsites", "depth"), skip = 1)
ind.depth.gg <- ggplot(ind_depth, aes(depth)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3)
ind.depth.gg <- ind.depth.gg + theme_light()
ggsave(filename = file.path(here::here(), "results/vcf_QC/raw/RAW_indDepth_FILE.pdf"),
       plot = ind.depth.gg, 
       width = 8,
       height = 4,
       units = c("in"), bg = "white")

# Plot histogram of proportion missing data per individual
ind_miss  <- read_delim(file.path(here::here(), "vcf_analyses/vcf_stats/raw/FILE_raw.imiss"), delim = "\t",
                        col_names = c("ind", "ndata", "nfiltered", "nmiss", "fmiss"), skip = 1)
ind.miss.gg <- ggplot(ind_miss, aes(fmiss)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3)
ind.miss.gg <- ind.miss.gg + theme_light()
ggsave(filename = file.path(here::here(), "results/vcf_QC/raw/RAW_indMiss_FILE.pdf"),
       plot = ind.miss.gg, 
       width = 8,
       height = 4,
       units = c("in"), bg = "white")

# Plot histogram of inbreeding coefficients/heterozygosity
ind_het <- read_delim(file.path(here::here(), "vcf_analyses/vcf_stats/raw/FILE_raw.het"), delim = "\t",
                      col_names = c("ind","ho", "he", "nsites", "f"), skip = 1)
ind.het.gg <- ggplot(ind_het, aes(f)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3)
ind.het.gg <- ind.het.gg + theme_light()
ggsave(filename = file.path(here::here(), "results/vcf_QC/raw/RAW_indHet_FILE.pdf"),
       plot = ind.het.gg, 
       width = 8,
       height = 4,
       units = c("in"), bg = "white")