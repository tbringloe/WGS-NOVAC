#plot allele frequencies
library(tidyverse)
library(ggplot2)
library(here)
library(gplots)
#Get populations to loop analysis over
POP_SPEC <- stringr::str_split(read_file(file.path(here::here(), "vcf_analyses/admixture/pop-species.4.list")), pattern = ",")[[1]]
#Loop bringing in data by population, cleaning up for fixed frequencies, and plotting allele frequency distribution
for(l in POP_SPEC) {
  nam <- paste(l,"_AFS.pdf", sep = "")
  tmp1=paste(l,".frq", sep = "")
  tmp2 <- read_delim(file.path(here::here(), "vcf_analyses/popstats/allele_frequencies",tmp1), delim = "\t", col_names = c("chr", "pos", "nalleles", "nchr", "a1", "a2"), skip = 1)
  tmp3 <- subset(tmp2, a1 > 0)
  tmp4 <- subset(tmp3, a1 < 1)
  pdf(file=file.path(here::here(), "results/nuclear_stats/",nam), width=8, height=2)
  a <- print(ggplot(tmp4, aes(a1)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3))
  a + theme_light() + ylim(0,8)
  dev.off()
}
#locally working code
#for(l in POP_SPEC) {
#  nam <- paste(l,"_AFS.pdf", sep = "")
#  tmp1=paste(l,".frq", sep = "")
#  tmp2 <- read_delim(file.path("vcf_analyses/popstats/allele_frequencies",tmp1), delim = "\t", col_names = c("chr", "pos", "nalleles", "nchr", "a1", "a2"), skip = 1)
#  tmp3 <- subset(tmp2, a1 > 0)
#  tmp4 <- subset(tmp3, a1 < 1)
#  pdf(file=file.path("vcf_analyses/popstats/allele_frequencies/",nam), width=8, height=2)
#  a <- print(ggplot(tmp4, aes(a1)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3))
#  a + theme_light() + ylim(0,8)
#  dev.off()
#}
#Loop bringing in data by population and plotting heterozygosity values
#clean up environment
rm(list = ls())
#bring in population list again
POP_SPEC <- stringr::str_split(read_file(file.path(here::here(), "vcf_analyses/admixture/pop-species.4.list")), pattern = ",")[[1]]
#loop over importing data for each population and adding column for population
het <- read_delim(file.path(here::here(), "vcf_analyses/popstats/heterozygosity/FILE.het"), delim = "\t", col_names = c("INDV", "O(HOM)", "E(HOM)", "N_SITES", "F"), skip = 1)
  loc <- rep(NA, length(het$INDV))
  for(l in POP_SPEC) {
  loc[grep(l, het$INDV)] <- l
  }
  het <- as_tibble(data.frame(het, loc))
#Plot data
pdf(file=file.path(here::here(), "results/nuclear_stats/FILE_Het.pdf"), width = 6, height = 3)
p <- ggplot(het, aes(x=loc, y=F, fill=loc)) + geom_boxplot()
p
dev.off()
#locally working code
#  het <- read_delim(file.path("vcf_analyses/popstats/heterozygosity/FILE.het"), delim = "\t", col_names = c("INDV", "O(HOM)", "E(HOM)", "N_SITES", "F"), skip = 1)
#  loc <- rep(NA, length(het$INDV))
#  for(l in POP_SPEC) {
#  loc[grep(l, het$INDV)] <- l
#  }
#  het <- as_tibble(data.frame(het, loc))
#Plot data
#pdf("FILE_Het.pdf", width = 6, height = 3)
#p <- ggplot(het, aes(x=loc, y=F, fill=loc)) + geom_boxplot()
#p
#dev.off()
#clean up environment
rm(list = ls())
#bring in population list again
POP_SPEC <- stringr::str_split(read_file(file.path(here::here(), "vcf_analyses/admixture/pop-species.4.list")), pattern = ",")[[1]]
#Import TajimaD data
TajimaD <- read_delim(file.path(here::here(), "vcf_analyses/popstats/TajimaD/FILE.Tajima.D"), delim = "\t", col_names = c("CHROM", "BIN_START", "N_SNPS", "TajimaD", "loc"), skip = 1)
#remove NaNs then get average values
TajimaD2 <- subset(TajimaD, TajimaD!="NaN")
TajimaD_scores <- aggregate(TajimaD2[, 4], list(TajimaD2$loc), mean)
colnames(TajimaD_scores)[1] <- "loc"
write.csv(TajimaD_scores, file=file.path(here::here(), "vcf_analyses/popstats/TajimaD/FILE_TajimaD_mean.csv"))
#Plot data
pdf(file=file.path(here::here(), "results/nuclear_stats/FILE_TajimaD.pdf"), width = 6, height = 3)
p <- ggplot(TajimaD, aes(x=loc, y=TajimaD, fill=loc)) + geom_violin() + geom_point(data=TajimaD_scores, shape=24, size=2)
p
dev.off()
#clean up environment
rm(list = ls())
#bring in population list again
POP_SPEC <- stringr::str_split(read_file(file.path(here::here(), "vcf_analyses/admixture/pop-species.4.list")), pattern = ",")[[1]]
#Import nucleotide diversity data
pi <- read_delim(file.path(here::here(), "vcf_analyses/popstats/nucleotide_diversity/FILE.sites.pi"), delim = "\t", col_names = c("CHROM", "POS", "Pi", "loc"), skip = 1)
#get mean site pi values, removing "-nan"s first
pi2 <- subset(pi, Pi!="-nan")
pi_site_scores <- aggregate(pi2[, 3], list(pi$loc), mean)
colnames(pi_site_scores)[1] <- "loc"
write.csv(pi_site_scores, file=file.path(here::here(), "vcf_analyses/popstats/nucleotide_diversity/FILE_pi_site_mean.csv"))
#Plot data
pdf(file=file.path(here::here(), "results/nuclear_stats/FILE_sites_pi.pdf"), width = 6, height = 3)
p <- ggplot(pi, aes(x=loc, y=Pi, fill=loc)) + geom_violin() + geom_point(data=pi_site_scores, shape=24, size=2)
p
dev.off()
#clean up environment
rm(list = ls())
#Import windowed nucleotide diversity data
pi <- read_delim(file.path(here::here(), "vcf_analyses/popstats/nucleotide_diversity/FILE.windowed.pi"), delim = "\t", col_names = c("CHROM", "BIN_START", "BIN_END", "N_VARIANT", "Pi", "loc"), skip = 1)
#correct for 0 values excluded by vcftools
#get mean site pi values, removing "-nan"s first
pi_window_scores <- aggregate(pi[, 5], list(pi$loc), mean)
colnames(pi_window_scores)[1] <- "loc"
write.csv(pi_window_scores, file=file.path(here::here(), "vcf_analyses/popstats/nucleotide_diversity/FILE_pi_site_mean.csv"))
#Plot data; set ylim value depending on scaling of values
pdf(file=file.path(here::here(), "results/nuclear_stats/FILE_window_pi.pdf"), width = 6, height = 3)
p <- ggplot(pi, aes(x=loc, y=Pi, fill=loc)) + geom_violin() + geom_point(data=pi_window_scores, shape=24, size=2) + ylim(0,0.025)
p
dev.off()
#clean up environment
rm(list = ls())
#Import relatedness AJK data
relat <- read_delim(file.path(here::here(), "vcf_analyses/popstats/relatedness/FILE.relatedness"), delim = "\t", col_names = c("INDV1", "INDV2", "RELATEDNESS_AJK"), skip = 1)
individuals <- read.table("sample.list2")
ind <- paste0(individuals[,1])
#colnames(relat)[1] <- "INDV1"
#colnames(relat)[2] <- "INDV2"
#colnames(relat)[3] <- "RELATEDNESS_AJK"
relat = relat[-1,]
relat.wide <- pivot_wider(data=relat, names_from = INDV2, values_from = RELATEDNESS_AJK)
relat.matrix <- apply(as.matrix(relat.wide), 2, as.numeric)
pdf(file=file.path(here::here(), "results/nuclear_stats/FILE_relatedness_AJK.pdf"), width = 10, height = 10)
p <- heatmap.2(relat.matrix, trace="none", Rowv=NA, Colv=NA, cexRow=0.6,cexCol = 0.6, labRow = ind, labCol = ind, col= colorRampPalette(c("lemonchiffon", "lemonchiffon", "lemonchiffon", "yellow", "red", "magenta", "midnightblue", "black"))(70))
p
dev.off()
#Import relatedness PHI data
relat2 <- read_delim(file.path(here::here(), "vcf_analyses/popstats/relatedness/FILE.relatedness2"), delim = "\t", col_names = c("INDV1", "INDV2", "N_AaAa", "N_AAaa", "N1_Aa", "N2_Aa", "RELATEDNESS_PHI"), skip = 1)
relat2 <- subset(relat2, select = -c(N_AaAa, N_AAaa, N1_Aa, N2_Aa))
#colnames(relat2)[1] <- "INDV1"
#colnames(relat2)[2] <- "INDV2"
#colnames(relat2)[3] <- "RELATEDNESS_PHI"
relat2 = relat2[-1,]
relat2.wide <- pivot_wider(data=relat2, names_from = INDV2, values_from = RELATEDNESS_PHI)
relat2.matrix <- apply(as.matrix(relat2.wide), 2, as.numeric)
pdf(file=file.path(here::here(), "results/nuclear_stats/FILE_relatedness_PHI.pdf"), width = 10, height = 10)
p <- heatmap.2(relat2.matrix, trace="none", Rowv=NA, Colv=NA, cexRow=0.6,cexCol = 0.6, labRow = ind, labCol = ind, col= colorRampPalette(c("lemonchiffon", "lemonchiffon", "lemonchiffon", "yellow", "red", "magenta", "midnightblue", "black"))(70))
p
dev.off()