# Dependencies exist from 01_vcfQC r scripts
library(ggplot2)
library(tidyverse)
library(here)
#load eigenvec values
pca <- read_table2(file.path(here::here(), "vcf_analyses/plink/", "FILE.eigenvec"), col_names = FALSE)
POP_SPEC <- stringr::str_split(read_file(file.path(here::here(), "vcf_analyses/admixture/pop-species.4.list")), pattern = ",")[[1]]
eigenval <- scan(file.path(here::here(), "vcf_analyses/plink/", "FILE.eigenval"))
#tidy up pca data
pca <- pca[,-1]
colnames(pca)[1] <- "ind"
colnames(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-1))
#add location information for the PCA plot
loc <- rep(NA, length(pca$ind))
for(l in POP_SPEC){
loc[grep(l, pca$ind)] <- l
}
pca <- as_tibble(data.frame(pca, loc))
#get percentages of variation explained by vectors
pve <- data.frame(PC = 1:20, pve = eigenval/sum(eigenval)*100)
#plot to visualise variance explained by vector
a <- ggplot(pve, aes(PC, pve)) + geom_bar(stat="identity")
#create vector image
pdf(file=file.path(here::here(), "results/pca/", "FILE_PCA_vectors.pdf"))
#insert ggplot code
a + ylab("Percentage variance explained") + theme_light()
dev.off()
#plot pca with % values on axes
b <- ggplot(pca, aes(PC1, PC2, col = loc)) + geom_point(size = 1)
b <- b + coord_equal() + theme_light()
b <- b + xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)"))
#create vector image
pdf(file=file.path(here::here(), "results/pca/", "FILE_PCA.pdf"))
#insert ggplot code
b
dev.off()