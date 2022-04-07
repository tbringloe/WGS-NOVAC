#install dependencies
#install.packages("PopGenome")
#load dependencies
library(PopGenome)
library(here)
#load alignments
Mito <- readData(file.path(here::here(), "results/organellar/Popgenome/mito"))
Plastid <- readData(file.path(here::here(), "results/organellar/Popgenome/chloro"))
#define populations
POP_SPEC <- stringr::str_split(read_file(file.path(here::here(), "vcf_analyses/admixture/pop-species.2.list")), pattern = ",")[[1]]
for(l in POP_SPEC) {
  nam <- paste(l, sep = "")
  tmp1=paste(l,".indv", sep = "")
  tmp2 <- stringr::str_split(read_file(tmp1), pattern = ",")[[1]]
  assign(nam, tmp2)
}
Mito <- set.populations(Mito,list(POPS))
Plastid <- set.populations(Plastid,list(POPS))
#check populations are correctly assigned
Mito@populations
Plastid@populations
#run tests and display values for each population
#Tajima's D test for neutrality
Mito <- neutrality.stats(Mito)
Mito_TajimasD <- Mito@Tajima.D
write.csv2(Mito_TajimasD, file=file.path(here::here(), "results/organellar/mito/FILE_Mito_mapped_TajimasD.csv"))
Plastid <- neutrality.stats(Plastid)
Plastid_TajimasD <- Plastid@Tajima.D
write.csv2(Plastid_TajimasD, file=file.path(here::here(), "results/organellar/plastid/FILE_Plastid_mapped_TajimasD.csv"))
#Nucleotide diversity, pi; note, nucleotide diversity is calculated within and between populations, but needs to be divided by number of sites to yield pi (average number of nucleotide difference between ind. per site)
Mito <- F_ST.stats(Mito)
Mito_pi <- Mito@nuc.diversity.within/Mito@n.sites
write.csv2(Mito_pi, file=file.path(here::here(), "results/organellar/mito/FILE_Mito_mapped_Pi.csv"))
Plastid <- F_ST.stats(Plastid)
Plastid_pi <- Plastid@nuc.diversity.within/Plastid@n.sites
write.csv2(Plastid_pi, file=file.path(here::here(), "results/organellar/plastid/FILE_Plastid_mapped_Pi.csv"))
#done