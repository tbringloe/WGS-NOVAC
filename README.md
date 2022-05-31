# Whole Genome Sequencing Nuclear and Organellar Variant Caller for Biodiversity Discovery

__Main author:__  Trevor T. Bringloe  
__Affiliation:__  University of Melbourne
__2nd Affiliation:__  Fisheries and Oceans Canada (DFO)   
__Group:__        NA   
__Location:__     New Brunswick, Canada  
__Affiliated publication:__  
__Contact:__      e-mail: trevor.bringloe@unimelb.edu.au | tel: (506)-259-2288


- [Objective](#objective)
- [Summary](#summary)
- [Status](#status)
- [Example results](#example-results)
- [Contents](#contents)
  + [Provided files](#provided-files)
- [Methods](#methods)
  + [Raw read QC and organellar assembly](#raw-read-qc-and-organellar-assembly)
  + [Read mapping to nuclear genome and VCF compilation](#read-mapping-to-nuclear-genome-and-vcf-compilation)
  + [Phylogenetic and population genomic analyses](phylogenetic-and-population-genomic-analyses)
- [Requirements](#requirements)
- [Caveats](#caveats)
- [Uncertainty](#uncertainty)
- [Acknowledgements](#acknowledgements)
- [References](#references)


## Objective
Make whole genome sequencing a new standard in biology for species discovery, and phylogenetic/population inferences.

This workflow will assemble organellar genomes and call genome-wide nuclear variant positions, perform filtering to retain high quality variant positions, output plots to assess the downstream effects of filtering criteria, and execute basic phylogenomic and population genomic analyses. The user specifies key criteria upfront, which can be adjusted for multiple runs of the workflow, ensuring data exploration and reproducibility.


## Summary
The age of DNA barcoding has transformed the field of biology, revealing remarkable levels of cryptic diversity, unexpected phylogeographic distributions, and novel evolutionary insights. Sequencing whole genomes, that is all the genomic information present in a set of specimens/samples, represents the next major step in species detection and inferring evolutionary relationships, yet only a handful of studies have leveraged this approach for these purposes. Steep learning curves and standardization of bioinformatic workflows present major barriers for the uptake of whole genome sequencing. Here, we introduce a workflow that inputs sequence data and a reference genome, and outputs fully assembled mitochondrial and chloroplast genomes (for photosynthetic organisms), and a file of nuclear-wide variant positions across samples (Single Nucleotide Polymorphisms and indels). The workflow operates through the  command-line interface, wherein users specify a set of key parameters related to assembly, read mapping, and quality control for the retention of nuclear variant positions (e.g. min/max coverage, r2 values for linkage disequilibrium, ect). Among the key outputs are organellar and nuclear phylogenies, a phylogenetic network of nuclear variant positions, admixture and PCA plots, and basic population statistics such as levels of diversity and inbreeding coefficients. The workflow has been used to detect 100,000s to 1,000,000s of variant positions in global kelp datasets, several orders of magnitude greater than current sanger sequencing and reduced genomic capture methods (e.g. RADseq). The workflow has revealed species and population level insights, including widespread hybridizations among species, a novel Arctic lineage, and high latitude glacial refugial populations in the North Atlantic. Facilitating the accessibility of bioinformatic workflows will be imperative to the transition to whole genome sequencing in biology, which in turn promises to reveal new species, holobiome associations, and functional insights.

## Example results
Among the provided files is an html that details results from the workflow analysing Gulf of Alaska kelp populations. This is not published work, it is meant to guide user expectations in terms of what the workflow is capable of. Note, some additional analyses (LD plots) are not yet part of the workflow, and figures were modified using inkscape (i.e. this does not represent raw output).

## Status
**In-development**
**added May 26, 2022: new GQ filtering option** The basic workflow is working, but the user is still expected to have enough expertise to troubleshoot and customize the workflow for their own datasets. Future enhancements will potentially include: dependency management, likely by bringing the workflow into snakemake; automated annotation of organellar genomes; phasing of nuclear variant positions; stream for long read data input; incorporate LD decay plots to inform plink parameter decisions; add reference genome assembly workflow; correction for variants mapped to paralogous genes.


## Contents
The workflow is currently contained in a slurm script, which is submitted to a shared high performance computing environment for job queueing. The command lines can, however, be run on a linux-based private server.

Several directories are also provided to manage files during the run. The script should be run from the root folder. The R scripts leverage the "Here" package, which swaps in your system file paths so as to avoid hang ups on novel systems. The various folder pathways have various purposes, such as store workflow scripts, hold temporary folders, final output for various analyses, ect. Details are provided below.

### Provided files
**01-07.XXX.R** - these are R scripts used by the workflow to plot results

**READme.txt** - contains important instructions regarding how to setup the environment prior to analysis

**config_mito.txt & config_plastid.txt** - these are modified config files used by NOVOPlasty for de novo organellar genome assembly

**TAXA_BLOCK** - a template used to modify nexus files before analysis by splitstree (if the user is using this option)

**to_exe_XXX** - configuration files also used by splitstree if the user is running this in the HPC environment

**WGS_NOVAC_20iv22.slurm** - the main workflow slurm script submitted via sbatch in an HPC environment

**07_Rmarkdown_GoA_prelim_summary_20iv22.html** - example of results produced by the workflow. Not published nor peer reviewed


## Methods
The workflow automates the compilation of variant positions across genomic compartments, leveraging preexisting programs. The user must open the slurm script and establish all the variable settings outlined at the beginning of the workflow. The HPC job should also be configured to match some of these setting, i.e. number of threads requested. A modest amount of memory should also be requested depending on the scale of the datasets, and can be increased if memory runs out and the job is killed. Below are example inputs that appear at the beginning of the slurm script.

```
################################################User inputs#############################################################
##General information for running workflow
threads=32 #be sure to specify this in the slurm script
EXAMPLE_SAMPLE=01_HP-9_HalibutPoint #put an example sample name here for nuclear variant compilation step (uses this sample to extract contig names)
FILE=GoA_test #an informative prefix to carry through analyses

##User inputs following file settings
reference_nuclear_genome=alaria_v3.fasta
reference_mitochondrial_genome=Alaria_marginata_HP-3_M_GoA.fasta
reference_chloroplast_genome=Alaria_marginata_HP-3_P_GoA.fasta
repeat_regions=Alaria_repeats.bed # optional file to specify repeat regions to exclude from VCF file
mito_seed=Alaria_COI.fasta
plastid_seed=Alaria_rbcL.fasta
ploidy=2 # specify ploidy of the read datasets for the mpileup step of the nuclear variant positions; if the dataset is a mix of haploid and diploid individuals, you will need to specify this manually in the samples_list file used at the bcftools mpileup step

##User inputs settings for NOVOPlasty
Number_reads_to_extract=5000000 #because organellar genomes are typically high copy, we need a fraction of total read data. Try specifying more reads if genome does not assemble as a circular contig
Genome_range_mito=35000-45000 #range in kbp
Genome_range_plastid=125000-140000 #range in kbp
Kmer=55 #you can try increasing this value if circular genome is not produced and coverage is high
Read_length=150
Insert_size=200

##Trimmomatic settings
Trailing=10
Headcrop=15
Average_quality=20
Min_length=75

##Mapping parameters for bowtie2; 0.6=up to 10% divergence in mapping high quality reads in end-to-end mode, 0.3=up to 5%, 0.12=up to 2%; bowtie2 manual states minimum-score function f to f(x) = 0 + -0.6 * x, where x is read length.
map_param_nuc=0.6
map_param_mito=0.3
map_param_chloro=0.12

##User inputs following filtering settings for compiling nuclear VCF file
min_cov=15
max_cov=100
allelic_balance_low=0.2
allelic_balance_high=5
minor_allele_frequency=0.02
min_GQ=30 #filters individual genotypes according to their probability of being correct, follows a phred scale, 30=1/1000 chance of SNP calling error
SNPgap=1 #filters genotypes adjacent to indels, which are enriched with errors (indels are removed)
min_Q=30 #follows phred scale, 30=1/1000 chance of SNP calling error
max_missing=0.9 # value is proportion of present data needed to keep a SNP site, so 0.9=10% missingness

##User specifies following parameters for removal of linked SNPs using PLINK
plink_r2=0.25
plink_window_size=25 #in kb

##User specifies values of k for ADMIXTURE analysis
k_low=2
k_high=3
pops=HalibutPoint,KayakBeach #here, list populations in order they should appear in admixture plot, seperated by a comma

##User specifies the window and step sizes in bp for vcftools pop gen analyses
window=50000 # in bp; it is recommended to set this high; vcftools will output results for windows with variant(s) present, upward biasing estimates due to excluded 0 value windows. Setting a larger window will minimize this bias.
step_size=5000 # in bp
```

Once the user is happy with these choices, the slurm script must be submitted to an HPC server. Alternatively, the individual line codes can be run on a private, linux-based server. If the user wants to run specific sections, avoid rerunning, or remove particular commmands from the script, they can put a # at the beginning of relevant lines, or delete altogether (providing they are not necessary). It is recommended to break up the job, running computationally intensive portions (calling variant positions) with many threads before running seperate jobs for the analyses (which generally do not use multiple threads).

```
sbatch WGS_NOVAC_20iv22.slurm
```

The workflow, in terms of the methods leveraged, are described here. For further details on methodology and benchmarking of the specific programs used, see [References](#references), which include links to GitHub/webpages.
### Raw read QC and organellar assembly
The user must first provide a single column text file called sample.list in the workflow root directory containing the sample IDs (i.e. 01_SampleID_Pop/Species) used as the prefix in the fastq files.

The pipeline starts by performing fastqc/multiqc on the raw read files. The files are then trimmed according to user specified criteria. Fastqc and multiqc are performed again on the trimmed reads, and should be checked over and compared to the raw QC output. Once the reads are trimmed, the workflow will move onto organellar variant calling using a read mapping approach, but the user can also optionally specify de novo assembly using NOVOPlasty (Dierckxsens et al. 2017); this program extends a user provided seed sequence based on coverage (i.e. high frequency kmers are used to extend the seed), until the assembled genome begins to overlap, confirming circularity. Seperate folders and configuration files are created for each sample. Carefully check over the final contig(s), if the assembly is not complete, this can potentially be improved two ways: 1) specify more reads to subsample for organellar assembly (improving coverage), or 2) specify a higher Kmer value, which may help bridge breaks in the assembly (you might have to increase the number of subreads as this will decrease kmer coverage). Performance will depend on the nature of the organellar genome. Complex samples with lots of different genomes, or repeat patterns, will lead to breaks in the assembly. This can create a lot of manual work trying to piece the organellar genome together. The former option (mapping approach) is suggested for downstream analyses, the latter (de novo assembly) is suggested for publishing organellar genomes.

Providing organellar assembly was successful, the genomes can be aligned using MAUVE (Darling et al. 2004). Phylogenetic analyses are performed using the read mapping results, and using RAxML (Stamatakis et al. 2014). Basic population statistics are also performed on the read mapping results, using PopGenome (Pfeifer et al. 2014).

### Read mapping to nuclear genome and VCF compilation
The trimmed reads are mapped to the reference genome using end-to-end mode in bowtie2 (Langmead and Salzberg 2012); an alternative version of the workflow is being developed to input long read data. A single bowtie2 parameter is specified by the user, which indicates a cutoff in percent divergence to map to the reference genome. For phylogenetic work mapping across multiple species, this parameter should be set high (perhaps as high as 20% divergence or more, depending on how distantly related the species are expected to be relative to the reference genome), for population level work, this parameter should be set low (2-5% divergence). The % divergence is approximate because the bowtie2 parameter being set is dynamic, since it factors in read quality information (poor base calls count less towards whether the read is mapped or not). Please refer to the bowtie2 manual to see exactly how mapping scores are calculated.

The workflow will loop through mapping through all the samples. It will also convert the resulting SAM files to BAM using samtools (Li et al. 2009), and sort the files for compilation using bcftools (Danacek et al. 2021). At this point, a clever one liner a colleagued pointed out to me is used to call variant positions across samples on a contig by contig basis (accessed here: https://www.ecseq.com/support/ngs-snippets/how-to-run-time-consuming-data-analysis-processes-in-parallel-on-unix-systems). This allows the workflow to compile VCF files in parallel (a VCF file each for multiple contigs at the same time, rather than one by one in a single file), which frees up what was orignally a significant chokepoint in the workflow. Once VCF files are compiled for all the contigs, they are merged into a single VCF, and all the temporary files are deleted.

The VCF file at this point represents the raw dataset. The workflow will generate several figures to help evaluate QC of the raw VCF file. The raw VCF is then filtered according to user specified criteria using a combination of bcftools (Danacek et al. 2021) and vcftools (Danacek et al. 2011). I have provided some guidance on how to set filtering parameters in the [Uncertainty](#uncertainty) section below. The workflow is divided into two streams, one with a minor allele frequency filter applied (justified in the Uncertainty section), the other without. Another set of QC figures are produced for each stream, and should be carefully studied and compared to the raw vcfQC files. These vcf files represent checkpoints with accompanying fasta files retained for the user.

A final filter is applied for linkage disequilibrium using Plink (Purcell et al. 2007). This too represents a checkpoint where a new vcf and fasta file is created. Several analyses are conducted on the various vcf files produced.

### Phylogenetic and population genomic analyses
Basic population statisitics are performed on the VCF file not filtered for minor allele frequency. The statistics are performed through vcftools (Danacek et al. 2011), and include the use of a sliding window to calculate nucleotide diversity (pi), and Tajima's D. Inbreeding coefficients (heterozygosity) is also calculated for all individuals. For datasets compiled with invariant positions, pixy is used for distance calculations (Korunes and Samuk 2021); note, preferred method as this corrects for bias introduced by missingness). If the analysis is at the species level, this output should be ignored. Phylogenetic analyses are performed on the VCF file filtered for linkage disequilibrium (LD). These analyses include an initial calculation of eigenvectors and eigenvalues using Plink, with a corresponding PCA plotted using R. Admixture analyses (Alexander and Lange 2011) also leverage the Plink output, and an R script plots ancestry proportions for all individuals at the various levels of k specified by the user. A fasta file of the LD pruned dataset is analysed using RAxML to produce ML tree output. 

Splitstree (Huson and Bryant 2006) can be used to plot uncorrected genomic distances as a network, which can be further scrutinized for shared genomic information across samples (i.e. potential hybridizations). As the program requires a GPU, it is recommended to run this on a standard computer.

There is optional command strings the user can modify to derive coverage information across genomic regions specified by a bed file(s). This could be useful if the user is interested in whether exons have been lost or expanded in particular genomes.

## Requirements
**Data: paired-end short read data with the following formatting: 01-SampleID_Population/species_1/2.fq.gz, 02_SampleID_Population/species_1/2.fq.gz, 03_SampleID_Population/species_1/2.fq.gz...;** it is critical the order of the samples follows the order you wish specimens to appear in the admixture plots, that is, grouped meaningfully by population or species. The naming format allows the workflow to derive species/population metadata, which is leveraged in some of the analyses.

**A reference genome** stored as a fasta file.

**Mitochondrial and plastid seed fasta files** (i.e. a gene in the target, or a closely related organism) stored in NOVOPlasty/seed_files

**Optional BED file containing repeat regions of the reference genome**, to be filtered from the VCF files, stored in vcf_analyses/repeats.

**Program dependencies:** fastqc v0.11.9, multiQC v1.9, trimmomatic v0.39, NOVOPlasty v4.2, seqtk v1.3, perl v.5.34.0, bowtie2 v2.4.2, samtools v1.13, bcftools v1.12, vcftools v0.1.16, r v4.0.4, admixture v1.3.0, plink v2.00, python v3.9.5, 

**R packages:** Here, tidyverse, rprojroot, ggplot2, optparse, PopGenome, gplots

## Caveats
The fasta files of the organellar sequences based on read mapping **should not** be published. Missing data will be treated as the reference allele, leading to errors in the sequence (these are accepted for downstream analysis, as mapping allows for efficient alignment of variant positions). De novo assembled genomes should be considered for publication purposes.

Nucleotide diversity and genetic distances in the nuclear datasets are severely biased because missing data is improperly treated as the reference allele by vcftools. VCF files with invariant position data should be used to calculate nucleotide diversity and genetic distances using Pixy (Korunes and Samuk 2021). This is currently implimented into the workflow, but the user must modify the slurm script to enable this stream (mpileup and call step).

I find standard practices and justifications for filtering parameters related to nuclear wide variants basically intractable. Some make sense, like having decent read depth, but others, such as minor allele frequency, appear to be born out of an adundance of caution, and may not be justifiable from a biological perspective. See my thoughts below with specific parameters in the Uncertainty section. Because these decisions are individualized across datasets, some results are not directly comparable across analyses. As in, one cannot explicitly compare results A with results B; ideally, the two datasets must be analysed together, under the same conditions, as a single workflow. For instance, exact values of inbreeding coefficients cannot (and should not) be forwarded as an objective truth comparable to other measures, since those values depend, in part, on the filtering choices made. They are informative for comparing one species or population relative to another within the dataset. Results should always be interpreted within the context of the data analysed and decisions made during that process.

A lot of analyses that use vcf as input assume sites without a variant position called represent the allele of the reference genome. This is a limitation born out of trying to efficiently store genomic information across samples (assumptions are made so information on all positions does not need to be stored or analysed). This has important downstream implications, however, which the community has yet to fully come to grips with. Consider all the meaninful biological data that is tossed in an effort to curate the vcf datasets to ensure no artifacts are present. You are not simply tossing that data, you are actually changing those sites back to the reference allele. This is why filtering decisions have such an important impact on vcf datasets, and why no two vcf files are directly comparable (even if samples were mapped to the same reference genome). It's an unfortunate limitation, but until we can efficiently store the full genomes, or somehow incorporate uncertainty without recording this as the reference allele, we are stuck. There is a great paper that tries to fix this, but I need to track it down. I will add here later.

The quality and reliability of the workflow inherently depends on the quality of the reference genome. The user should ensure the reference genome is as clean as possible from non-target taxa (otherwise you could potentially call SNPs in the non-target taxa, assuming they are present in all samples fed into the workflow). I suggest a blobtools (https://github.com/DRL/blobtools) approach to assessing the purity of the reference genome, and to use coverage, GC content, and taxonomic information to filter non-target contigs. Generating a solid reference genome is a logistical challenge in-and-of itself. I have colleagues with fantastic workflows for genome assembly, and I will attempt to stich those workflows into this one in future enhancements.

The workflow assumes a consistent ploidy level across all individuals. If this is not the case (i.e. datasets are a mixture of haploid and diploid individuals), this must be specified in the sample_list file used by bcftools at the mpileup step of the workflow. This cannot be automated, so the user will have to add ploidy manually to the file as a second tab-delimited column.

The workflow also assumes individuals are not related. If they are, this has critical implications for population genomic analyses. I heatmap of relatedness is among the future additions to the workflow to help users evaluate this aspect of their datasets. Note, the relatedness graphs need to be corrected to analyse individuals within population, rather than across. The current workflow biases values towards higher relatedness by considering all the populations together.

The workflow does not yet correct for variants mapped to paralogous genes, which is expected to artificially boost heterozygosity. For now, user may choose to map to single copy genes as the reference genome. Corrections are currently in development, specifically a filtering step for excess heterozygosity. Command strings are now available to identify regions of gene expansion based on coverage information; if the user identifies putative expanded regions of the genome, these can be filtered using the vcftools --exclude-bed command.

## Uncertainty
The workflow has been tested using the exact versions in the above dependencies. No effort has been made to test other versions, or enhance compatibility of dependencies (modules are purged and loaded as needed). This needs to be managed, and is an area of active work (e.g. exploring bringing workflow into snakemake). I cannot guarantee the command line arguments will work with different versions of the above programs (in some cases, such as bcftools, I can confirm they won't). Be wary using different versions with unknown behavior, particularly for the variant filtering procedures.

The workflow is meant to automate and dampen learning curves as much as possible for the user, but important decisions must be made regarding filtering criteria for genome-wide nuclear variants. Some guidance from personal experience are provided here for the user to consider, but it is no substitute for exploring individual datasets and making informed decisions regarding these parameters.

*User inputs following filtering settings for compiling VCF file - Presets and justification, which are not gospel*

**min_cov=15** # 15 is a good target for minimum read coverage for calling a SNP. There is a certain percentage chance the base call for a read will be incorrect. In fact, across all reads mapped to an entire genome, there should be a lot of those artifacts presents in the raw VCF file. Stacking reads snuffs out the chance of calling a SNP when in fact it was just a bad base call on the read. Imagine stacking a 1/100 chance of an incorrect base call; calling an incorrect SNP evaporates at higher read depth. It may be justifiable in some cases to call bases with far less reads than 15, I've seen as little as 5-8 in the literature. There is a movement for low coverage genotyping approaches because it saves a lot on sequencing costs (literally 1/3 of the costs if targetting 5x vs 15x). HOWEVER, I still advise at least 15x coverage. Remember, that sequencing effort is halved when calling heterozygous sites (i.e. two alleles). It's a garbage-in garbage-out philosophy with this parameter, and my opinion is one can severely compromise confidence in calling variant positions by skimping on coverage. This is arguably the most important parameter to decide on, and should be made before any sequencing is done, never in retrospect to salvage a dataset. If you are going for lower coverage, do the math and present a compelling case. For anyone working with large charismatic genomes, I'm sorry, you did this to yourself. 

**max_cov=100** # There is a lot more freedom with this parameter. My suggestion is to evaluate coverage across samples (site coverage and individual coverage) to inform an upper mean threshold present in the read dataset. If a handful of samples deviate a lot from the mean towards higher coverage values, you can always subsample reads from those samples. By getting a reasonably even distribution of coverage across samples, you can set this threshold to strategically eliminate high copy regions of the genome (which could lead to mapping errors and spurious variant calls; note, an option to specify a BED file is also provided to remove repeat regions for this same reason). Setting an upper coverage threshold also removes sites correspinding to the abundant organellar genomes (if these have not already been removed from the reference genome provided). Part of the strength of the pipeline is teasing apart genomic signal across nuclear and organellar compartments, so make sure you are not conflating those signals (though I suspect organellar signal would be utterly dwarfed by the nuclear signal if mapping to a complete genome; nonetheless, it is good practice to keep these tidy and seperate). You can evaluate organellar coverage in the NOVOPlasty output log; remember, the full read datasets are subsampled prior to organellar assembly, so the raw coverage provided in the log files should be scaled up accordingly. 

**allelic_balance_low=0.2** # Plenty of heterozygous sites will be called if mapping reads from diploid specimens to a nuclear genome. The confidence of calling a heterozygous site depends, in part, on how balanced the reads are between the two alleles. The expectation is that the reads should be evenly divided between the two alleles. Anything highly skewed from this 50:50 distribution could indicate contamination of the sample. The portion of the pipeline that calls variant positions also integrates a lot on information into that call; there is supposed to be a concensus call based on the information in the mapped reads. At some point, the algorithm must tip in favour of a heterozygous call; hopefully that mostly coincides with biologically relevant information, but sometimes it will just be random because we are dealing with millions of potential variant sites. The user must make a decision, however informed, of where that threshold lies for when reads are balanced enough to call a heterozygous site. You don't want to be too punishing, as you will remove meaningful sites if you are conservative, but you want to eliminate artifacts. I think the 1:5 ratio achieves a good balance here. Remember, this is where you will be thankful for the generous read coverage; allelic balance will be prone to more randomness at lower coverage, making it harder to retain sites with confidence. 

**allelic_balance_high=5** # Same justification as above. This simply refers to skewed ratios favouring the alternative heterozygous allele.

**minor_allele_frequency=0.02** # This is potentially the greatest source of contention when filtering variant positions. If calling artifacts are entering your dataset, they will likely be at extremely low frequency across all samples, and are probably one off spurious calls. Put differently, biologically meaningful alleles will likely appear more than once, which leans into replication to parse out fact from fiction. For instance, a single spurious allele call at a site among 24 diploid individuals works out to an allele frequency of 0.02 for that one allele (allele occurs 1 in a diploid organism, out of potentially 48 [2x24] alleles). The idea here is to cull one offs as probable calling errors. You should calculate the threshold of minor allele frequency retention based on that one off rule. If you have 114 diploid organism, that frequency has to be 0.004 to occur once among 228 potential alleles. Of course, all new mutations start at a frequency of 1 among a pool of alleles, that is, an extremely low minor allele frequency. A small proportion of these alleles will go on to fixation in species/population; most will be lost. Regardless of the outcome,  there should always be a large standing pool of extremely rare alleles in your datasets. This is why the wofkflow divides into two streams, one with this parameter applied, the other without, and why the basic population statistics are applied to the VCF file withouth a minor allele frequency filter applied. It doesn't make sense to me to cull genomic diversity before measuring it, and with the other filtering controls in place, I'm confident datasets with these low frequency alleles can be safely interpreted.

**min_Q and min_GQ=30** # These scores are a measure of the precision or accuracy of the variant calls. This score follows a phred scale, such that scores of 30=1/1000 chance of variant calling error. That means, by keeping sites with a Q score of at least 30, you restrict your analysis to the prospect that, at worst, a given genotype/site has a 1 in 1000 chance of being wrong. The GQ parameter works at the individual genotype level (GQ format field in vcf; probability the variant call for a given specimens at a given position is an error), whereas Q works at the site level (QUAL column in vcf; probability of an alternate allele present at a given site). When calling millions of sites, that could mean thousands of incorrect variant positions. You can view the distribution in scores as part of the vcfQC output. Most of the scores likely far exceed 30 (the upper limit is 999, which would work out to a vanishingly small chance of an error). Really, this is just a step to weed out really poor calls, or genotypes with competing probabilities. What goes into the Q score for genotypes? According to Danecek et al. 2021 "(from the mpileup step which summarizes read alignment information for all positions on the reference) Genotype likelihoods are then calculated, representing how consistent is the observed
data with the possible diploid genotypes. The calculation takes into account mapping qualities of the reads, base qualities and the probability of local misalignment, per-base alignment quality (BAQ)". Mapping qualities are recorded by bowtie2, such that is reads map equally well to multiple positions across the reference genome, they are assigned a low mapping score (an excellent deep dive can be found here: https://biofinysics.blogspot.com/2014/05/how-does-bowtie2-assign-mapq-scores.html). The site quality score (QUAL) would be a downstream summary of the likelihood of a non reference allele based on all the calls (and presumably likelihood information) across all the samples. Filtering variants in the dataset for all this information therefore incorporates upstream information on read quality and mapping performance.

**SNPgap** # Indels are expected to be enriched with genotyping errors, and therefore the user may want to filter genotypes adjacent to indels. The value is set to 1, a conservative approach would be to remove SNPs within 3 bp of indels, but this comes at the expense of increased missingness in the dataset.

**max_missing=0.9** # This value represents the proportion of present data needed to keep a SNP site, so 0.9=10% missingness. You can thank vcftools for that backwards logic. There is no clear threshold for how much missing data should be tolerated before a site is tossed. As with all these parameters, it's a balance between data retention and confidence in that data. My gut feeling is that 20% (0.8) is an acceptable value, but the user should probably explore several thresholds to see what impact, if any, this has on results.

**plink_r2=0.15 and plink_window_size=25 # in kb** # Filtering for linkage disequilibrium is a key assumption when evaluating population structure or phylogenetic signal. Linked positions might lead to an uneven distribution of phylogenetic signal across the genome, with accelerated evolution in some areas of the genome leading to "spikes" in signal related to non-random evolutionary events; pruning for linked sites is a way of evening the playing field, and not interpreting non-random evolutionary events (such as selective forces) for random ones (like population structure that accumulates over time through random mutations and genetic drift). Interpreting that "random" signal is crucial for making inferences about population history, because it is shaped by things like population expansion or bottlenecks. Setting an r2 value (i.e. how correlated sites must be to be pruned) is not straightforward, and neither is the window size. The best approach is to plot the decay in linkage disequilibiurm using PopLDdecay (https://github.com/BGI-shenzhen/PopLDdecay), and use this to inform the choices made here. Implimenting LDdecay plots remains a priority for the workflow.

## Acknowledgements
All the inspirational students, postdocs, mentors, and forum junkies across the globe who contributed to my own learning journey in bioinformatics, of which this workflow is a direct result.

## References
**Please see useful pages here:**

Scripts provided by joanam: https://github.com/speciationgenomics/scripts

Scripts for xmfa to fasta conversion by kjolley: https://github.com/kjolley/seq_scripts

Very useful guide for filtering and visualizing SNPs (from which some of this workflow is derived): https://speciationgenomics.github.io/filtering_vcfs/ 

Mermaid wofkflow diagram: https://mermaid.live/

Clever one liner for VCF compilation step: https://www.ecseq.com/support/ngs-snippets/how-to-run-time-consuming-data-analysis-processes-in-parallel-on-unix-systems

Deep dive on bowtie2's mapping quality score: https://biofinysics.blogspot.com/2014/05/how-does-bowtie2-assign-mapq-scores.html

**Citations for programs:**

Alexander, D. H., Lange, K. (2011). Enhancements to the ADMIXTURE algorithm for individual ancestry estimation. BMC Bioinformatics, 12, 1-6. https://dalexander.github.io/admixture/download.html

Danacek, P., Auton, A., Goncalo, A., Albers, C. A., Banks, E., DePristo, M. A., Handsaker, R., Lunter, G., Marth, G., Sherry, S. T., McVean, G., Durbin, R. & 1000 Genomes Project Analysis Group. (2011). The variant call format and VCFtools. Bioinformatics 27: 2156-8. http://vcftools.sourceforge.net/

Danacek, P., Bonfield, J. K., Liddle, J., Marshall, J., Ohan, V., Pollard, M. O., Whitwham. A., Keane, T., McCarthy, S. A., Davies, R. M., & Li, H. (2021). Twelve years of SAMtools and BCFtools. Gigascience, 10, giab008. https://github.com/samtools/bcftools.

Darling A. C., Mau B., Blattner F. R., Perna N.T. 2004. Mauve: multiple alignment of conserved genomic sequence with rearrangements. Genome Research 14 :1394-1403. doi:10.1101/gr.2289704

Dierckxsens, N., Mardulyn, P., & Smits, G. (2017). NOVOPlasty: de novo assembly of organelle genomes from whole genome data. Nucleic Acids Research, 45, e18. https://github.com/ndierckx/NOVOPlasty

Huson, D. H., & Bryant, D. (2006). Application of phylogenetic networks in evolutionary studies. Molecular Biology and Evolution, 23, 254-67. https://software-ab.informatik.uni-tuebingen.de/download/splitstree4/welcome.html

Korunes, K. L., Samuk, K. (2021). Pixy: unbiased estimation of nucleotide diversity and divergence in the precense of missing data. Molecular Ecology Resources. 21: 1359-1368.

Langmead, B., & Salzberg, S. L. (2012). Fast gapped-read alignment with Bowtie 2. Nature Methods, 9, 357-9. http://bowtie-bio.sourceforge.net/bowtie2/index.shtml

Li, H., Handsaker, B., Wysoker, A., Fennell, T., Ruan, J., Homer, N., Marth, G., Abecasis, G., Durbin, R., & 1000 Genome Project Data Processing Subgroup. (2009). The sequence alignment/map format and SAMtools. Bioinformatics 25: 2078-9. http://www.htslib.org/

Pfeifer, B., Wittelsbürger, U., Ramos-Onsins, S. E., & Lercher, M. J. (2014). PopGenome: An efficient swiss army knife for population genomic analyses in R. Molecular Biology and Evolution, 31, 1929-1936.https://cran.r-project.org/web/packages/PopGenome/index.html

Purcell, S., Neale, B., Todd-Brown, K., Thomas, L., Ferreira, M. A. R., Bender, D., Maller, J., Sklar, P., de Bakker, P. I. W., Daly, M. J., & Sham, P. C. (2007). PLINK: a tool set for whole-genome association and population based linkage analyses. American Journal of Human Genetics, 81, 559-575. https://www.cog-genomics.org/plink/

Stamatakis, A. (2014). RAxML version 8: a tool for phylogenetic analysis and post-analysis of large phylogenies. Bioinformatics, 30, 1312-3. https://cme.h-its.org/exelixis/web/software/raxml/


## Workflow diagrams
## Organellar workflow
![Workflow_organellar_v2](https://user-images.githubusercontent.com/79611349/162233966-e421dccf-93b6-4ca9-ad8c-581f22ffe743.jpg)
## Nuclear workflow
![Workflow_nuclear_v2](https://user-images.githubusercontent.com/79611349/162234013-3d141c60-f267-4e48-b4e7-04471024e9bd.jpg)
