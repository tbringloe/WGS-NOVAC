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
- [Contents](#contents)
  + [Subsections within contents](#subsections-within-contents)
- [Methods](#methods)
  + [Subsections within methods](#subsections-within-methods)
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


## Status
In-development


## Contents
Describe the contents of the repository. Are there multiple scripts or directories? What are there purpose and how do they relate to each other?
The workflow is currently contained in a slurm script, which is specific to queing jobs in shared high performance computing environments. The command lines can, however, be run on a linux-based private server.

Several directories are also provided. The script should be run from the root folder. The R scripts leverage the "Here" package, which swaps in your system file paths so as to avoid hang ups on novel systems. The various folder pathways have various purposes, such as store workflow scripts, hold temporary folders, final output for various analyses, ect. Details are provided below.

### Subsections within contents
IDX/-index folder for the reference genome, used by bowtie2 for mapping reads.

NOVOPlasty/ - contains subfolders and NOVOPlasty perl script that assemble organellar genomes (not my script, see dependencies below). The main workflow script will create new subfolders within the relevant folder stream (mito or plastid) to store input and output files for individual sample assemblies.

R_code/ - contains several R scripts that generate plots at various workflow checkpoints. These include vcfQC scripts, which plot aspects of the raw and filtered VCF files the user should scrutinize for quality control purposes, including distributions in quality score, read depth/site and /individual, missingness/site and /individual, minor allele frequencies, and heterozygosity. R scripts are also provided to plot admixture (provided by GitHub user joanam; speciationgenomics/scripts) and PCA results.

read_files/ - read files are deposited in read_files/raw.

reference_genome/ - a reference genome is stored here

sorted_bam/ - sorted BAM files are stored here during the workflow

tmp/ - some temporary files are stored here during the workflow before being deleted

vcf/ - final vcf files are stored here

vcf_analyses/ - contains several subfolders where analytical output is stored

WGS_NOVAC_14iii22.slurm - the main workflow slurm script submitted via sbatch in an HPC environment


## Methods
What methods were used to achieve the purpose? This should be detailed as possible.
### Subsections within methods
Often useful to organise the methods under multiple subheadings.


## Requirements
Data: paired-end short read data with the following formatting: 01-SampleID_Population/species_1/2.fq.gz, 02_SampleID_Population/species_1/2.fq.gz, 03_SampleID_Population/species_1/2.fq.gz...; it is critical the order of the samples follows the order you wish specimens to appear in the admixture plots, that is, grouped meaningfully by population or species. The naming format allows the workflow to derive species/population metadata, which is leveraged in some of the analyses.

A reference genome stored as a fasta file.

Mitochondrial and plastid seed fasta files (i.e. a gene in the target, or a closely related organism) stored in NOVOPlasty/seed_files

Optional BED file containing repeat regions of the reference genome, to be filtered from the VCF files, stored in vcf_analyses/repeats.

Dependencies: fastqc v0.11.9, multiQC v1.9, trimmomatic v0.39, NOVOPlasty v4.2, seqtk v1.3, perl v.5.34.0, bowtie2 v2.4.2, samtools v1.13, bcftools v1.12, vcftools v0.1.16, r v4.0.4, admixture v1.3.0, plink v2.00, python v3.9.5, 

R packages: Here, tidyverse, rprojroot, ggplot2, optparse, PopGenome

## Caveats
I find standard practices and justifications for filtering parameters related to nuclear wide variants basically intractable. Some make sense, like having decent read depth, but others, such as minor allele frequency, appear to be born out of an adundance of caution, and may not be justifiable from a biological perspective. See my thoughts below with specific parameters in the Uncertainty section. Because these decisions are individualized across datasets, some results are not directly comparable across analyses. As in, one cannot explicitly compare results A with results B; ideally, the two datasets must be analysed together, under the same conditions, as a single workflow. For instance, exact values of inbreeding coefficients cannot (and should not) be forwarded as an objective truth comparable to other measures, since those values depend, in part, on the filtering choices made. They are informative for comparing one species or population relative to another within the dataset. Results should always be interpreted within the context of the data analysed and decisions made during that process.

A lot of analyses that use vcf as input assume sites without a variant position called represent the reference allele. This is a limitation born out of trying to efficiently store genomic information across samples. This likely has important downstream implications, however, which the community has yet to fully come to grips with. Consider all the meaninful biological data that is tossed in an effort to curate the vcf datasets to ensure no artifacts are present. You are not simply tossing that data, you are actually changing those sites back to the reference allele. This is why filtering decisions have such an important impact on vcf dataset, and why no two vcf files are directly comparable (even if samples were mapped to the same reference genome). It's an unfortunate limitation, but until we can store the full genomes, or somehow incorporate uncertainty without recording this as the reference allele, we are stuck. There is a great paper that tries to fix this, but I need to track it down. I will add here later.


## Uncertainty
The workflow has been tested using the exact versions in the above dependencies. No effort has been made to test other versions, or enhance compatibility of dependencies (modules are purged and loaded as needed). This needs to be managed, and is an area of active work (e.g. exploring bringing workflow into snakemake). I cannot guarantee the command line arguments will work with different versions of the above programs (in some cases, such as bcftools, I can confirm they won't). Be wary using different versions with unknown behavior, particularly for the variant filtering procedures.

The workflow is meant to automate and dampen learning curves as much as possible for the user, but important decisions must be made regarding filtering criteria for genome-wide nuclear variants. Some guidance from personal experience are provided here for the user to consider, but it is no substitute for exploring individual datasets and making informed decisions regarding these parameters.

##User inputs following filtering settings for compiling VCF file## Presets and justification, which are not gospel ##

min_cov=15 # 15 is a good target for minimum read coverage for calling a SNP. There is a certain percentage chance the base call for a read will be incorrect. In fact, across all the reads mapped to an entire genome, there should be a lot of those artifacts presents in the raw VCF file. Stacking reads snuffs out the chance of calling a SNP when in fact it was just a bad base call on the read. Imagine stacking a 1/100 chance of an incorrect base call; calling an incorrect SNP evaporates at higher read depth. It may be justifiable in some cases to call bases with far less reads than 15, I've seen as little as 5-8 in the literature. There is a movement for low coverage genotyping approaches because it saves a lot on sequencing costs (literally 1/3 of the costs if targetting 5x vs 15x). HOWEVER, I still advise at least 15x coverage. Remember, that sequencing effort is halved when calling heterozygous sites (i.e. two alleles). It's a garbage-in garbage-out philosophy with this parameter, and my opinion is one can severely compromise confidence in calling variant positions by skimping on coverage. This is arguably the most important parameter to decide on, and should be made before any sequencing is done, never in retrospect to salvage a dataset. If you are going for lower coverage, do the math and present a compelling case. For anyone working with large charismatic genomes, I'm sorry, you did this to yourself. 

max_cov=100 # There is a lot more freedom with this parameter. My suggestion is to evaluate coverage across samples (site coverage and individual coverage) to inform an upper mean threshold present in the read dataset. If a handful of samples deviate a lot from the mean, you can always subsample reads from those samples. By getting a reasonably even distribution of coverage across samples, you can set this threshold to strategically eliminate high copy regions of the genome (which could lead to mapping errors and spurious variant calls), and the abundant organellar genomes (if these have not already been removed from the reference genome provided). Part of the strength of the pipeline is teasing apart genomic signal across nuclear and organellar compartments, so make sure you are not conflating those signals (though I suspect organellar signal would be utterly dwarfed by the nuclear signal if mapping to a complete genome; nonetheless, it is good practice to keep these tidy and seperate). You can evaluate organellar coverage in the NOVOPlasty output log; remember, the full read datasets are subsampled prior to organellar assembly, so the raw coverage provided in the log files should be scaled up accordingly. 

allelic_balance_low=0.2 # Plenty of heterozygous sites will be called if mapping reads from diploid specimens to a nuclear genome. The confidence of calling a heterozygous site depends, in part, on how balanced the reads are between the two alleles. The expectation is that the reads should be evenly divided between the two alleles. Anything highly skewed from this 50:50 distribution could indicate contamination of the sample. The portion of the pipeline that calls variant positions also integrates a lot on information into that call; there is supposed to be a concensus call based on the information in the mapped reads. At some point, the algorithm must tip in favour of a heterozygous call; hopefully that mostly coincides with biologically relevant information, but sometimes it will just be random because we are dealing with millions of potential variant sites. The user must make a decision, however informed, of where that threshold lies for when reads are balanced enough to call a heterozygous site. You don't want to be too punishing, as you will remove meaningful sites if you are conservative, but you want to eliminate artifacts. I think the 1:5 ratio achieves a good balance here. Remember, this is where you will be thankful for the generous read coverage; allelic balance will be prone to more randomness at lower coverage, making it harder to retain sites with confidence. 

allelic_balance_high=5 # Same justification as above. This simply refers to skewed ratios favouring the alternative heterozygous allele.

minor_allele_frequency=0.02 # This is potentially the greatest source of contention when filtering variant positions. If calling artifacts are entering your dataset, they will likely be at extremely low frequency across all samples, and are probably one off spurious calls. Put differently, biologically meaningful alleles will likely appear more than once, which leans into replication to parse out fact from fiction. For instance, a single spurious allele call at a site among 24 diploid individuals works out to an allele frequency of 0.02 for that one allele (allele occurs 1 in a diploid organism, out of potentially 48 [2x24] alleles). The idea here is to cull one offs as probable calling errors. You should calculate what the threshold of minor allele frequency retention based on that one off rule. If you have 114 diploid organism, that frequency has to be 0.004 to occur once among 228 potential alleles. Of course, all new mutations start at a frequency of 1 among a pool of alleles, that is, an extremely low minor allele frequency. A small proportion of these alleles will go on to fixation in species/population; most will be lost. Regardless of the outcome,  there should always be a large standing pool of extremely rare alleles in your datasets. This is why the wofkflow divides into two streams, one with this parameter applied, the other without, and why the basic population statistics are applied to the VCF file withouth a minor allele frequency filter applied. It doesn't make sense to me to cull genomic diversity before measuring it, and with the other filtering controls in place, I'm confident datasets with these low frequency alleles can be safely interpreted.

min_Q=30 # This score follows a phred scale, such that scores of 30=1/1000 chance of SNP calling error. That means, by keeping sites with a Q score of at least 30, you accept that 1 in every 1000 variant sites could be flat out wrong. When calling millions of sites, that could mean thousands of incorrect variant positions. You can view the distribution in scores as part of the vcfQC output. Most of the scores likely far exceed 30 (the upper limit is 999, which would work out to a vanishingly small chance of an error). Really, this is just a step to weed out really poor calls. What goes into the Q score? Not entirely sure, when I tried to track this down it appeared to factor in read quality information, mapping scores, and I assume read depth and consistency is also factored in.

max_missing=0.9 # This value represents the proportion of present data needed to keep a SNP site, so 0.9=10% missingness. You can thank vcftools for that backwards logic. There is no clear threshold for how much missing data should be tolerated before a site is tossed. As with all these parameters, it's a balance between data retention and confidence in that data. My gut feeling is that 20% (0.8) is an acceptable value, but the user should probably explore several thresholds to see what impact, if any, this has on results.

## Acknowledgements
All the inspirational students, postdocs, mentors, and forum junkies across the globe who contributed to my own learning journey in bioinformatics, of which this workflow is a direct result.

## References
Please see GitHub pages here for:

NOVOPlasty: https://github.com/ndierckx/NOVOPlasty

Scripts provided by joanam: https://github.com/speciationgenomics/scripts

Very useful guide for filtering and visualizing SNPs (from which some of this workflow is derived): https://speciationgenomics.github.io/filtering_vcfs/ 

Mermaid wofkflow diagram: https://mermaid.live/

## Workflow diagram
![Workflow_mermaid](https://user-images.githubusercontent.com/79611349/158207461-b82bde6c-e6e2-4b1d-92a2-6e225bf78519.jpg)

