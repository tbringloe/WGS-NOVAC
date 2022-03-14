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
Data: paired-end short read data with the following formatting: 01-SampleID_Population/species_1/2.fq.gz, 02_SampleID_Population/species_1/2.fq.gz, 03_SampleID_Population/species_1/2.fq.gz...; a reference genome stored as a fasta file; mitochondrial and plastid seed fasta files (i.e. a gene in the target, or a closely related organism) stored in NOVOPlasty/seed_files; optional BED file containing repeat regions of the reference genome, to be filtered from the VCF files, stored in vcf_analyses/repeats.

Dependencies: fastqc v0.11.9, multiQC v1.9, trimmomatic v0.39, NOVOPlasty v4.2, seqtk v1.3, perl v.5.34.0, bowtie2 v2.4.2, samtools v1.13, bcftools v1.12, vcftools v0.1.16, r v4.0.4, admixture v1.3.0, plink v2.00, python v3.9.5, 

R packages: Here, tidyverse, rprojroot, ggplot2, optparse, PopGenome

## Caveats
I find standard practices and justifications for filtering parameters related to nuclear wide variants basically intractable. Some make sense, like having decent read depth, but others, such as minor allele frequency, appear to be born out of an adundance of caution, and may not be justifiable to impliment from a biological perspective. See my thoughts below with specific parameters in Uncertainty section. Because these decisions are individualized across datasets, some results are not directly comparable across analyses. As in, one cannot compare results A with results B; ideally, the two datasets must be analysed together, under the same conditions, as a single workflow. For instance, exact values of inbreeding coefficient cannot (and should not) be forwarded as an objective truth comparable to other measures, since those values depend, in part, on the filtering choices made. Results should always be interpreted within the context of the data analysed and decisions made during that process.


## Uncertainty
The workflow has been tested using the exact versions in the above dependencies. No effort has been made to test other versions, or enhance compatibility of dependencies (modules are purged and loaded as needed). This needs to be managed, and is an area of active work (e.g. exploring bringing workflow into snakemake). I cannot guarantee the command line arguments will work with different versions of the above programs (in some cases, such as bcftools, I can confirm they won't). Be wary using different versions with unknown behavior, particularly for the variant filtering procedures.

The workflow is meant to automate and dampen learning curves as much as possible for the user, but important decisions must be made regarding filtering criteria for genome-wide nuclear variants. Some guidance from personal experience are provided here for the user to consider, but it is no substitute for exploring individual datasets and making informed decisions regarding these parameters.

##User inputs following filtering settings for compiling VCF file## Presets and justification ##
min_cov=15 # 15 is a good target for minimum read coverage for calling a SNP. There is a certain percentage chance the base call for a read will be incorrect. In fact, across all the reads mapped to an entire genome, there should be a lot of those artifacts presents in the raw VCF file. Stacking reads snuffs out the chance of calling a SNP when in fact it was just a bad base call on the read. Imaging stacking a 1/100 chance of an incorrect base call; calling an incorrect SNP evaporates at higher read depth. In fact, it may be justifiable in some cases to call bases with far less reads than 15, I've seen as little as 5-8 in the literature. There is a movement for low coverage genotyping approaches because it saves a lot on sequencing costs (literally 1/3 of the costs if targetting 5x vs 15x). HOWEVER, I still advise at least 15x coverage. That sequencing effort is halved when calling heterozygous sites (i.e. two alleles). It's a garbage in-garbage out philosophy with this parameter, and my opinion is one can severely compromise certainty in variant positions by skimping on coverage. This is arguably the most important parameter to decide on, and should be made before any sequencing is done, never in retrospect to salvage a dataset. For anyone working with large genomes, I feel your pain.
max_cov=100
allelic_balance_low=0.2
allelic_balance_high=5
minor_allele_frequency=0.02
min_Q=30 #follows phred scale, 30=1/1000 chance of SNP calling error
max_missing=0.9 # value is proportion of present data needed to keep a SNP site, so 0.9=10% missingness

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

