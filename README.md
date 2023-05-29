# Whole Genome Sequencing Variant Caller for Biodiversity Discovery

__Main author:__  Trevor T. Bringloe  
__Affiliation:__  Fisheries and Oceans Canada (DFO)   
__Group:__        NA   
__Location:__     New Brunswick, Canada  
__Affiliated publication:__  
__Contact:__      e-mail: tbringloe@gmail.com | tel: (506)-259-2288


- [Objective](#objective)
- [Summary](#summary)
- [Status](#status)
- [Contents](#contents)
  + [Provided files](#provided-files)
- [Methods](#methods)
  + [Raw read QC](#raw-read-qc)
  + [Organellar assembly](#organellar-assembly)
  + [Read mapping to nuclear genome and VCF compilation](#read-mapping-to-nuclear-genome-and-vcf-compilation)
  + [Genomic analyses to infer structure](genomic-analyses-to-infer-structure)
- [Requirements](#requirements)
- [Caveats](#caveats)
- [Uncertainty](#uncertainty)
- [Acknowledgements](#acknowledgements)
- [References](#references)


## Objective
Make whole genome sequencing a new standard in biology for species discovery, and phylogenetic/population inferences.

This workflow will assemble organellar genomes and call genome-wide nuclear variant positions, perform filtering to retain high quality variant positions, output plots to assess the downstream effects of filtering criteria, and execute basic phylogenomic and population genomic analyses. The user specifies key criteria and workflow pathways upfront based on data types and objectives, which can be adjusted for multiple runs of the workflow, ensuring data exploration and reproducibility.


## Summary
The age of DNA barcoding has transformed the field of biology, revealing remarkable levels of cryptic diversity, unexpected phylogeographic distributions, and novel evolutionary insights. Sequencing whole genomes, that is all the genomic information present in a set of specimens/samples, represents the next major step in species detection and inferring evolutionary relationships. Here, we refined a preexissting workflow that inputs sequence data and a reference genome, and outputs fully assembled mitochondrial and chloroplast genomes (for photosynthetic organisms), and a file of nuclear-wide variant positions across samples (Single Nucleotide Polymorphisms). The workflow is a bash script with forks for low vs high coverage datasets, and variant only or variant+invariant datasets. Users establish the environment(i.e. executables) and specify variables before launching the script. The workflow will also stop at key steps to ensure user intervention for quality control purposes. Among the key outputs are organellar genomic sequences (assembled in parallel for a set of samples), PCA plots, and admixture results. The workflow has been used to detect 100,000s to 1,000,000s of variant positions in global datasets in mammalian populations, but an be generally applied to any organism. Our forthcoming work shows high correspondence with ddRADseq datasets, but with greater fidelity to detect genetic clusters.

## Status
**Development and testing completed**
The workflow is refined from earlier work, with newly added forks to the analysis for low and high coverage options, and invariant+variant datasets. The user is still expected to have enough expertise to troubleshoot and customize the workflow for their own datasets. Future improvements will reflect user feedback.

## Contents
The workflow is currently contained in a bash, which can be run on any linux environment on a private server. Adjustments can be made to to run the script on shared clusters with job queueing; the user would simply submit the job as a slurm script loading required modules and running the command 'bash MOBELS_lcWGS_workflow_25v23_master.sh'. The bash script also contains instructions on other dependencies (i.e. R packages) and file/folder structure to establish prior to running.

### Provided files
**01_vcfQC_RAW.R** - R script for plotting vcf QC prior to filtering

**02_vcfQC_FILTERED.R** - R script for plotting vcf QC post filtering

**03_PCA_plots.R** - R script for plotting PCA results

**04_admixture_plots.R** - R script for plotting admixture results, reproduced and modified with permission from joanam: https://github.com/speciationgenomics/scripts

**ill_adap_uni.fa** - optional fasta file of illumina universal adapters for filtering during read trimming

**lcWGS_config_mito.txt** - optional configuration file fo assembling mitochondrial genomes with NOVOPlasty

**lcWGS_config_plastid.txt** - optional configuration file for assembling chloroplast genomes with NOVOPlasty

**MOBELS_lcWGS_workflow_25v23_master.sh** - bash script to execute the workflow

**MOBELS_plots.R** - alternative R script to generate PCA and admixture plots, courtesy of Audrey Bourret, also available here: XXX.


## Methods
The workflow automates the compilation of variant positions across genomic compartments, leveraging preexisting programs. The user must open the bash script, follow instructions on establishing the environment, and set all variables outlined at the beginning of the workflow. Below are the instructions and disclaimer that appear at the beginning of the slurm script.

```
###########################################################################
############################  Instructions  ###############################
###########################################################################
##Note there are several aspects of the workflow users will need to tailor or consider for things to run smoothly, namely:
##1. User must specify presets listed below, including executable pathways if needed
##2. The user must establish the following naming convention: OrderToAppearInAdmixturePlots.GeographicIdentifier(PCAlabels).SampleYear(optional).Species(optional).ExtractionID(optional).SampleID.MappedGenome(optional, dl=delphinapterusleucas).n_s.bam. Example: 7.RES.2007.BELUGA.ADN_20_05497.S_20_04096.dl.n_s.bam
##3. Several files must be preped before hand, including a reference genome and a bed file of repeat regions to exclude deposited in a folder reference_genome/ created in the launching directory; user can specify to no remove sites from repeat regions, in which case a bed file is not needed
##4. Deposit raw reads and ill_adap_uni.fa into a newly generated raw_reads folder with suffix convention: <sampleID>_R1.fastq.gz; <sampleID>_R2.fastq.gz
##5. The workflow will not work in one run, you need to inspect read quality and rerun the script after specifying you are happy to proceed, see presets below
##6. If assembling organellar genomes, users will need to download NOVOPlasty: https://github.com/ndierckx/NOVOPlasty and deposit in launching directory (files should appear in folder NOVOPlasty-master/)
##7. Change permissions of NOVOPlasty/ files to ensure you can execute the perl script
##8. Users must also create the folder NOVOPlasty-master/seed_files where they should deposit seed files.
##9. Users must also deposit provided NOVOPlasty config file lcWGS_config_mito.txt and lcWGS_config_plastid.txt into NOVOPlasty-master/
#10. Users should inspect the multiqc mapped reads report in results/bowtie2_output and ensure they are satisfied with results
#11. Users must make an R_code folder where they should deposit provided R scripts
#12. There are fringe cases where the workflow will break. For instance, user intervention is needed if datasets represent a mix of haploid and diploid individuals (see ploidy below)
#13. Users may decide to compile a full sites vcf. Fair warning, files may be enormous if working with a large (i.e. >1Gbp) genome. Users could specify to compile based on a list of chromosomes/scaffolds/contigs in this case, see below
#14. Users may decide to compile the vcf for specified chromosome/scaffolds/contigs and whether to concatenate these results. User must provide a file seqs_list with one chromosome/scaffold/contig per line to compile
#15. A variant site filter for heterozygosity is not implemented here. This should be considered when mapping to a non-target species or mapping different species to a single reference genome (i.e. sample specific duplications mapping to same non-duplicated region of reference genome inflate heterozygosity) 
#16. Users will have to investigate variant position filtering and modify script as needed. Explanations for justifications are provided throughout, see presets and command line arguments
#17. The option to filter individuals from the vcf is not provided here. Users may consider removing low quality samples using: bcftools view <input>.vcf.gz -s^ <sample1,sample2,sample3...> --threads <threads> -Oz -o vcf/<output>_raw.vcf.gz. See bcftools manual, a list can be provided if a large number of samples are to be excluded (-S <list_file, one sample per line>). Ensure the output replaces "$FILE"_raw.vcf.gz before restarting the workflow at the vcf filtering step
#18. The user should preinstall the following R packages: tidyverse, here, rprojroot, ggplot2, optparse, RColorBrewer, viridis
#19. The user is advised to run the workflow using nohup to allow the script to run without an open terminal. Input the following three commands: nohup bash MOBELS_lcWGS_workflow_15iii23.sh; ctrl+x; bg. Terminal output will write to nohup.out. Monitor resource usage with top or htop commands

###########################################################################
#############################  Disclaimer  ################################
###########################################################################
#I cannot guarantee the workflow will work. If the workflow hangs up for some 
#reason, which it probably will, you will have to specify a checkpoint to restart 
#where things left off. This is important earlier in the workflow so as not to 
#reproduce costly computation during read trimming, mapping, and vcf compilation.
#The user is also responsible for understanding system capabilities and limitations.
#Depending on the amount of data being processed, disk space and RAM could be
#a limiting resource. These limitations vary widely and are project specific.

#Use the workflow mindfully.
```

Below are the user presets to modify in the bash script.
```
###########################################################################
#############################  Checkpoints ################################
###########################################################################
#The user should specify where to pick up the workflow here, whether after inspecting read or vcf quality or continuing after a hangup
trim_raw_reads=0 # specify yes [1] or no [0], yes will trim reads and stop the workflow to allow users to inspect read quality, no will pick up the workflow post trimming if user is happy with read quality
assemble_mitogenome=0 # specify yes [1] or no [0] to assemble using NOVOPlasty
assemble_plastid=0 # specify yes [1] or no [0] to assemble using NOVOPlasty
read_mapping=0 # specify yes [1] or no [0] to map reads or skip to vcf compilation, respectively 
vcf_compilation=1 # specify yes [1] or no [0], yes will compile vcf then stop the workflow to allow users to inspect raw vcf quality plots, which should inform filtering parameters
vcf_filtration=1 # specify yes [1] or no [0] to filter vcf or skip to running PCA and admixture, respectively

###########################################################################
############################  User presets  ###############################
###########################################################################
#Server presets, add executable pathways as needed. If command is available without pathway leave as is
threads=5 # number of available cpus for parallization of workflow tasks
java=java
trimmomatic=/media/genyoda/Extra_Storage/Projets/Data_Trevor/Trimmomatic-0.39/trimmomatic-0.39.jar
fastqc=fastqc
multiqc=multiqc
seqtk=seqtk
NOVOPlasty=NOVOPlasty4.3.1.pl # User specifies executable with version numbers here, which must appear in NOVOPlasty-master/NOVOPlastyX.X.X.pl from the job launch directory
bowtie2=/media/genyoda/Extra_Storage/Projets/Data_Trevor/bowtie2-2.4.5-linux-x86_64/bowtie2
samtools=samtools
bcftools=/media/genyoda/Extra_Storage/Projets/Data_Trevor/bcftools-1.16/bcftools
vcftools=vcftools
plink=/media/genyoda/Extra_Storage/Projets/Data_Trevor/plink_linux_x86_64_20230116/plink
admixture=admixture
R=R
Rscript=Rscript

#Project presets, including files
FILE=MOBELS_S_20_00703_Beluga_Chromosome3 # an informative prefix to carry through analyses
reference_nuclear_genome=S_20_00703_MOBELS_0202.22xii22.FINAL_scaffolds #Index name for reference genome, built using bowtie2. Should be <TextToExtractHere>.fasta
adapter=ill_adap_uni.fa # a fasta file with adapters to scan and remove during read trimming. A generic file is provided but can be swapped out by user
repeat_regions=beluga_00703_repeats.bed # repeat regions to exclude from VCF file
EXAMPLE_SAMPLE=2.AtlanticCanada.AT002 # an example sample name used in variant compilation step (uses this sample to extract contig names)
ploidy=2 # ploidy of the read datasets for the mpileup step; if the dataset is a mix of haploid and diploid individuals, you will need to specify this manually in the samples_list file used at the bcftools mpileup step

#forks in the workflow to specify
fastqc_rawreads=0 # specify yes [1] or no [0], this is to save a time consuming step but it may help to assess raw read QC to inform trimmomatic presets
delete_rawreads=0 # specify yes [1] or no [0], this will save space if raw reads are no longer needed following trimming
full_sites_vcf=1 # specify yes [1] or no [0] to compile an all-sites vcf that includes non-variant positions (important for distance based estimates such as nucleotide diversity and Fst estimates, which should be calculated using pixy)
compile_seqs_list=1 # specify yes [1] or no [0] if the user is providing the list seqs_list in launch directory of chromosomes/scaffolds/contigs to compile individually. Workflow will generate vcfs for each specified sequence and a concatenated dataset to continue working with
delete_sortedbams=0 # specify yes [1] or no [0], this will save space if sorted bams are no longer needed; only specify yes if user is confident dataset won't need to be recompiled and priority is disk space
filter_highcoverage=0 # specify yes [1] or no [0] if user has compiled at high coverage, adds extra filtering steps for coverage, allelic balance, genotype scores, and site scores (requires bcftools 1.13 or higher)
remove_repeat_sites=1 # specify yes [1] or no [0] if user has provided a bed file to exclude variant positions from repeat regions

#trimmomatic presets
trailing=20 # trim 3' end of reads until quality reaches specified value
headcrop=15 # remove specified number of basepairs from beginning of reads, typically needed as illumina base calling is not properly calibrated early in reads, but inspect raw read fastqc files
avgqual=25 # average quality of read required to retain
minlen=75 # minimum lenght of read to retain

#NOVOPlasty presets
Number_reads_to_extract=7500000 # organellar genomes are typically high copy, so a fraction of total read data will be needed. This will significantly reduced RAM usage. Try specifying more reads if genome does not assemble as a circular contig; plastid genomes will likely assemble into multiple fragments due to a common inverted repeat region
Genome_range_mito=30000-35000 # range in kbp
Genome_range_plastid=125000-140000 # range in kbp
Kmer=55 # try increasing this value if circular genome is not produced and coverage is high
Read_length=150 # tailor according to read data
Insert_size=200 # tailor according to read data
mito_seed=Alaria_marginata_cox1.fasta # config file set to not extend seed, so this can be a closely related species
plastid_seed=Alaria_marginata_rbcL.fasta # config file set to not extend seed, so this can be a closely related species

#Mapping parameter for bowtie2; 0.6=up to 10% divergence in mapping high quality reads in end-to-end mode, 0.3=up to 5%, 0.12=up to 2%; bowtie2 manual states minimum-score function f to f(x) = 0 + -0.6 * x, where x is read length.
map_param_nuc=0.6

#User inputs for filtering variant position artifacts with lcWGS datasets; see https://github.com/tbringloe/WGS-NOVAC for explanations on filtering options
minor_allele_frequency=0.004 # calculate based on dataset features and presets, including any signal unique to a sample(s) and missingness threshold (missing sites are not calculated into MAF by vcftools)
SNPgap=5 # filters genotypes within specified distance (bp) to indels, which are potentially enriched with errors (indels are removed downstream)
max_missing=0.8 # proportion of present data needed to keep a SNP site by vcftools, 0.9=10% missingness

#User inputs for filtering variant positions with high coverage WGS datasets; see https://github.com/tbringloe/WGS-NOVAC for explanations on filtering options; still need a filter for heterozygosity (i.e. duplications mapping to same non-duplicated region of reference genome) 
min_MQ=30 # minimum mapping quality score used at mpileup step, reads with a mapping quality score less than specified value are not considered, follows phred scale
min_BQ=30 # minimum base quality score used at mpileup step, bases with a quality score less than specified value are not considered, follows phred scale
min_cov=15 # minimum coverage to keep a variant position, should eliminate low confidence calls, though calling at low coverage appears to perform well. Could impact estimates of heterozygous sites
max_cov=100 # maximum coverage to keep a variant position, should eliminate high copy regions potentially enriched with mapping artifacts
allelic_balance_low=0.2 # ratio of reads matching reference and alternate call [0.2=1 reference to 5 alternate], the lower the ratio the more likely the position should be the alternate allele
allelic_balance_high=5 # ratio of reads matching reference and alternate call [5=5 reference to 1 alternate], the higher the ratio the more likely the position should be the reference allele
min_GQ=30 # filters individual genotypes according to their probability of being correct, follows a phred scale, 30=1/1000 chance of an allele calling error
min_Q=30 # filters sites according to the probability there is an alternate allele present, follows a phred scale, 30=1/1000 chance of SNP calling error

#User specifies following parameters for removal of linked SNPs using PLINK
plink_r2=0.25 # correlation coefficient above which to randomly remove sites within sliding window
plink_window_size=50 # sliding window size in kb

#User specifies values of k for ADMIXTURE analysis
k_low=2 # lowest number of ancestral populations to test
k_high=5 # highest number of ancestral populations to test
RBpalette=Set3 # RColourBrewer palette to use for admixture plots, see https://r-graph-gallery.com/38-rcolorbrewers-palettes; make sure there are enough colours for k_high
pops=CBS,FRB,JAM,NEH,NHB,NHS,NWH,REB,RES,SAN,SEH,SHS,SLE,SWH,UNG,NA #here, list populations in order they should appear in admixture plot, should be same terms used in naming convention
```

Once the user is happy with these choices, the user can run the script.

```
bash MOBELS_lcWGS_workflow_25v23_master.sh
```

The workflow, in terms of the methods leveraged, are described here. For further details on methodology and benchmarking of the specific programs used, see [References](#references), which include links to GitHub/webpages.
### Raw read QC
The workflow starts with optional fastqc/multiqc on the raw read files. The files are then trimmed using [trimmomatic](http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf) according to user specified criteria and fastqc and multiqc are performed on the trimmed reads for the user to evaluate. The workflow will exit until the user specifies they are satisfied with read quality in the user presets.

### Organellar assembly
The workflow has an optional step of assembling organellar genome(s) using [NOVOPlasty](https://github.com/ndierckx/NOVOPlasty) (Dierckxsens et al. 2017); this program extends a user provided seed sequence based on coverage (i.e. high frequency kmers are used to extend the seed), until the assembled genome begins to overlap, confirming circularity. The workflow first subsamples the amount of reads for assembly, and then assembles organellar genomes for all samples in parallel, greatly reducing linear time requirements. Carefully check over the final contig(s), if the assembly is not complete, this can potentially be improved two ways: 1) specify more reads to subsample for organellar assembly (improving coverage), or 2) specify a higher Kmer value, which may help bridge breaks in the assembly (you might have to increase the number of subreads as this will decrease kmer coverage). Performance will depend on the nature of the organellar genome. Complex samples with many different genomes at relatively comparable coverage, or repeat patterns, will lead to breaks in the assembly.

### Read mapping to nuclear genome and VCF compilation
The trimmed reads are mapped to the reference genome using end-to-end alignments in [bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml) (Langmead and Salzberg 2012). A single bowtie2 parameter is specified by the user, which indicates a cutoff in percent divergence to map to the reference genome. For phylogenetic work mapping across multiple species, this parameter should be set high (perhaps as high as 20% divergence or more, depending on how distantly related the species are expected to be relative to the reference genome), for population level work, this parameter could be set lower (2-5% divergence) if the user is concerned about non-target organisms mapping. The % divergence is approximate because the bowtie2 parameter being set is dynamic, since it factors in read quality information (poor base calls count less towards whether the read is mapped or not). Please refer to the [bowtie2 manual](https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml) to see exactly how scores are calculated.

The workflow will loop through mapping each sample. It will also convert the resulting SAM files to BAM using [samtools](http://www.htslib.org/doc/samtools.html) (Li et al. 2009), and sort the files for compilation using [bcftools](https://samtools.github.io/bcftools/bcftools.html) (Danacek et al. 2021). The workflow calls variant positions across samples on a contig by contig basis, before concatenating into one vcf, reducing linear time of the workflow.

The VCF file at this point represents the raw dataset. The workflow will generate several figures to help evaluate QC of the raw VCF file and exit. The user must use the output to inform parameters decisions on vcf filtration, then restart the workflow. The raw VCF is then filtered according to user specified criteria using a combination of [bcftools](https://samtools.github.io/bcftools/bcftools.html) (Danacek et al. 2021) and [vcftools](https://vcftools.sourceforge.net/man_latest.html) (Danacek et al. 2011). I have provided some guidance on how to set filtering parameters in the [Uncertainty](#uncertainty) section below.

### Genomic analyses to infer structure
The primary purpose of the workflow is to distill WGS datasets into a vcf file for analysis. Basic analyses to infer structure are provided, but users will undoubtedly have a wider array of analyses in mind. This includes the calculation of principal components through [plink](https://www.cog-genomics.org/plink/1.9/strat#pca), and admixture proportions using [ADMIXTURE](https://dalexander.github.io/admixture/admixture-manual.pdf). For guidance on how to (not over) interpret admixture proportions, see [Lawson et al. 2018.](https://www.nature.com/articles/s41467-018-05257-7).

For datasets compiled with invariant positions, [pixy](https://pixy.readthedocs.io/en/latest/) should be used for distance calculations (Korunes and Samuk 2021) in oder to correct for missing genotype bias.

## Requirements
**Data: paired-end short read data with the following formatting: OrderToAppearInAdmixturePlots.GeographicIdentifier(PCAlabels).SampleYear(optional).Species(optional).ExtractionID(optional).SampleID.MappedGenome(optional).R1.fastq.gz. Example: 7.RES.2007.BELUGA.ADN_20_05497.S_20_04096.dl.n_s.bam** it is critical the order of the samples follows the order you wish specimens to appear in the admixture plots, that is, grouped meaningfully by population or species. The naming format allows the workflow to derive species/population metadata, which is leveraged in some of the analyses. Alternatively, the naming convention is less critical if leveraging alternative MOBELS_plots.R code provided here. 

**A reference genome** stored as a fasta file.

**Mitochondrial and plastid seed fasta files** (i.e. a gene in the target, or a closely related organism) stored in NOVOPlasty/seed_files

**Optional BED file containing repeat regions of the reference genome**, to be filtered from the VCF files, stored in reference_genome/.

**Program dependencies:** fastqc v0.11.9, multiQC v1.9, trimmomatic v0.39, NOVOPlasty v4.2, seqtk v1.3, perl v.5.34.0, bowtie2 v2.4.2, samtools v1.13, bcftools v1.12, vcftools v0.1.16, r v4.0.4, admixture v1.3.0, plink v2.00 

**R packages:** tidyverse, here, rprojroot, ggplot2, optparse, RColorBrewer, viridis

## Caveats
Nucleotide diversity and genetic distances in the nuclear datasets are severely biased because missing data is improperly treated as the reference allele by vcftools. VCF files with invariant position data should be used to calculate nucleotide diversity and genetic distances using [pixy](https://pixy.readthedocs.io/en/latest/) (Korunes and Samuk 2021). Pixy is not directly integrated into the workflow, so users will have to run that program independently.

Standard practices and justifications for filtering parameters related to nuclear wide variant calling are not obvious and likely highly dependent on your chosen study system. Because these decisions are individualized across datasets, users should cautiously interpret values across analyses (e.g. distance based measurements such as Fst, inbreeding coefficients, ect). Results should always be interpreted within the context of the data analysed and decisions made during that process.

The quality and reliability of the workflow inherently depends to some degree on the quality of the reference genome. The user should ensure the reference genome is as clean as possible from non-target taxa. I suggest a [blobtools](https://github.com/DRL/blobtools) (https://github.com/DRL/blobtools) approach to assessing the purity of the reference genome, and to use coverage, GC content, and taxonomic information to filter non-target contigs. Generating a solid reference genome is a logistical challenge in-and-of itself. For a suggested workflow generating a reference assembly, see XXXlink to assembly github pageXXX.

The workflow assumes a consistent ploidy level across all individuals. If this is not the case (i.e. datasets are a mixture of haploid and diploid individuals), this must be specified in the sample_list file used by bcftools at the mpileup step of the workflow. This cannot be automated, so the user will have to add ploidy manually to the file as a second tab-delimited column and modify the script accordingly.

The workflow is missing two potentially important filtering steps users should consider incorporating independently. Related individuals shouold be removed; this can be calculated using --relatedness and --relatedness2 flags in [vcftools](https://vcftools.sourceforge.net/man_latest.html). The workflow also does not account for high site heterozygosity, which could indicate reads are mapped to paralogous regions, resulting in incorrect heterozygous calls. Several steps already implemented should help mitigate this: the removal of repetitive regions, and setting a high mapping quality score (30+) should ensure sites are called without reads mapping to multiple locations of the genome. In the case of multiple regions of a sample mapping to a single region in the reference genome, this is better addressed by removing high coverage regions, or highly heterozygous sites, which can be calculated using the [populations](https://catchenlab.life.illinois.edu/stacks/comp/populations.php) module in stacks, and subsequently removed using the --exclude-positions flag in [vcftools](https://vcftools.sourceforge.net/man_latest.html). A polite reminder is provided in the workflow.
```
echo use populations module of stacks to calculate heterozygosity per site defining all individuals are single population
  echo populations -V <in>.vcf.gz -M <population_structure>.txt -t <threads> --out-path .
  echo decide on a threshold for acceptable site heterozygosity and remove sites above threshold
  echo vcftools --gzvcf <in>.vcf.gz --exclude-positions <heterozygous_sites>.txt --recode --out <out>
```
The workflow has been tested using the exact versions in the above dependencies. No effort has been made to test other versions, or enhance compatibility of dependencies.

Finally, my impression is that main patterns are highly robust to the parameter decisions made by the user. The user should nontheless take the time to convince themselves of this.

## Uncertainty

The workflow is meant to automate and dampen learning curves as much as possible for the user, but important decisions must be made regarding filtering criteria for genome-wide nuclear variants. Some guidance from personal experience are provided here for the user to consider, but it is no substitute for exploring individual datasets and making informed decisions regarding these parameters. Note that most of the filtering parameters below apply to high coverage WGS datasets.

*User inputs following filtering settings for compiling VCF file - Presets and justification*

**min_cov=15** # 15 is a good target for minimum read coverage for calling a SNP. There is a certain percentage chance the base call for a read will be incorrect, and this chance of an incorrect call, when upscaled to the entire genome, can results 1000s of incorrect SNP calls. Stacking reads reduces the chance of calling an incorrect SNP. That being said, I've had good success calling SNPs in large mammalian genomes using low coverage (~4-6x), results that were validated with higher coverage ddRADseq datasets. It is for this reason workflow pathways were incorporated to facilitate lcWGS datasets. Note that the current workflow uses variant calling positions, but approaches integrating genotype-likelihoods are gaining traction. See [Lou et al. 2021](https://onlinelibrary.wiley.com/doi/10.1111/mec.16077).

**max_cov=100** # Users should evaluate coverage across samples (site coverage and individual coverage) to inform an upper mean threshold present in the read dataset. If a handful of samples deviate a lot from the mean towards higher coverage values, you can always subsample reads from those samples. By getting a reasonably even distribution of coverage across samples, you can set this threshold to strategically eliminate high copy regions of the genome (which could lead to mapping errors and spurious variant calls). Setting an upper coverage threshold also removes sites corresponding to the abundant organellar genomes (if these have not already been removed from the reference genome provided).

**allelic_balance_low=0.2** # Heterozygous sites will be called if mapping reads from diploid specimens to a nuclear genome. The confidence of calling a heterozygous site depends, in part, on how balanced the reads are between the two alleles. The expectation is that the reads should be evenly divided between the two alleles. Anything highly skewed from this 50:50 distribution could indicate contamination of the sample (low number of reads supporting an alternate allele). The user must make a decision, however informed, of where that threshold lies for when reads are balanced enough to call a heterozygous site. The 1:5 ratio is somewhat arbitraty, attempting to balance the removal of artifacts while not losing real signal. Generous read coverage is especially important for this filtering; allelic balance will be prone to more randomness at lower coverage, making it harder to retain sites with confidence.

**allelic_balance_high=5** # Same justification as above. This simply refers to skewed ratios favouring the alternative heterozygous allele.

**minor_allele_frequency=0.05** # Two reasons exist to filter low frequency alleles from SNP datasets. One, they may represent artifacts; some level of replication of the genotype across individuals should be a prerequisite before keeping a given genotype. Second, low frequency alleles provide little information in the way of species or population structure (i.e. contributes noise), and can therefore be removed to unobscure genomic signal. The threshold will depend on the scope of the dataset and user needs. If the user is interested in keeping genotypes specific to only a few individuals (e.g. maybe a low sample size for a population or species), the threshold will have to be lowered accordingly to allow those SNPs into the dataset. The user should also remember all new mutations start at a low frequency. The absence of low frequency alleles can be indicative of population history, such as a recent bottleneck. When calculating nucleotide diversity or Tajimas D, the user mayshould consider skipping this filter.

**min_MQ, min_BQ, min_GQ, and min_Q=30** # These scores are a measure of the expected probability of a correct mapping/base/genotype/variant site call. This score follows a phred scale, such that scores of 30=1/1000 chance of an error. MQ and BQ operate at the read level, and are applied when calling variant positions. The MQ parameters refers to mapping scores; a score of 30 would offer a high degree of confidence the read mapped to a single location of the genome (1/1000 chance the read does not uniquely map). BQ refers to base call, such that only bases with less than a 1/1000 chance of being incorrect are considered when calling variant positions. The GQ and Q parameters are applied after variant positions are called. The GQ parameter works at the individual genotype level (GQ format field in vcf; probability the variant call for a given specimen at a given position is an error), whereas Q works at the site level (QUAL column in vcf; probability of an alternate allele present at a given site). Altogether, if all Q parameters are set to 30, a user could state that variant positions were called with reads with a less than 1/1000 chance of a mapping or base call error, and that variant genotypes/positions have a less than 1/1000 chance of a calling error. 

**SNPgap** # Indels are expected to be enriched with genotyping errors, and therefore the user may want to filter genotypes adjacent to indels. Five seems like a reasonably conservative value, but if not informed in any way.

**max_missing=0.9** # This value represents the proportion of non-missing data needed to keep a SNP site, where 0.9=10% missingness. The amount of missing data the user is willing to tolerate depends on the type of analyses being conducted. When working with population genomics and a reference genome for the target species, the amount of missing gemotypes is typically quite low (<10%). Missing data is expected to increase when using a reference genome from a species other than the target. If conducting phylogenetic analyses, setting a high threshold (e.g. 0% missingness) would help to ensure conserved homologous sites are analysed.

**plink_r2=0.15 and plink_window_size=25 # in kb** # Filtering for linkage disequilibrium is a key assumption when evaluating population structure or phylogenetic signal. Linked positions might lead to an uneven distribution of phylogenetic signal across the genome, with accelerated evolution in some areas of the genome leading to "spikes" in signal related to non-random evolutionary events; pruning for linked sites is a way of normalizing genomic signal, and not interpreting non-random evolutionary events (such as selective forces) for random ones (like population structure that accumulates over time through random mutations and genetic drift). Setting an r2 value (r=correlation coefficient) is not straightforward, and neither is the window size. The best approach is to plot the decay in linkage disequilibiurm independently using [PopLDdecay](https://github.com/BGI-shenzhen/PopLDdecay) (https://github.com/BGI-shenzhen/PopLDdecay), and use this to inform the choices made here.

**Heterozygosity=0.6** #  There is no clear cutoff for when to remove sites with higher than expected levels of heterozygosity. As with allelic balance, the idea here is to find a value that removes more noise than signal, which will depend in part on the knowledge and intuition of the user. This filter must be applied independently of the workflow presentd here.

**Relatedness=0.25** # Setting a kinship coefficient threshold to 0.25 ensures siblings and and parent offspring relationships are removed from the dataset. Setting a threshold depends on end-user needs, and whther related individuals are likely to impact interpretation of patterns. This filter must be applied independently of the workflow presentd here.

## Acknowledgements
All the inspirational students, postdocs, mentors, DFO research scientists, and forum junkies across the globe who contributed to my own learning journey in bioinformatics, of which this workflow is a direct result.

## References
**Please see useful pages here:**

Scripts provided by joanam: https://github.com/speciationgenomics/scripts

Very useful guide for filtering and visualizing SNPs (from which some of this workflow is derived): https://speciationgenomics.github.io/filtering_vcfs/ 

Mermaid wofkflow diagram: https://mermaid.live/

Clever one liner for VCF compilation step: https://www.ecseq.com/support/ngs-snippets/how-to-run-time-consuming-data-analysis-processes-in-parallel-on-unix-systems

Deep dive on bowtie2's mapping quality score: https://biofinysics.blogspot.com/2014/05/how-does-bowtie2-assign-mapq-scores.html

**Citations for programs:**

Alexander, D. H., Lange, K. (2011). Enhancements to the ADMIXTURE algorithm for individual ancestry estimation. BMC Bioinformatics, 12, 1-6. https://dalexander.github.io/admixture/download.html

Danacek, P., Auton, A., Goncalo, A., Albers, C. A., Banks, E., DePristo, M. A., Handsaker, R., Lunter, G., Marth, G., Sherry, S. T., McVean, G., Durbin, R. & 1000 Genomes Project Analysis Group. (2011). The variant call format and VCFtools. Bioinformatics 27: 2156-8. http://vcftools.sourceforge.net/

Danacek, P., Bonfield, J. K., Liddle, J., Marshall, J., Ohan, V., Pollard, M. O., Whitwham. A., Keane, T., McCarthy, S. A., Davies, R. M., & Li, H. (2021). Twelve years of SAMtools and BCFtools. Gigascience, 10, giab008. https://github.com/samtools/bcftools.

Dierckxsens, N., Mardulyn, P., & Smits, G. (2017). NOVOPlasty: de novo assembly of organelle genomes from whole genome data. Nucleic Acids Research, 45, e18. https://github.com/ndierckx/NOVOPlasty

Huson, D. H., & Bryant, D. (2006). Application of phylogenetic networks in evolutionary studies. Molecular Biology and Evolution, 23, 254-67. https://software-ab.informatik.uni-tuebingen.de/download/splitstree4/welcome.html

Korunes, K. L., Samuk, K. (2021). Pixy: unbiased estimation of nucleotide diversity and divergence in the precense of missing data. Molecular Ecology Resources. 21: 1359-1368.

Langmead, B., & Salzberg, S. L. (2012). Fast gapped-read alignment with Bowtie 2. Nature Methods, 9, 357-9. http://bowtie-bio.sourceforge.net/bowtie2/index.shtml

Lawson, D, J., van Dorp, L., Falush, D. 2018. A tutorial on how to not over-interpret STRUCTURE and ADMIXTURE bar plots. Nature Communications, 9, 3258.

Li, H., Handsaker, B., Wysoker, A., Fennell, T., Ruan, J., Homer, N., Marth, G., Abecasis, G., Durbin, R., & 1000 Genome Project Data Processing Subgroup. (2009). The sequence alignment/map format and SAMtools. Bioinformatics 25: 2078-9. http://www.htslib.org/

Lou, R. N., Jacobs, A., Wilder, A. P., Therkildsen, N. O. 2021 A beginner's guide to low-coverage whole genome sequencing for population genomics. Molecular Ecology, 30, 5966-5993.

Purcell, S., Neale, B., Todd-Brown, K., Thomas, L., Ferreira, M. A. R., Bender, D., Maller, J., Sklar, P., de Bakker, P. I. W., Daly, M. J., & Sham, P. C. (2007). PLINK: a tool set for whole-genome association and population based linkage analyses. American Journal of Human Genetics, 81, 559-575. https://www.cog-genomics.org/plink.
