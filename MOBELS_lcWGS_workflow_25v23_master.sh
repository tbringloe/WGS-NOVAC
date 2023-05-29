#!/bin/bash
###############Contact information####################
#Trevor T. Bringloe, PhD
#Department of Fisheries and Ocean, Research Scientist
#Direction des sciences démersales et benthiques
#850 route de la Mer
#Mont-Joli (QC)
#G5H 3Z4
#506-259-2288

##This is a workflow to map whole genome reads to a reference genome at low coverage, call SNPs, and perform basic population analyses. For workflow with more filtering criteria (i.e. mapping at high coverage), see https://github.com/tbringloe/WGS-NOVAC
##Originally designed and used on kelp, adapted for beluga and narwhal, and whatever organism is your passion##

#,,,,,,,,,,,,,,,,,,,,,,:,::::,:,,,,,,,,,,,,,,,,,,,,,,,,\sss,,,,,,,,,,sss/,,
#,,,,,,,,,,,,,,,,,,,,,,,:,::,:,,,,,,,,,,,,,,,,,,,,,,,,,,\w888-\,,,/8888/,,,
#,,,,,,,,,,,,,,,,,,,,,,,,::::,,,,,,,,,,,,,,,,,,,,,,,,,,,,\8888j|8/8888/,,,,
#,,,,,,,,,,,,,,,,,,,,,,,__::__,,,,,,,,,,,,,,,,,,,,,,,,,,,,\8888888888/,,,,,
#,,,,,,,,,,,,,,,,,/88888888888888oooooo.....,,,,,,,,,,,,,,,,'j8888j',,,,,,,
#,,,,,,,,,,,,,,,,,/8888888888888888888888888........,,,,,,,,,,]88[,,,,,,,,,
#,,,,,,,,,,,,,,,,,|88888888888888888888888888888888888888888888888],,,,,,,,
#,,,,,,,,,,,,,,,,,|88888888888-~-888888888888888888888888888888888/,,,,,,,,
#,,,,,,,,,,,,,,,,,|88888888:=[[O]]=:88888888888888888888888888888/,,,,,,,,,
#<<<<<<<<<<<<<<<<</888888888sssssss8888888888888888888888888888D',,,,,,,,,,
#,,,,,,,,,,,,,,,,,'''...DFO...>:888888888888888888888888888888D,,,,,,,,,,,,
#,,,,,,,,,,,,,,,,,,,,,,,,.....nn/8888888888888888\\888888aa'''',,,,,,,,,,,,
#,,,,,,,,,,,,,,,,,,,,,,,......nn''88o'''''''''''''aaao'''''''',,,,,,,,,,,,,
#,,,,,,,,,,,,,,,,,,,,,,,.....n88o'',,,,,,,,,,,,,,,''aao,,,,,,,,,,,,,,,,,,,,
#,,,,,,,,,,,,,,,,,,,,,,,,,,[jj',,,,,,,,,,,,,,,,,,,,,,'ao,,,,,,,,,,,,,,,,,,,
#,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,

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

###########################################################################
##############################  Workflow  #################################
###########################################################################

#####Start workflow from beginning, or skip to read mapping if reads are already trimmed
if [ $trim_raw_reads  -eq  1 ]
then
  ######Make directories for worflow
  echo making directories
  mkdir bowtie2_output
  mkdir IDX
  mkdir trimmed_reads
  mkdir results
  mkdir results/read_QC
  mkdir results/read_QC/raw
  mkdir results/read_QC/trimmed
  mkdir results/bowtie2_output
  mkdir results/vcf_QC
  mkdir results/vcf_QC/raw
  mkdir results/vcf_QC/filtered
  mkdir results/pca
  mkdir results/admixture
  mkdir vcf_analyses
  mkdir vcf_analyses/vcf_stats
  mkdir vcf_analyses/vcf_stats/raw/
  mkdir vcf_analyses/vcf_stats/filtered/
  mkdir vcf_analyses/plink
  mkdir vcf_analyses/admixture
  mkdir NOVOPlasty-master/read_files
  mkdir Mitogenomes
  mkdir Plastomes
  mkdir tmp
  mkdir vcf
  mkdir sorted_bams

  #create sample list of raw reads for processing
  cd raw_reads
  ls *_R1.fastq.gz > sample.list
  sed 's/_R1.fastq.gz//g' sample.list > sample.list2
  rm sample.list
  mv sample.list2 sample.list
  mv sample.list ..
  cd ..

  ######Quality check and trim raw reads

  #Optional check for raw read quality
  if [ $fastqc_rawreads  -eq  1 ]
  then
    echo checking raw read quality
    cat sample.list | xargs -I {} -n 1 -P $threads sh -c "$fastqc raw_reads/{}_R1.fastq.gz raw_reads/{}_R2.fastq.gz -o results/read_QC/raw"
    multiqc results/read_QC/raw -o results/read_QC
    mv results/read_QC/multiqc_report.html results/read_QC/raw_reads_multiqc_report.html
  else
    echo not checking QC of raw read files
  fi

  #Loop trimmomatic over raw reads
  echo trimming reads
  cat sample.list | while read line
  do 
  $java -jar $trimmomatic PE -threads 32 raw_reads/"$line"_R1.fastq.gz raw_reads/"$line"_R2.fastq.gz trimmed_reads/"$line"_R1_trimmed.fastq.gz trimmed_reads/"$line"_R1_tup.fastq.gz trimmed_reads/"$line"_R2_trimmed.fastq.gz trimmed_reads/"$line"_R2_tup.fastq.gz ILLUMINACLIP:raw_reads/$adapter:2:30:10 TRAILING:$trailing HEADCROP:$headcrop AVGQUAL:$avgqual MINLEN:$minlen

  #optional deletion of raw read files
  if [ $delete_rawreads  -eq  1 ]
  then
    rm "$line"_R1.fastq.gz
    rm "$line"_R2.fastq.gz
  else
    echo not deleting raw read files
  fi

  rm trimmed_reads/"$line"_R1_tup.fastq.gz
  rm trimmed_reads/"$line"_R2_tup.fastq.gz
  done

  #Check for trimmed read quality
  echo checking trimmed read quality
  cat sample.list | xargs -I {} -n 1 -P $threads sh -c "$fastqc trimmed_reads/{}_R1_trimmed.fastq.gz trimmed_reads/{}_R2_trimmed.fastq.gz -o results/read_QC/trimmed"
  $multiqc results/read_QC/trimmed -o results/read_QC
  mv results/read_QC/multiqc_report.html results/read_QC/trimmed_reads_multiqc_report.html
  echo stopping workflow to allow user to check read quality
  exit 0
else
  echo User is satisfied with read quality, moving onto organellar assembly and/or read mapping
fi

if [ $assemble_mitogenome  -eq  1 ]
then
  #####Assemble mitogenomes
  echo assembling mitogenomes
  #Extract reads for downstream assembly
  echo extracting reads for downstream processing
  cat sample.list | xargs -I {} -n 1 -P $threads sh -c "$seqtk sample -s100 trimmed_reads/{}_R1_trimmed.fastq.gz $Number_reads_to_extract > NOVOPlasty-master/read_files/{}_sub_R1.fastq"
  cat sample.list | xargs -I {} -n 1 -P $threads sh -c "$seqtk sample -s100 trimmed_reads/{}_R2_trimmed.fastq.gz $Number_reads_to_extract > NOVOPlasty-master/read_files/{}_sub_R2.fastq"

  #setup NOVOPlasty environment to run multiple samples in parallel, including specifying parameters in config file
  echo setting up NOVOPlasty environments
  cat sample.list | while read line
  do
  mkdir NOVOPlasty-master/$line
  cp NOVOPlasty-master/lcWGS_config_mito.txt NOVOPlasty-master/lcWGS_config_mito_"$line".txt
  sed -i "s/PROJECT_NAME/$line/g" NOVOPlasty-master/lcWGS_config_mito_"$line".txt
  sed -i "s/GENOME_RANGE_MITO/$Genome_range_mito/g" NOVOPlasty-master/lcWGS_config_mito_"$line".txt
  sed -i "s/KMER/$Kmer/g" NOVOPlasty-master/lcWGS_config_mito_"$line".txt
  sed -i "s/MITO_SEED/$mito_seed/g" NOVOPlasty-master/lcWGS_config_mito_"$line".txt
  sed -i "s/READ_LENGTH/$Read_length/g" NOVOPlasty-master/lcWGS_config_mito_"$line".txt
  sed -i "s/INSERT_SIZE/$Insert_size/g" NOVOPlasty-master/lcWGS_config_mito_"$line".txt
  sed -i "s/FORWARD/"$line"_sub_R1.fastq/g" NOVOPlasty-master/lcWGS_config_mito_"$line".txt
  sed -i "s/REVERSE/"$line"_sub_R2.fastq/g" NOVOPlasty-master/lcWGS_config_mito_"$line".txt
  mv NOVOPlasty-master/lcWGS_config_mito_"$line".txt NOVOPlasty-master/$line
  cp NOVOPlasty-master/filter_reads.pl NOVOPlasty-master/$line
  cp -a NOVOPlasty-master/"$NOVOPlasty" NOVOPlasty-master/$line
  done
  #Run NOVOPlasty in parallel for multiple samples
  echo running NOVOPlasty in parallel for all samples
  cat sample.list | xargs -I {} -n 1 -P $threads sh -c "NOVOPlasty-master/{}/$NOVOPlasty -c NOVOPlasty-master/{}/lcWGS_config_mito_{}.txt"
  #Output is in working directory, so move to relevant sample folder
  echo Moving NOVOplasty output
  cat sample.list | while read line
  do
  sed -i "s/Contig1/$line/g" Circularized_assembly_*_"$line"_mito.fasta
  mv log_"$line"_mito.txt NOVOPlasty-master/"$line"
  mv contigs_tmp_"$line"_mito.txt NOVOPlasty-master/"$line"
  mv Circularized_assembly_*_"$line"_mito.fasta Mitogenomes
  rm NOVOPlasty-master/read_files/"$line"_sub_R*.fastq
  done
  echo if sample mitogenome does not appear in Mitogenome folder check relevant sample folder in NOVOPlasty-master/
else
  echo skipping mitogenome assembly
fi
if [ $assemble_plastid  -eq  1 ]
then
  ####Assemble plastid genomes
  echo assembling plastomes
  #Extract reads for downstream assembly
  echo extracting reads for downstream processing
  cat sample.list | xargs -I {} -n 1 -P $threads sh -c "$seqtk sample -s100 trimmed_reads/{}_R1_trimmed.fastq.gz $Number_reads_to_extract > NOVOPlasty-master/read_files/{}_sub_R1.fastq"
  cat sample.list | xargs -I {} -n 1 -P $threads sh -c "$seqtk sample -s100 trimmed_reads/{}_R2_trimmed.fastq.gz $Number_reads_to_extract > NOVOPlasty-master/read_files/{}_sub_R2.fastq"

  #setup NOVOPlasty environment to run multiple samples in parallel, including specifying parameters in config file
  echo setting up NOVOPlasty environments
  cat sample.list | while read line
  do
  mkdir NOVOPlasty-master/$line
  cp NOVOPlasty-master/lcWGS_config_plastid.txt NOVOPlasty-master/lcWGS_config_plastid_"$line".txt
  sed -i "s/PROJECT_NAME/$line/g" NOVOPlasty-master/lcWGS_config_plastid_"$line".txt
  sed -i "s/GENOME_RANGE_PLASTID/$Genome_range_plastid/g" NOVOPlasty-master/lcWGS_config_plastid_"$line".txt
  sed -i "s/KMER/$Kmer/g" NOVOPlasty-master/lcWGS_config_plastid_"$line".txt
  sed -i "s/PLASTID_SEED/$plastid_seed/g" NOVOPlasty-master/lcWGS_config_plastid_"$line".txt
  sed -i "s/READ_LENGTH/$Read_length/g" NOVOPlasty-master/lcWGS_config_plastid_"$line".txt
  sed -i "s/INSERT_SIZE/$Insert_size/g" NOVOPlasty-master/lcWGS_config_plastid_"$line".txt
  sed -i "s/FORWARD/"$line"_sub_R1.fastq/g" NOVOPlasty-master/lcWGS_config_plastid_"$line".txt
  sed -i "s/REVERSE/"$line"_sub_R2.fastq/g" NOVOPlasty-master/lcWGS_config_plastid_"$line".txt
  mv NOVOPlasty-master/lcWGS_config_plastid_"$line".txt NOVOPlasty-master/$line
  cp NOVOPlasty-master/filter_reads.pl NOVOPlasty-master/$line
  cp -a NOVOPlasty-master/"$NOVOPlasty" NOVOPlasty-master/$line
  done
  #Run NOVOPlasty in parallel for multiple samples
  echo running NOVOPlasty in parallel for all samples
  cat sample.list | xargs -I {} -n 1 -P $threads sh -c "NOVOPlasty-master/{}/$NOVOPlasty -c NOVOPlasty-master/{}/lcWGS_config_plastid_{}.txt"
  #Output is in working directory, so move to relevant sample folder
  echo Moving NOVOplasty output
  cat sample.list | while read line
  do
  sed -i "s/Contig1/$line/g" Circularized_assembly_*_"$line"_plastid.fasta
  mv log_"$line"_plastid.txt NOVOPlasty-master/"$line"
  mv contigs_tmp_"$line"_plastid.txt NOVOPlasty-master/"$line"
  mv Circularized_assembly_*_"$line"_plastid.fasta Plastomes
  rm NOVOPlasty-master/read_files/"$line"_sub_R*.fastq
  done
  echo if sample plastome does not appear in Plastome folder check relevant sample folder in NOVOPlasty-master/
else
  echo skipping plastome assembly
fi

#####Map reads to reference genome
if [ $read_mapping  -eq  1 ]
then
  echo starting read mapping to call variant positions
  #Build bowtie2 index for mapping providing user is satisfied with read quality checks
  echo building bowtie2 index for read mapping
  "$bowtie2"-build reference_genome/"$reference_nuclear_genome".fasta IDX/"$reference_nuclear_genome"

  #Map reads to reference genome using bowtie2, pipe into samtools to convert output to sorted bam, then index
  echo mapping reads to reference genome and converting to sorted bam format
  cat sample.list | while read line
  do
  $bowtie2 --score-min L,0,-"$map_param_nuc" --no-unal -p $threads -x IDX/"$reference_nuclear_genome" -1 trimmed_reads/"$line"_R1_trimmed.fastq.gz -2 trimmed_reads/"$line"_R2_trimmed.fastq.gz -S "$line"_"$reference_nuclear_genome".sam 2> results/bowtie2_output/"$line"_"$reference_nuclear_genome".sd.out
  $samtools view -S -@ "$threads" -b "$line"_"$reference_nuclear_genome".sam | \
  $samtools sort -@ "$threads" -o "$line"_"$reference_nuclear_genome".n_s.bam -
  $samtools index -@ "$threads" "$line"_"$reference_nuclear_genome".n_s.bam
  rm "$line"_"$reference_nuclear_genome".sam
  done

  #Generate multiqc report of bowtie2 output
  echo generating multiqc report for read mapping
  #$multiqc results/bowtie2_output -o results/bowtie2_output
  #mv results/bowtie2_output/multiqc_report.html results/bowtie2_output/mapped_reads_multiqc_report.html
  echo inspect multiqc report for mapped reads in results/bowtie2_output
else
  echo skipping read mapping
fi

#####Compile vcf
if [ $vcf_compilation  -eq  1 ]
then
  ls -v *.n_s.bam > sample_list1
  ls -v *.n_s.bam > sample_list2
  sed -i "s/$/\    $ploidy/" sample_list2
  #index the reference genome for mpileup
  $samtools faidx reference_genome/"$reference_nuclear_genome".fasta
  if [ $compile_seqs_list  -eq  1 ]
  then
    if [ $full_sites_vcf  -eq  1 ]
    then
      echo compiling full site vcfs
      cat seqs_list | xargs -I {} -n 1 -P $threads sh -c "$bcftools mpileup --threads 1 -Ou -f reference_genome/$reference_nuclear_genome.fasta -q "$min_MQ" -Q "$min_BQ" -a AD,ADF,ADR,DP -r {} -b sample_list1 | $bcftools call --threads 1 -S sample_list2 -f GQ,GP -m -Oz > vcf/tmp.{}.vcf.gz"
      #list temporary vcf files
      ls vcf/tmp* > list_tmp
      concat into single vcf
      $bcftools concat --threads $threads -f list_tmp -Oz -o vcf/"$FILE"_raw.vcf.gz
      rm list_tmp
      mv *n_s.bam sorted_bams
      mv *n_s.bam.bai sorted_bams
    else
      echo compiling variant only vcfs
      cat seqs_list | xargs -I {} -n 1 -P $threads sh -c "$bcftools mpileup --threads 1 -Ou -f reference_genome/$reference_nuclear_genome.fasta -q "$min_MQ" -Q "$min_BQ" -a AD,ADF,ADR,DP -r {} -b sample_list1 | $bcftools call --threads 1 --variants-only -S sample_list2 -f GQ,GP -mv -Oz > vcf/tmp.{}.vcf.gz"
      #list temporary vcf files
      ls vcf/tmp* > list_tmp
      #concat into single vcf
      $bcftools concat --threads $threads -f list_tmp -Oz -o vcf/"$FILE"_raw.vcf.gz
      rm list_tmp
      mv *n_s.bam sorted_bams
      mv *n_s.bam.bai sorted_bams
    fi
  else
    if [ $full_sites_vcf  -eq  1 ]
    then
      #mpileup and call on each contig, specified parallel pileups by creating xargs list of each contig. Use samtools view to obtain a list of all contig names, and input each as a xargs command, with -P as number of parallel jobs
      echo compiling full site vcf
      $samtools view -@ "$threads" -H "$EXAMPLE_SAMPLE"_"$reference_nuclear_genome".n_s.bam | grep "\@SQ" | sed 's/^.*SN://g' | cut -f 1 | xargs -I {} -n 1 -P $threads sh -c "$bcftools mpileup --threads 1 -Ou -f reference_genome/$reference_nuclear_genome.fasta -q "$min_MQ" -Q "$min_BQ" -a AD,ADF,ADR,DP -r {} -b sample_list1 | $bcftools call --threads 1 -S sample_list2 -f GQ,GP -m -Oz > tmp/tmp.{}.vcf.gz"
      #list temporary vcf files
      ls tmp/tmp* > list_tmp
      #concat into single vcf
      $bcftools concat --threads $threads -f list_tmp -Oz -o vcf/"$FILE"_raw.vcf.gz
      rm tmp/tmp*
      rm list_tmp
      mv *n_s.bam sorted_bams
      mv *n_s.bam.bai sorted_bams
    else
      #mpileup and call on each contig, specified parallel pileups by creating xargs list of each contig. Use samtools view to obtain a list of all contig names, and input each as a xargs command, with -P as number of parallel jobs
      echo compiling variant only vcf
      $samtools view -@ "$threads" -H "$EXAMPLE_SAMPLE"_"$reference_nuclear_genome".n_s.bam | grep "\@SQ" | sed 's/^.*SN://g' | cut -f 1 | xargs -I {} -n 1 -P $threads sh -c "$bcftools mpileup --threads 1 -Ou -f reference_genome/$reference_nuclear_genome.fasta -q "$min_MQ" -Q "$min_BQ" -a AD,ADF,ADR,DP -r {} -b sample_list1 | $bcftools call --threads 1 --variants-only -S sample_list2 -f GQ,GP -mv -Oz > tmp/tmp.{}.vcf.gz"
      #list temporary vcf files
      ls tmp/tmp* > list_tmp
      #concat into single vcf
      $bcftools concat --threads $threads -f list_tmp -Oz -o vcf/"$FILE"_raw.vcf.gz
      rm tmp/tmp*
      rm list_tmp
      mv *n_s.bam sorted_bams
      mv *n_s.bam.bai sorted_bams
    fi
  fi

  #optional deletion of sorted bam files
  if [ $delete_sortedbams  -eq  1 ]
  then
    rm -R sorted_bams
  else
    echo not deleting sorted bams
  fi

  #Check raw vcf quality to inform filtering parameters
  ##QC plots of the raw vcf dataset; aspects of the plots can and should be altered in the R script; please restart pipeline from here if initial thresholds can be optimized
  echo generating raw vcf quality plots
  $vcftools --gzvcf vcf/"$FILE"_raw.vcf.gz --depth --out vcf_analyses/vcf_stats/raw/"$FILE"_raw
  $vcftools --gzvcf vcf/"$FILE"_raw.vcf.gz --site-mean-depth --out vcf_analyses/vcf_stats/raw/"$FILE"_raw
  $vcftools --gzvcf vcf/"$FILE"_raw.vcf.gz --freq2 --max-alleles 2 --out vcf_analyses/vcf_stats/raw/"$FILE"_raw
  $vcftools --gzvcf vcf/"$FILE"_raw.vcf.gz --site-quality --out vcf_analyses/vcf_stats/raw/"$FILE"_raw
  $vcftools --gzvcf vcf/"$FILE"_raw.vcf.gz --missing-indv --out vcf_analyses/vcf_stats/raw/"$FILE"_raw
  $vcftools --gzvcf vcf/"$FILE"_raw.vcf.gz --missing-site --out vcf_analyses/vcf_stats/raw/"$FILE"_raw
  $vcftools --gzvcf vcf/"$FILE"_raw.vcf.gz --het --out vcf_analyses/vcf_stats/raw/"$FILE"_raw
  cp R_code/01_vcfQC_RAW.R R_code/01_vcfQC_"$FILE"_RAW.R
  sed -i "s/FILE/$FILE/g" R_code/01_vcfQC_"$FILE"_RAW.R
  $R CMD BATCH R_code/01_vcfQC_"$FILE"_RAW.R
  rm R_code/01_vcfQC_"$FILE"_RAW.R
  echo user should inspect raw vcf plots in results/vcf_QC/raw and modify filtering presets accodingly before restarting the workflow
  exit 0
else
  echo skipping vcf compilation
  echo user is satisfied with filtering presets, moving onto filtering vcf
fi

#####Filter vcf
if [ $vcf_filtration  -eq  1 ]
then
  if [ $filter_highcoverage  -eq  1 ]
  then
    if [ $remove_repeat_sites  -eq  1 ]
    then
      #removing sites from repetitive regions
      echo filtering high coverage vcf
      echo removing sites from repetitive regions
      $vcftools --gzvcf vcf/"$FILE"_raw.vcf.gz --exclude-bed reference_genome/$repeat_regions --recode --out vcf/"$FILE"_no_repeats
      echo setting heterozygous calls biased towards reference to reference allele
      $bcftools +setGT vcf/"$FILE"_no_repeats.recode.vcf -- -t q -i 'GT="het" & FMT/AD[:0]/FMT/AD[:1]>'"$allelic_balance_high" -n 0/0 -Oz -o vcf/"$FILE"_tmp1.vcf.gz
      rm vcf/"$FILE"_no_repeats.recode.vcf
    else
      echo not removing sites from repeat regions
      echo setting heterozygous calls biased towards reference to reference allele
      $bcftools +setGT vcf/"$FILE"_raw.vcf.gz -- -t q -i 'GT="het" & FMT/AD[:0]/FMT/AD[:1]>'"$allelic_balance_high" -n 0/0 -Oz -o vcf/"$FILE"_tmp1.vcf.gz
    fi
    #set heterozygous alleles with allelic imbalance biased towards alternate allele (e.g. if set to <0.2, then 1 reference to each 5 alternate allele) to alternate allele; note untested and only compatible with bcftools/1.13 and later
    echo setting heterozygous calls biased towards alternate to alternate allele
    $bcftools +setGT vcf/"$FILE"_tmp1.vcf.gz -- -t q -i 'GT="het" & FMT/AD[:0]/FMT/AD[:1]<'"$allelic_balance_low" -n c:'1/1' -Oz -o vcf/"$FILE"_tmp2.vcf.gz
    rm vcf/"$FILE"_tmp1.vcf.gz
    #set alleles with low read depth to missing
    echo setting calls with low coverage to missing
    $bcftools +setGT vcf/"$FILE"_tmp2.vcf.gz -- -t q -i 'FMT/DP<'"$min_cov" -n ./. -Oz -o vcf/"$FILE"_tmp3.vcf.gz
    rm vcf/"$FILE"_tmp2.vcf.gz
    #set alleles with high read depth to missing; raw QC plots can be used to inform parameter settings
    echo setting calls with high coverage to missing
    $bcftools +setGT vcf/"$FILE"_tmp3.vcf.gz -- -t q -i 'FMT/DP>'"$max_cov" -n ./. -Oz -o vcf/"$FILE"_tmp4.vcf.gz
    rm vcf/"$FILE"_tmp3.vcf.gz
    #set alleles with low GQ (genotyping score; phred scale) to missing
    echo setting calls with low genotyping scores to missing
    $bcftools filter vcf/"$FILE"_tmp4.vcf -i 'FMT/GQ>'"$min_GQ" -S . -Oz -o vcf/"$FILE"_tmp5.vcf.gz
    rm vcf/"$FILE"_tmp4.vcf.gz
    #set alleles within specified distance of indel to missing
    echo setting calls withing specified distance to indels to missing
    $bcftools filter vcf/"$FILE"_tmp5.vcf.gz --SnpGap "$SNPgap":'indel' -S . -Oz -o vcf/"$FILE"_tmp6.vcf.gz
    rm vcf/"$FILE"_tmp5.vcf.gz
    #remove indels from dataset and keep sites with 1-2 alleles (between reference and alternate)
    echo removing indels and alleles with greater than greater than 2 alternate alleles
    $bcftools view vcf/"$FILE"_tmp6.vcf.gz --exclude-types indels --min-alleles 1 --max-alleles 2 -Oz -o vcf/"$FILE"_tmp7.vcf.gz
    rm vcf/"$FILE"_tmp6.vcf.gz
    if [ $full_sites_vcf  -eq  1 ]
    then
      $vcftools --gzvcf vcf/"$FILE"_tmp7.vcf.gz --minQ "$min_Q" --maf "$minor_allele_frequency" --max-missing "$max_missing" --recode --out vcf/"$FILE"_filtered
      $bcftools view vcf/"$FILE"_filtered.recode.vcf -Oz -o vcf/"$FILE"_filtered.vcf.gz
      rm vcf/"$FILE"_filtered.recode.vcf
    else
      $vcftools --gzvcf vcf/"$FILE"_tmp7.vcf.gz --min-alleles 2 --max-alleles 2 --minQ "$min_Q" --maf "$minor_allele_frequency" --max-missing "$max_missing" --recode --out vcf/"$FILE"_filtered
      $bcftools view vcf/"$FILE"_filtered.recode.vcf -Oz -o vcf/"$FILE"_filtered.vcf.gz
      rm vcf/"$FILE"_filtered.recode.vcf
    fi  
    echo check output to determine number of positions masked or filtered at each step
  else
    #remove repeat regions using specified BED file
    echo filtering low coverage vcf
    if [ $full_sites_vcf  -eq  1 ]
    then
      $vcftools --gzvcf vcf/"$FILE"_raw.vcf.gz --exclude-bed reference_genome/$repeat_regions --maf "$minor_allele_frequency" --max-missing "$max_missing" --recode --out vcf/"$FILE"_tmp1
    else
      $vcftools --gzvcf vcf/"$FILE"_raw.vcf.gz --exclude-bed reference_genome/$repeat_regions --min-alleles 2 --max-alleles 2 --maf "$minor_allele_frequency" --max-missing "$max_missing" --recode --out vcf/"$FILE"_tmp1
    fi
    #set -S to . to switch genotypes to a missing value, 0 for reference allele
    $bcftools filter vcf/"$FILE"_tmp1.recode.vcf --SnpGap "$SNPgap":'indel' -S 0 -o vcf/"$FILE"_tmp2.vcf.gz -Oz
    rm vcf/*_tmp1.recode.vcf
    $bcftools view vcf/"$FILE"_tmp2.vcf.gz --exclude-types indels --min-alleles 1 --max-alleles 2  -Oz -o vcf/"$FILE"_filtered.vcf.gz
    rm vcf/*_tmp2.vcf.gz
  fi

  #Check filtered vcf quality to inform filtering parameters
  echo generating filtered vcf quality plots
  $vcftools --gzvcf vcf/"$FILE"_filtered.vcf.gz --depth --out vcf_analyses/vcf_stats/filtered/"$FILE"_filtered
  $vcftools --gzvcf vcf/"$FILE"_filtered.vcf.gz --site-mean-depth --out vcf_analyses/vcf_stats/filtered/"$FILE"_filtered
  $vcftools --gzvcf vcf/"$FILE"_filtered.vcf.gz --freq2 --max-alleles 2 --out vcf_analyses/vcf_stats/filtered/"$FILE"_filtered
  $vcftools --gzvcf vcf/"$FILE"_filtered.vcf.gz --site-quality --out vcf_analyses/vcf_stats/filtered/"$FILE"_filtered
  $vcftools --gzvcf vcf/"$FILE"_filtered.vcf.gz --missing-indv --out vcf_analyses/vcf_stats/filtered/"$FILE"_filtered
  $vcftools --gzvcf vcf/"$FILE"_filtered.vcf.gz --missing-site --out vcf_analyses/vcf_stats/filtered/"$FILE"_filtered
  $vcftools --gzvcf vcf/"$FILE"_filtered.vcf.gz --het --out vcf_analyses/vcf_stats/filtered/"$FILE"_filtered
  cp R_code/02_vcfQC_FILTERED.R R_code/02_vcfQC_"$FILE"_vcfQC_FILTERED.R
  sed -i "s/FILE/$FILE/g" R_code/02_vcfQC_"$FILE"_vcfQC_FILTERED.R
  $R CMD BATCH R_code/02_vcfQC_"$FILE"_vcfQC_FILTERED.R
  rm R_code/02_vcfQC_"$FILE"_vcf_stats_FILTERED.R
  echo user should check filtered vcf plots and, if necessary, modify filtering presets and rerun workflow
else
  echo skipping vcf filtration
  echo user is satisfied with filtering results, moving onto PCA and admixture analyses
  echo user should consider adding another filtering step to remove highly heterozygous sites, with example commands provided below
  echo use populations module of stacks to calculate heterozygosity per site defining all individuals are single population
  echo populations -V <in>.vcf.gz -M <population_structure>.txt -t <threads> --out-path .
  echo decide on a threshold for acceptable site heterozygosity and remove sites above threshold
  echo vcftools --gzvcf <in>.vcf.gz --exclude-positions <heterozygous_sites>.txt --recode --out <out>
fi

#####Generate PCA plots
#Extract unlinked loci for PCA and ADMIXTURE analyses
echo running plink to prune linked sites
#identify linked loci, --indep-pairwise sets parameters of interest: window size for calculating r2 in kb, step size for moving along reference scaffolds, and r2 limit, above which loci are pruned
$plink --vcf vcf/"$FILE"_filtered.vcf.gz --double-id --allow-extra-chr --set-missing-var-ids @:#\$1,\$2 --indep-pairwise "$plink_window_size"'kb' 1 "$plink_r2" --out vcf_analyses/plink/"$FILE"
#Extract unlinked loci for ADMIXTURE analyses
#If looking to bootstrap admixture analyses, add--recode 12 and use .ped and .map as input for admixture
$plink --vcf vcf/"$FILE"_filtered.vcf.gz --double-id --allow-extra-chr --set-missing-var-ids @:#\$1,\$2 --extract vcf_analyses/plink/"$FILE".prune.in --make-bed --pca --out vcf_analyses/plink/"$FILE"
#convert plink output to vcf for phylogenetic analysis of unlinked SNPs
echo outputing pruned vcf
$plink --bfile vcf_analyses/plink/"$FILE" --allow-extra-chr --recode vcf --out vcf_analyses/plink/"$FILE".LD_pruned
$bcftools view vcf_analyses/plink/"$FILE".LD_pruned -Oz -o vcf_analyses/plink/"$FILE".LD_pruned.vcf.gz
#fixes columns names so admixture doesn't expect human chromosomes
awk '{$1=0;print $0}' vcf_analyses/plink/"$FILE".bim > vcf_analyses/plink/"$FILE".bim.tmp
mv vcf_analyses/plink/"$FILE".bim.tmp vcf_analyses/plink/"$FILE".bim
#Get sample IDs, with some clunky intermediate files
$bcftools query -l vcf/"$FILE"_filtered.vcf.gz > vcf_analyses/admixture/vcf.list.1
sed 's/\./\ /g' vcf_analyses/admixture/vcf.list.1 > vcf_analyses/admixture/sample.list.2
sed 's/_n_s.bam//' vcf_analyses/admixture/sample.list.2 > vcf_analyses/admixture/sample.list.3
cut -f 2 -d ' ' vcf_analyses/admixture/sample.list.3 > vcf_analyses/admixture/sample.list.4
##plot PCA of results
##Get list of species/populations to plot
sort vcf_analyses/admixture/sample.list.4 | uniq > vcf_analyses/admixture/pop-species.1.list
#A quirky workaround to get R looping through populations
paste -d, -s vcf_analyses/admixture/pop-species.1.list > vcf_analyses/admixture/pop-species.2.list
cat vcf_analyses/admixture/pop-species.2.list > vcf_analyses/admixture/pop-species.3.list
sed ' 1 s/.*/&,NA/' vcf_analyses/admixture/pop-species.2.list > vcf_analyses/admixture/pop-species.4.list
#PCA plot using plink output
echo generating PCA plots for PCs 1-6
cp R_code/03_PCA_plots.R R_code/03_"$FILE"_PCA_plots.R
sed -i "s/FILE/$FILE/g" R_code/03_"$FILE"_PCA_plots.R
R CMD BATCH R_code/03_"$FILE"_PCA_plots.R
rm R_code/03_"$FILE"_PCA_plots.R

##Admixture
echo running admixture
#Add -B flag to bootstrap admixture analysis
for i in $(seq $k_low $k_high)
do
$admixture --cv vcf_analyses/plink/"$FILE".bed $i -j"$threads" > vcf_analyses/admixture/log${i}.out
mv *.P vcf_analyses/admixture
mv *.Q vcf_analyses/admixture
done
#collect cross-validation errors to determine best value of k
awk '/CV/ {print $3,$4}' vcf_analyses/admixture/*out | cut -c 4,7-20 > vcf_analyses/admixture/"$FILE".cv.error
awk '{split($1,name,"."); print $1,name[2]}' vcf_analyses/plink/${FILE}.nosex > vcf_analyses/admixture/"$FILE".list
#Get sample IDs, with some clunky intermediate files
$bcftools query -l vcf/"$FILE"_filtered.vcf.gz > vcf_analyses/admixture/vcf.list.1
sed 's/\./ /g' vcf_analyses/admixture/vcf.list.1 > vcf_analyses/admixture/sample.list.2
awk '{print $2}' vcf_analyses/admixture/sample.list.2 > vcf_analyses/admixture/sample.list.3
paste vcf_analyses/admixture/vcf.list.1 vcf_analyses/admixture/sample.list.3 > vcf_analyses/admixture/ind.list
#plot ancestry proportions at different values of k (admixture results)
echo The R script provided here was written by Joana Meier, September 2019, modified by Trevor Bringloe March 2022
echo https://github.com/speciationgenomics/scripts/blob/master/plotADMIXTURE.r
$Rscript --vanilla R_code/04_admixture_plots.R -p $FILE -i vcf_analyses/admixture/ind.list -k "$k_high" -r $RBpalette -l $pops
echo admixture plots 'done'
echo workflow completed :)
echo see end of bash script for commands to run pixy
exit 0

##Pixy
#optional commands to fetch nucleotide diversity and fst values from an invariant+variant vcf
#must setup conda environment <pixy> before running commands and create populations file for defining structure
tabix vcf/"$FILE"_invariant+variant_filtered.vcf.gz
conda activate pixy
pixy --stats pi dxy fst --vcf vcf/"$FILE"_invariant+variant_filtered.vcf.gz --populations MOBELS_pops_pixy_16iv23.txt --window size 50000 --n_cores 16 --output_prefix "$FILE"

###########################################################################
############################  END WORKFLOW  ###############################
###########################################################################
