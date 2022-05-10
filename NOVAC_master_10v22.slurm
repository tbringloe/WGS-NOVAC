#!/bin/bash
##Commands to produce a Variant Call Format (VCF) file to store SNPs+indels for downstream analysis
##Script designed to process forward and reverse short read (150 bp) sequences #If multiple reads files are present, concatenate into single files using some variation of cat <sample_ID>_<*wild card for different libraries>_1.fastq.gz > <sample_ID>.fastq # Note NOVOPlasty does not process single end sequences
##User must provide a reference genome for the target species located /reference_genome, or a closely related species if tailoring for phylogeographic analysis
##User must provide seed sequences for mitochondrial and plastid sequences located in NOVOPlasty/seed_files; note, these don't have to be from the target species, settings are set to not directly extend seed sequence.
##User must provide a single column sample list of sample IDs (sample.list)
##Script was written with the intention of running as a slurm script # can we wrap this so it runs on private lab servers?
##!!!Script assumes the following naming convention for read files, ordering the samples by population in a relevant order: 01-SampleID_Population/species_1/2.fq.gz, 02_SampleID_Population/species_1/2.fq.gz, 03_SampleID_Population/species_1/2.fq.gz...
##Analysis assumes sequencing effort is sufficient in all specimens; check raw QC of vcf file, and optionally include a step to remove individuals from the analysis
##Dependencies are extensive; NOVOPlasty; Perl; java; trimmomatic; bowtie2; samtools; bcftools; vcftools; Plink; Admixture; R; CRAN packages PopGenome, qqman, ggplot2, tidyverse, RColorBrewer; PopDecay; RAiSD;
#######################################SLURM parameters#########################################################
#####################swap this framework out for respective HPC environment#####################################
# Created by the University of Melbourne job script generator for SLURM
# Mon Dec 09 2019 15:50:02 GMT+1100 (Australian Eastern Daylight Time)

# Partition for the job:
#SBATCH --partition=physical

# Multithreaded (SMP) job: must run on one node and the cloud partition
#SBATCH --nodes=1

# The name of the job:
#SBATCH --job-name="Variant_dataset_build_test"

# The project ID which this job should run under:
#SBATCH --account=

# Maximum number of tasks/CPU cores used by the job:
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=1

# The amount of memory in megabytes per process in the job:
#SBATCH --mem=50000

# Send yourself an email when the job:
# aborts abnormally (fails)
#SBATCH --mail-type=FAIL
# begins
#SBATCH --mail-type=BEGIN
# ends successfully
#SBATCH --mail-type=END

# Use this email address:
#SBATCH --mail-user=

# The maximum running time of the job in days-hours:mins:sec
#SBATCH --time=0-11:59:59

# check that the script is launched with sbatch
if [ "x$SLURM_JOB_ID" == "x" ]; then
   echo "You need to submit your job to the queuing system with sbatch"
   exit 1
fi

# Run the job from the directory where it was launched (default)

# The modules to load:

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
min_Q=30 #filters entire sites according to the probability there is an alternate allele present, follows phred scale, 30=1/1000 chance of SNP calling error
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
################################################PIPELINE############################################################################################################################################################
##Make sure your work environment is setup, see READme.txt
##Index reference genomes for mapping
module load gcc/10.2.0
module load bowtie2/2.4.2
bowtie2-build reference_genome/$reference_nuclear_genome IDX/reference_nuclear_genome
bowtie2-build reference_genome/$reference_mitochondrial_genome IDX/reference_mitochondrial_genome
bowtie2-build reference_genome/$reference_chloroplast_genome IDX/reference_chloroplast_genome

#Perform QC on raw files in parallel using xargs command
module purge
module load fastqc/0.11.9-java-11
module load trimmomatic/0.39-java-11.0.2
cat sample.list | xargs -I {} -n 1 -P $threads sh -c "fastqc read_files/{}_1.fq.gz read_files/{}_2.fq.gz -o results/read_QC/raw"
#trim raw read files based on specified criteria; trimmed files are deposited into tmp/trimming before moving into tmp as an annoying workaround, trimmomatic otherwise deletes and overwrites files in the loop
cat sample.list | xargs -I {} -n 1 -P $threads sh -c "java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE read_files/{}_1.fq.gz read_files/{}_2.fq.gz tmp/trimming/{}_1_tp.fq.gz tmp/trimming/{}_1_tup.fq.gz tmp/trimming/{}_2_tp.fq.gz tmp/trimming/{}_2_tup.fq.gz TRAILING:$Trailing HEADCROP:$Headcrop AVGQUAL:$Average_quality MINLEN:$Min_length"
cat sample.list | while read line
do
mv tmp/trimming/"$line"_1_tp.fq.gz tmp/"$line"_trimmed_1.fq.gz
mv tmp/trimming/"$line"_2_tp.fq.gz tmp/"$line"_trimmed_2.fq.gz
done
cat sample.list | xargs -I {} -n 1 -P $threads sh -c "fastqc tmp/{}_trimmed_1.fq.gz tmp/{}_trimmed_2.fq.gz -o results/read_QC/trimmed"

#Summarize QC, be sure to evaluate
module purge
module load foss/2019b
module load multiqc/1.9-python-3.7.4
multiqc results/read_QC/raw -o results/read_QC/raw
multiqc results/read_QC/trimmed -o results/read_QC/trimmed

#################################Optional stream to assemble the organellar genomes using NOVOPlasty###############################################################################################################
##Mito
#Extract subreads in parallel for assembly, since we don't typically need all the reads to assemble the organelles; seqtk outputs uncompressed reads, which we need anyways for NOVOPlasty; -s flag specifies a random seed, which we need to grab same forward and reverse reads
module purge
module load foss/2019b
module load seqtk/1.3
cat sample.list | xargs -I {} -n 1 -P $threads sh -c "seqtk sample -s100 tmp/{}_trimmed_1.fq.gz "$Number_reads_to_extract" > tmp/{}_sub_1.fq"
cat sample.list | xargs -I {} -n 1 -P $threads sh -c "seqtk sample -s100 tmp/{}_trimmed_2.fq.gz "$Number_reads_to_extract" > tmp/{}_sub_2.fq"
#setup NOVOPlasty environment to run multiple samples in parallel
cat sample.list | while read line
do
#setup file paths for analysis
mkdir NOVOPlasty/mito/$line
#Replace specified parameters in config file
cp NOVOPlasty/config_mito.txt NOVOPlasty/config_mito_"$line".txt
sed -i "s/PROJECT_NAME/"$line"/g" NOVOPlasty/config_mito_"$line".txt
sed -i "s/GENOME_RANGE_MITO/$Genome_range_mito/g" NOVOPlasty/config_mito_"$line".txt
sed -i "s/KMER/$Kmer/g" NOVOPlasty/config_mito_"$line".txt
sed -i "s/MITO_SEED/$mito_seed/g" NOVOPlasty/config_mito_"$line".txt
sed -i "s/READ_LENGTH/$Read_length/g" NOVOPlasty/config_mito_"$line".txt
sed -i "s/INSERT_SIZE/$Insert_size/g" NOVOPlasty/config_mito_"$line".txt
sed -i "s/FORWARD/"$line"_sub_1.fq/g" NOVOPlasty/config_mito_"$line".txt
sed -i "s/REVERSE/"$line"_sub_2.fq/g" NOVOPlasty/config_mito_"$line".txt
mv NOVOPlasty/config_mito_"$line".txt NOVOPlasty/mito/$line
cp NOVOPlasty/filter_reads.pl NOVOPlasty/mito/$line
cp NOVOPlasty/LICENSE NOVOPlasty/mito/$line
cp -a NOVOPlasty/NOVOPlasty4.2.pl NOVOPlasty/mito/$line
cp NOVOPlasty/README.md NOVOPlasty/mito/$line
done
#Run NOVOPlasty in parallel for multiple samples
module purge
module load gcccore/11.2.0
module load perl/5.34.0
cat sample.list | xargs -I {} -n 1 -P $threads sh -c "NOVOPlasty/mito/{}/NOVOPlasty4.2.pl -c NOVOPlasty/mito/{}/config_mito_{}.txt"
#Output is in working directory, so move to relevant sample folder
cat sample.list | while read line
do
mv *$line* NOVOPlasty/mito/$line
done

##Plastid## Can remove this portion of the script if not a photosynthesic organism
#setup NOVOPlasty environment to run multiple samples in parallel
cat sample.list | while read line
do
#setup file paths for analysis
mkdir NOVOPlasty/plastid/$line
#Replace specified parameters in config file
cp NOVOPlasty/config_plastid.txt NOVOPlasty/config_plastid_"$line".txt
sed -i "s/PROJECT_NAME/"$line"/g" NOVOPlasty/config_plastid_"$line".txt
sed -i "s/GENOME_RANGE_PLASTID/$Genome_range_plastid/g" NOVOPlasty/config_plastid_"$line".txt
sed -i "s/KMER/$Kmer/g" NOVOPlasty/config_plastid_"$line".txt
sed -i "s/PLASTID_SEED/"$plastid_seed"/g" NOVOPlasty/config_plastid_"$line".txt
sed -i "s/READ_LENGTH/$Read_length/g" NOVOPlasty/config_plastid_"$line".txt
sed -i "s/INSERT_SIZE/$Insert_size/g" NOVOPlasty/config_plastid_"$line".txt
sed -i "s/FORWARD/"$line"_sub_1.fq/g" NOVOPlasty/config_plastid_"$line".txt
sed -i "s/REVERSE/"$line"_sub_2.fq/g" NOVOPlasty/config_plastid_"$line".txt
mv NOVOPlasty/config_plastid_"$line".txt NOVOPlasty/plastid/$line
cp NOVOPlasty/filter_reads.pl NOVOPlasty/plastid/$line
cp NOVOPlasty/LICENSE NOVOPlasty/plastid/$line
cp -a NOVOPlasty/NOVOPlasty4.2.pl NOVOPlasty/plastid/$line
cp NOVOPlasty/README.md NOVOPlasty/plastid/$line
done
#Run NOVOPlasty
cat sample.list | xargs -I {} -n 1 -P $threads sh -c "NOVOPlasty/plastid/{}/NOVOPlasty4.2.pl -c NOVOPlasty/plastid/{}/config_plastid_{}.txt"
#Output is in working directory, so move to relevant sample folder
cat sample.list | while read line
do
mv *$line* NOVOPlasty/plastid/$line
done

###############################################Stream to call variant positions on organellar genomes#####################################################################################################################################
#mitochondrial variant positions
module purge
module load gcc/10.2.0
module load bowtie2/2.4.2
cat sample.list | while read line
do
bowtie2 --score-min L,0,-"$map_param_mito" --no-unal -p $threads -x IDX/reference_mitochondrial_genome -1 tmp/"$line"_sub_1.fq -2 tmp/"$line"_sub_2.fq -S tmp/"$line"_m.sam
done

#Sort, convert, and index Sequence Alignment Map to binary format for compiling step; steps are broken into multiple loops because environment dependencies are not compatible
module load samtools/1.12
module load bcftools/1.12
cat sample.list | while read line
do
samtools view -S -@ "$threads" -b tmp/"$line"_m.sam > tmp/"$line"_m.bam
#File paths get incorporated into vcf file downstream, so sorted bams are deposited into root/launch directory and moved later
samtools sort -@ "$threads" -o "$line"_m_s.bam tmp/"$line"_m.bam
samtools index -@ "$threads" "$line"_m_s.bam
rm tmp/"$line"_m.bam
rm tmp/"$line"_m.sam
done

#Compile BAM files into VCF file; this step takes a long time depending on the size of the dataset, the process is done in parallel here by dividing the task up by contig
ls *_m_s.bam > sample_list1
ls *_m_s.bam > sample_list2
sed -i "s/$/\    1/" sample_list2
#index the reference genome for mpileup
samtools faidx reference_genome/$reference_mitochondrial_genome
#mpileup and call SNPs
bcftools mpileup --threads $threads -Ou -f reference_genome/$reference_mitochondrial_genome -a AD,ADF,ADR,DP -b sample_list1 | bcftools call --threads $threads --variants-only -S sample_list2 -f GQ,GP -mv -Oz > vcf/"$FILE"_raw_mito.vcf.gz

#Apply filtering criteria
module purge
module load gcc/10.3.0
module load bcftools/1.15
##Old option that works for earlier versions of bcftools (earlier than 1.13), but looses information by setting alternate alleles to missing; not recommended for organellar genomes as this will bias sites towards the reference genome, which will have implications in the downstream analyses. Be sure to inspect an alignment of the genomes for problematic mapping areas.
#set heterozygous alleles with allelic imbalance biased towards reference allele (e.g. if set to 5, then 5 reference to each 1 alternate allele) to reference
#bcftools +setGT vcf/"$FILE"_raw_mito.vcf.gz -- -t q -i 'GT="het" & FMT/AD[:0]/FMT/AD[:1]>'"$allelic_balance_high" -n 0/0 > vcf/"$FILE"_tmp1.vcf
#set heterozygous alleles with allelic imbalance biased towards alternate allele (e.g. if set to <0.2, then 1 reference to each 5 alternate allele) to missing
#bcftools +setGT vcf/"$FILE"_tmp1.vcf -- -t q -i 'GT="het" & FMT/AD[:0]/FMT/AD[:1]<'"$allelic_balance_low" -n ./. > vcf/"$FILE"_tmp2.vcf

#A better set of commands that switches imbalanced allelic depth for heterozygous sites to more represented homozygous genotype (reference or alternate), leverages new syntax options of bctools/1.13 and later, tested here on v.1.15
#set heterozygous alleles with allelic imbalance biased towards reference allele (e.g. if set to 5, then 5 reference to each 1 alternate allele) to reference allele
bcftools +setGT vcf/"$FILE"_raw_mito.vcf.gz -- -t q -i 'GT="het" & FMT/AD[:0]/FMT/AD[:1]>'"$allelic_balance_high" -n 0/0 > vcf/"$FILE"_tmp1.vcf
#set heterozygous alleles with allelic imbalance biased towards alternate allele (e.g. if set to <0.2, then 1 reference to each 5 alternate allele) to alternate allele; note untested and only compatible with bcftools/1.13 and later
bcftools +setGT vcf/"$FILE"_tmp1.vcf -- -t q -i 'GT="het" & FMT/AD[:0]/FMT/AD[:1]<'"$allelic_balance_low" -n c:'1/1' > vcf/"$FILE"_tmp2.vcf

#set alleles with low read depth to missing
bcftools +setGT vcf/"$FILE"_tmp2.vcf -- -t q -i 'FMT/DP<'"$min_cov" -n ./. > vcf/"$FILE"_tmp3.vcf
##Filter sites for Q values [quality scores]; sites with any missing data and indels are removed to ensure SNPs are kept in alignment for downstream analysis
#set alleles with low GQ (genotyping score; phred scale) to missing
bcftools +setGT vcf/"$FILE"_tmp3.vcf -- -t q -i 'FMT/GQ<'"$min_GQ" -n ./. > vcf/"$FILE"_tmp4.vcf

module purge
module load gcc/8.3.0
module load vcftools/0.1.16
vcftools --gzvcf vcf/"$FILE"_tmp4.vcf --max-missing 1.0 --max-alleles 2 --minQ "$min_Q" --remove-indels --recode --out vcf/"$FILE"_mito
rm /vcf/*tmp*

#generate fasta of concensus sequences for each sample
module purge
module load gcc/10.3.0
module load bcftools/1.15
#compress and final vcf files
bgzip vcf/"$FILE"_mito.recode.vcf
bcftools index vcf/"$FILE"_mito.recode.vcf.gz
ls *_m_s.bam > sample_list_mito
cat sample_list_mito | while read line
do
bcftools consensus -I -s $line -f reference_genome/$reference_mitochondrial_genome vcf/"$FILE"_mito.recode.vcf.gz > results/organellar/mito/mapping/"$line"_map_mito.fasta
done
#swap out sample IDs in fasta files for PopGenome analysis downstream
cat sample.list | while read line
do
sed -i 1s/^.*$/'>'"$line"/ results/organellar/mito/mapping/"$line"_m_s.bam_map_mito.fasta
done
cat results/organellar/mito/mapping/*.fasta > results/organellar/Popgenome/mito/"$FILE"_mito.fasta

#move sorted bam files to sorted_bam
mv *_m_s.bam sorted_bam
mv *m_s.bam.bai sorted_bam

##chloroplast variant positions
module purge
module load gcc/10.2.0
module load bowtie2/2.4.2
cat sample.list | while read line
do
bowtie2 --score-min L,0,-"$map_param_chloro" --no-unal -p $threads -x IDX/reference_chloroplast_genome -1 tmp/"$line"_sub_1.fq -2 tmp/"$line"_sub_2.fq -S tmp/"$line"_p.sam
done

#Sort, convert, and index Sequence Alignment Map to binary format for compiling step; steps are broken into multiple loops because environment dependencies are not compatible
module load samtools/1.12
module load bcftools/1.12
cat sample.list | while read line
do
samtools view -S -@ "$threads" -b tmp/"$line"_p.sam > tmp/"$line"_p.bam
#File paths get incorporated into vcf file downstream, so sorted bams are deposited into root/launch directory and moved later
samtools sort -@ "$threads" -o "$line"_p_s.bam tmp/"$line"_p.bam
samtools index -@ "$threads" "$line"_p_s.bam
rm tmp/"$line"_p.bam
rm tmp/"$line"_p.sam
done

#Compile BAM files into VCF file; this step takes a long time depending on the size of the dataset, the process is done in parallel here by dividing the task up by contig
ls *_p_s.bam > sample_list1
ls *_p_s.bam > sample_list2
sed -i "s/$/\    1/" sample_list2
#index the reference genome for mpileup
samtools faidx reference_genome/$reference_chloroplast_genome
#mpileup and call SNPs
bcftools mpileup --threads $threads -Ou -f reference_genome/$reference_chloroplast_genome -a AD,ADF,ADR,DP -b sample_list1 | bcftools call --threads $threads --variants-only -S sample_list2 -f GQ,GP -mv -Oz > vcf/"$FILE"_raw_chloro.vcf.gz

#Apply filtering criteria
module purge
module load gcc/10.3.0
module load bcftools/1.15
##Old option that works for earlier versions of bcftools (earlier than 1.13), but looses information by setting alternate alleles to missing; not recommended for organellar genomes as this will bias sites towards the reference genome, which will have implications in the downstream analyses. Be sure to inspect an alignment of the genomes for problematic mapping areas.
#set heterozygous alleles with allelic imbalance biased towards reference allele (e.g. if set to 5, then 5 reference to each 1 alternate allele) to reference
#bcftools +setGT vcf/"$FILE"_raw_chloro.vcf.gz -- -t q -i 'GT="het" & FMT/AD[:0]/FMT/AD[:1]>'"$allelic_balance_high" -n 0/0 > vcf/"$FILE"_tmp1.vcf
#set heterozygous alleles with allelic imbalance biased towards alternate allele (e.g. if set to <0.2, then 1 reference to each 5 alternate allele) to missing
#bcftools +setGT vcf/"$FILE"_tmp1.vcf -- -t q -i 'GT="het" & FMT/AD[:0]/FMT/AD[:1]<'"$allelic_balance_low" -n ./. > vcf/"$FILE"_tmp2.vcf

#A better set of commands that switches imbalanced allelic depth for heterozygous sites to more represented homozygous genotype (reference or alternate), leverages new syntax options of bctools/1.13 and later, tested here on v.1.15
#set heterozygous alleles with allelic imbalance biased towards reference allele (e.g. if set to 5, then 5 reference to each 1 alternate allele) to reference allele
bcftools +setGT vcf/"$FILE"_raw_chloro.vcf.gz -- -t q -i 'GT="het" & FMT/AD[:0]/FMT/AD[:1]>'"$allelic_balance_high" -n 0/0 > vcf/"$FILE"_tmp1.vcf
#set heterozygous alleles with allelic imbalance biased towards alternate allele (e.g. if set to <0.2, then 1 reference to each 5 alternate allele) to alternate allele; note untested and only compatible with bcftools/1.13 and later
bcftools +setGT vcf/"$FILE"_tmp1.vcf -- -t q -i 'GT="het" & FMT/AD[:0]/FMT/AD[:1]<'"$allelic_balance_low" -n c:'1/1' > vcf/"$FILE"_tmp2.vcf

#set alleles with low read depth to missing
bcftools +setGT vcf/"$FILE"_tmp2.vcf -- -t q -i 'FMT/DP<'"$min_cov" -n ./. > vcf/"$FILE"_tmp3.vcf
##Filter sites for Q values [quality scores] and missing data above specified threshold; sites with any missing data and indels are removed to ensure SNPs are kept in alignment for downstream analysis
#set alleles with low GQ (genotyping score; phred scale) to missing
bcftools +setGT vcf/"$FILE"_tmp3.vcf -- -t q -i 'FMT/GQ<'"$min_GQ" -n ./. > vcf/"$FILE"_tmp4.vcf

module purge
module load gcc/8.3.0
module load vcftools/0.1.16
vcftools --gzvcf vcf/"$FILE"_tmp4.vcf --max-missing 1.0 --max-alleles 2 --minQ "$min_Q" --remove-indels --recode --out vcf/"$FILE"_chloro
rm /vcf/*tmp*

#generate fasta of concensus sequences for each sample
#index and compress final vcf files
module purge
module load gcc/10.3.0
module load bcftools/1.15
bgzip vcf/"$FILE"_chloro.recode.vcf
bcftools index vcf/"$FILE"_chloro.recode.vcf.gz
ls *_p_s.bam > sample_list_chloro
cat sample_list_chloro | while read line
do
bcftools consensus -I -s $line -f reference_genome/$reference_chloroplast_genome vcf/"$FILE"_chloro.recode.vcf.gz > results/organellar/plastid/mapping/"$line"_map_chloro.fasta
done
#swap out sample IDs in fasta files for PopGenome analysis downstream
cat sample.list | while read line
do
sed -i 1s/^.*$/'>'"$line"/ results/organellar/plastid/mapping/"$line"_p_s.bam_map_chloro.fasta
done
cat results/organellar/plastid/mapping/*.fasta > results/organellar/Popgenome/chloro/"$FILE"_plastid.fasta

#move sorted bam files to sorted_bam
mv *_p_s.bam sorted_bam
mv *p_s.bam.bai sorted_bam
#Remove subsampled reads and sorted bam files
#rm tmp/*sub_1.fq
#rm tmp/*sub_2.fq
#rm sorted_bam/*bam

###########WARNING#####################WARNING###################WARNING###################WARNING#####################WARNING#####################WARNING#########################WARNING###############################WARNING######################
##############Do not publish organellar genome sequences based on read mapping, any sites without data will be force called to the reference genome, causing spurious calls, examine de novo assemblies for publishing genomes########################
#################This will happen for any regions with higher (than the specified) rate of sequence evolution, or otherwise potentially difficult to align portions of the genomes, please inspect the alignment of genomes###########################
#######################################Mapped genomes are conceptually ok for this workflow but do not represent the true sequences of the organellar genomes#########################################################################################
######################################################################################################################################################################################################################################################

#############################################################Map reads to nuclear reference genome, then convert to sorted bam file###################################################################################################################
module purge
module load gcc/10.2.0
module load bowtie2/2.4.2
cat sample.list | while read line
do
bowtie2 --score-min L,0,-"$map_param_nuc" --no-unal -p $threads -x IDX/reference_nuclear_genome -1 tmp/"$line"_trimmed_1.fq.gz -2 tmp/"$line"_trimmed_2.fq.gz -S tmp/"$line"_n.sam
done

#Sort, convert, and index Sequence Alignment Map to binary format for compiling step; steps are broken into multiple loops because environment dependencies are not compatible
module load samtools/1.12
module load bcftools/1.12
cat sample.list | while read line
do
samtools view -S -@ "$threads" -b tmp/"$line"_n.sam > tmp/"$line"_n.bam
#File paths get incorporated into vcf file downstream, so sorted bams are deposited into root/launch directory and moved later
samtools sort -@ "$threads" -o "$line"_n_s.bam tmp/"$line"_n.bam
#If you are struggling to manage space on the server, at this point you could delete all the trimmed and subbed reads, and all the intermediate alignment files. Default here is to not delete, in case there are errors in the workflow and the job must restart with these files in place.
#rm tmp/$line*
samtools index -@ "$threads" "$line"_n_s.bam
rm tmp/"$line"_n.bam
rm tmp/"$line"_n.sam
done

#Compile BAM files into VCF file; this step takes a long time depending on the size of the dataset, the process is done in parallel here by dividing the task up by contig
ls *_n_s.bam > sample_list1
ls *_n_s.bam > sample_list2
sed -i "s/$/\    "$ploidy"/" sample_list2
#index the reference genome for mpileup
samtools faidx reference_genome/$reference_nuclear_genome
#mpileup and call on each contig, specified parallel pileups by creating xargs list of each contig. Use samtools view to obtain a list of all contig names, and input each as a xargs command, with -P as number of parallel jobs
samtools view -@ "$threads" -H "$EXAMPLE_SAMPLE"_n_s.bam | grep "\@SQ" | sed 's/^.*SN://g' | cut -f 1 | xargs -I {} -n 1 -P $threads sh -c "bcftools mpileup --threads $threads -Ou -f reference_genome/$reference_nuclear_genome -a AD,ADF,ADR,DP -r {} -b sample_list1 | bcftools call --threads $threads --variants-only -S sample_list2 -f GQ,GP -mv -Oz > tmp/tmp.{}.vcf.gz"
#list temporary vcf files
ls tmp/tmp* > list_tmp
#concat into single vcf
bcftools concat --threads $threads -f list_tmp -Oz > vcf/"$FILE"_raw_nuc.vcf.gz
rm tmp/tmp*
rm list_tmp
mv *_n_s.bam sorted_bam
mv *n_s.bam.bai sorted_bam

#########################Congrats, the computationally intensive part is over, consider running the remaining script with less threads if the job hangs up or you need to optmize parameters#########################################

module purge
module load gcc/8.3.0
module load vcftools/0.1.16
##QC plots of the raw vcf dataset; aspects of the plots can and should be altered in the R script; please restart pipeline from here if initial thresholds can be optimized
vcftools --gzvcf vcf/"$FILE"_raw_nuc.vcf.gz --depth --out vcf_analyses/QC_RAW/$FILE
vcftools --gzvcf vcf/"$FILE"_raw_nuc.vcf.gz --site-mean-depth --out vcf_analyses/QC_RAW/$FILE
vcftools --gzvcf vcf/"$FILE"_raw_nuc.vcf.gz --freq2 --max-alleles 2 --out vcf_analyses/QC_RAW/$FILE
vcftools --gzvcf vcf/"$FILE"_raw_nuc.vcf.gz --site-quality --out vcf_analyses/QC_RAW/$FILE
vcftools --gzvcf vcf/"$FILE"_raw_nuc.vcf.gz --missing-indv --out vcf_analyses/QC_RAW/$FILE
vcftools --gzvcf vcf/"$FILE"_raw_nuc.vcf.gz --missing-site --out vcf_analyses/QC_RAW/$FILE
vcftools --gzvcf vcf/"$FILE"_raw_nuc.vcf.gz --het --out vcf_analyses/QC_RAW/$FILE
cp R_code/01_vcfQC_RAW.R R_code/01_vcfQC_"$FILE"_RAW.R
sed -i "s/FILE/$FILE/g" R_code/01_vcfQC_"$FILE"_RAW.R
module purge
module load foss/2020b
module load fosscuda/2020b
module load r/4.0.4
R CMD BATCH R_code/01_vcfQC_"$FILE"_RAW.R
rm R_code/01_vcfQC_"$FILE"_RAW.R

##Apply filtering criteria
module purge
module load gcc/10.3.0
module load bcftools/1.15
##Old option that works for earlier versions of bcftools (earlier than 1.13), but looses information by setting alternate alleles to missing
#set heterozygous alleles with allelic imbalance biased towards reference allele (e.g. if set to 5, then 5 reference to each 1 alternate allele) to reference
#bcftools +setGT vcf/"$FILE"_raw_nuc.vcf.gz -- -t q -i 'GT="het" & FMT/AD[:0]/FMT/AD[:1]>'"$allelic_balance_high" -n 0/0 > vcf/"$FILE"_tmp1.vcf
#set heterozygous alleles with allelic imbalance biased towards alternate allele (e.g. if set to <0.2, then 1 reference to each 5 alternate allele) to missing
#bcftools +setGT vcf/"$FILE"_tmp1.vcf -- -t q -i 'GT="het" & FMT/AD[:0]/FMT/AD[:1]<'"$allelic_balance_low" -n ./. > vcf/"$FILE"_tmp2.vcf

#A better set of commands that switches imbalanced allelic depth for heterozygous sites to more represented homozygous genotype (reference or alternate), leverages new syntax options of bctools/1.13 and later, tested here on v.1.15
#set heterozygous alleles with allelic imbalance biased towards reference allele (e.g. if set to 5, then 5 reference to each 1 alternate allele) to reference allele
bcftools +setGT vcf/"$FILE"_raw_nuc.vcf.gz -- -t q -i 'GT="het" & FMT/AD[:0]/FMT/AD[:1]>'"$allelic_balance_high" -n 0/0 > vcf/"$FILE"_tmp1.vcf
#set heterozygous alleles with allelic imbalance biased towards alternate allele (e.g. if set to <0.2, then 1 reference to each 5 alternate allele) to alternate allele; note untested and only compatible with bcftools/1.13 and later
bcftools +setGT vcf/"$FILE"_tmp1.vcf -- -t q -i 'GT="het" & FMT/AD[:0]/FMT/AD[:1]<'"$allelic_balance_low" -n c:'1/1' > vcf/"$FILE"_tmp2.vcf

#set alleles with low read depth to missing
bcftools +setGT vcf/"$FILE"_tmp2.vcf -- -t q -i 'FMT/DP<'"$min_cov" -n ./. > vcf/"$FILE"_tmp3.vcf
##set alleles with high read depth to missing; note raw QC plots can be used to inform parameter settings
bcftools +setGT vcf/"$FILE"_tmp3.vcf -- -t q -i 'FMT/DP>'"$max_cov" -n ./. > vcf/"$FILE"_tmp4.vcf
#set alleles with low GQ (genotyping score; phred scale) to missing
bcftools +setGT vcf/"$FILE"_tmp4.vcf -- -t q -i 'FMT/GQ<'"$min_GQ" -n ./. > vcf/"$FILE"_tmp5.vcf

#optional step to remove samples with low overall coverage, i.e. sequencing effort is insufficient in these individuals and they should be removed before downstream analysis; replace XXX with samples names in vcf file (will be the sorted bam file name, e.g. XXX_n_s.bam), --remove-indv XXX can be specified multiple times for individual samples; if using this command, remove the hashtag prior to the command line and change $FILE_tmp4.vcf to $FILE_tmp4.1.recode.vcf in the next command step
module purge
module load gcc/8.3.0
module load vcftools/0.1.16
#vcftools --gzvcf vcf/"$FILE"_tmp5.vcf --remove-indv XXX --recode --out vcf/"$FILE"_4.1
#optional step to remove sites from low complexity and repeat regions; bed file can be created from output of repeat analysis when annotating the reference genome; if using this command, remove the hashtag prior to the command line and change $FILE_tmp4.vcf or $FILE_tmp4.1.recode.vcf to $FILE_tmp4.1.recode.vcf or $FILE_tmp4.2.recode.vcf in the next command step
vcftools --gzvcf vcf/"$FILE"_tmp4.vcf --exclude-bed vcf/$repeat_regions --recode --out vcf/"$FILE"_tmp4.2
##Filter sites for Q values [quality scores] and missing data above specified threshold, and remove monomorphic sites and sites with >2 alleles
vcftools --gzvcf vcf/"$FILE"_tmp4.2.recode.vcf --min-alleles 1 --max-alleles 2 --max-missing "$max_missing" --minQ "$min_Q" --recode --out vcf/"$FILE"_noMAFapplied
##Use vcftools to remove sites with minor allele frequency specified by the user; note, some pop statistics will be performed on dataset not filtered for MAF
vcftools --gzvcf vcf/"$FILE"_noMAFapplied.recode.vcf --maf "$minor_allele_frequency" --recode --out vcf/"$FILE"_MAFapplied
#compress final vcf files
bgzip vcf/"$FILE"_noMAFapplied.recode.vcf
bgzip vcf/"$FILE"_MAFapplied.recode.vcf
rm /vcf/*tmp*

#QC plots of the filtered vcf dataset without MAF filter applied; aspects of the plots can and should be altered in the R script; please restart pipeline from here if initial thresholds can be optimized
vcftools --gzvcf vcf/"$FILE"_noMAFapplied.recode.vcf.gz --depth --out vcf_analyses/QC_FILTERED/noMAFapplied/"$FILE"_noMAF
vcftools --gzvcf vcf/"$FILE"_noMAFapplied.recode.vcf.gz --site-mean-depth --out vcf_analyses/QC_FILTERED/noMAFapplied/"$FILE"_noMAF
vcftools --gzvcf vcf/"$FILE"_noMAFapplied.recode.vcf.gz --freq2 --max-alleles 2 --out vcf_analyses/QC_FILTERED/noMAFapplied/"$FILE"_noMAF
vcftools --gzvcf vcf/"$FILE"_noMAFapplied.recode.vcf.gz --site-quality --out vcf_analyses/QC_FILTERED/noMAFapplied/"$FILE"_noMAF
vcftools --gzvcf vcf/"$FILE"_noMAFapplied.recode.vcf.gz --missing-indv --out vcf_analyses/QC_FILTERED/noMAFapplied/"$FILE"_noMAF
vcftools --gzvcf vcf/"$FILE"_noMAFapplied.recode.vcf.gz --missing-site --out vcf_analyses/QC_FILTERED/noMAFapplied/"$FILE"_noMAF
vcftools --gzvcf vcf/"$FILE"_noMAFapplied.recode.vcf.gz --het --out vcf_analyses/QC_FILTERED/noMAFapplied/"$FILE"_noMAF
cp R_code/02_vcfQC_noMAF_FILTERED.R R_code/02_vcfQC_"$FILE"_noMAF_FILTERED.R
sed -i "s/FILE/$FILE/g" R_code/02_vcfQC_"$FILE"_noMAF_FILTERED.R
module purge
module load foss/2020b
module load fosscuda/2020b
module load r/4.0.4
R CMD BATCH R_code/02_vcfQC_"$FILE"_noMAF_FILTERED.R
rm R_code/02_vcfQC_"$FILE"_noMAF_FILTERED.R

#QC plots of the filtered vcf dataset with MAF filter applied; aspects of the plots can and should be altered in the R script; please restart pipeline from here if initial thresholds can be optimized
module purge
module load gcc/8.3.0
module load vcftools/0.1.16
vcftools --gzvcf vcf/"$FILE"_MAFapplied.recode.vcf.gz --depth --out vcf_analyses/QC_FILTERED/MAFapplied/"$FILE"_MAF
vcftools --gzvcf vcf/"$FILE"_MAFapplied.recode.vcf.gz --site-mean-depth --out vcf_analyses/QC_FILTERED/MAFapplied/"$FILE"_MAF
vcftools --gzvcf vcf/"$FILE"_MAFapplied.recode.vcf.gz --freq2 --max-alleles 2 --out vcf_analyses/QC_FILTERED/MAFapplied/"$FILE"_MAF
vcftools --gzvcf vcf/"$FILE"_MAFapplied.recode.vcf.gz --site-quality --out vcf_analyses/QC_FILTERED/MAFapplied/"$FILE"_MAF
vcftools --gzvcf vcf/"$FILE"_MAFapplied.recode.vcf.gz --missing-indv --out vcf_analyses/QC_FILTERED/MAFapplied/"$FILE"_MAF
vcftools --gzvcf vcf/"$FILE"_MAFapplied.recode.vcf.gz --missing-site --out vcf_analyses/QC_FILTERED/MAFapplied/"$FILE"_MAF
vcftools --gzvcf vcf/"$FILE"_MAFapplied.recode.vcf.gz --het --out vcf_analyses/QC_FILTERED/MAFapplied/"$FILE"_MAF
cp R_code/02_vcfQC_MAF_FILTERED.R R_code/02_vcfQC_"$FILE"_MAF_FILTERED.R
sed -i "s/FILE/$FILE/g" R_code/02_vcfQC_"$FILE"_MAF_FILTERED.R
module purge
module load foss/2020b
module load fosscuda/2020b
module load r/4.0.4
R CMD BATCH R_code/02_vcfQC_"$FILE"_MAF_FILTERED.R
rm R_code/02_vcfQC_"$FILE"_MAF_FILTERED.R

##You now have a filtered variant call format dataset to use for phylogenetic and population genomics analyses; be sure to check filtered QC plots to determine if you are happy with the dataset##

####################################Admixture analysis to infer population structure at different values of k; always use unlinked (LD pruned) datasets for these analyses#############################################
#Extract unlinked loci for PCA and ADMIXTURE analyses
module purge
module load plink/1.9b_6.21-x86_64
#identify linked loci, --indep-pairwise sets parameters of interest: window size for calculating r2 in kb, step size for moving along reference scaffolds, and r2 limit, above which loci are pruned
plink --vcf vcf/"$FILE"_MAFapplied.recode.vcf.gz --double-id --allow-extra-chr --set-missing-var-ids @:#\$1,\$2 --indep-pairwise "$plink_window_size"'kb' 1 "$plink_r2" --out vcf_analyses/plink/$FILE
#Extract unlinked loci for ADMIXTURE analyses
plink --vcf vcf/"$FILE"_MAFapplied.recode.vcf.gz --double-id --allow-extra-chr --set-missing-var-ids @:#\$1,\$2 --extract vcf_analyses/plink/"$FILE".prune.in --make-bed --pca --out vcf_analyses/plink/$FILE
#convert plink output to vcf for phylogenetic analysis of unlinked SNPs
plink --bfile vcf_analyses/plink/$FILE --allow-extra-chr --recode vcf --out vcf_analyses/plink/"$FILE".LD_pruned
#fixes columns names so admixture doesn't expect human chromosomes
awk '{$1=0;print $0}' vcf_analyses/plink/$FILE.bim > vcf_analyses/plink/"$FILE".bim.tmp
mv vcf_analyses/plink/"$FILE".bim.tmp vcf_analyses/plink/"$FILE".bim
#run admixture k clusters
module load admixture/1.3.0
for i in $(seq $k_low $k_high)
do
admixture --cv vcf_analyses/plink/"$FILE".bed $i > vcf_analyses/admixture/log${i}.out
mv *.P vcf_analyses/admixture
mv *.Q vcf_analyses/admixture
done
#collect cross-validation errors to determine best value of k and prep file for analysis in R; haven't worked out list file yet, should be two column file with sample names and populations/species
awk '/CV/ {print $3,$4}' vcf_analyses/admixture/*out | cut -c 4,7-20 > vcf_analyses/admixture/"$FILE".cv.error
awk '{split($1,name,"."); print $1,name[2]}' vcf_analyses/plink/${FILE}.nosex > vcf_analyses/admixture/"$FILE".list
##plot ancestry proportions at different values of k (admixture results)
module load foss/2019b
module load gcc/8.3.0
module load bcftools/1.9
##Get sample IDs, with some clunky intermediate files
bcftools query -l vcf/"$FILE"_noMAFapplied.recode.vcf.gz > vcf_analyses/admixture/vcf.list.1
sed 's/_/\ /' vcf_analyses/admixture/vcf.list.1 > vcf_analyses/admixture/sample.list.2
sed 's/_n_s.bam//' vcf_analyses/admixture/sample.list.2 > vcf_analyses/admixture/sample.list.3
cut -f 2 -d ' ' vcf_analyses/admixture/sample.list.3 > vcf_analyses/admixture/sample.list.4
sed 's/_/\ /' vcf_analyses/admixture/sample.list.4 > vcf_analyses/admixture/sample.list.5
cut -f 2 -d ' ' vcf_analyses/admixture/sample.list.5 > vcf_analyses/admixture/sample.list.6
paste vcf_analyses/admixture/vcf.list.1 vcf_analyses/admixture/sample.list.6 > vcf_analyses/admixture/ind.list
##The R script provided here was written by Joana Meier, September 2019, modified by Trevor Bringloe March 2022, https://github.com/speciationgenomics/scripts/blob/master/plotADMIXTURE.r
module load foss/2020b
module load fosscuda/2020b
module load r/4.0.4
Rscript --vanilla R_code/03_admixture_plots.R -p $FILE -i vcf_analyses/admixture/ind.list -k "$k_high" -l $pops

##plot PCA of results
##Get list of species/populations to plot
sort vcf_analyses/admixture/sample.list.6 | uniq > vcf_analyses/admixture/pop-species.1.list
#A quirky workaround to get R looping through populations
paste -d, -s vcf_analyses/admixture/pop-species.1.list > vcf_analyses/admixture/pop-species.2.list
cat vcf_analyses/admixture/pop-species.2.list > vcf_analyses/admixture/pop-species.3.list
sed ' 1 s/.*/&,NA/' vcf_analyses/admixture/pop-species.2.list > vcf_analyses/admixture/pop-species.4.list
echo "." | tee -a vcf_analyses/admixture/pop-species.3.list
#PCA plot through R using output through plink
cp R_code/04_PCA_plots.R R_code/04_"$FILE"_PCA_plots.R
sed -i "s/FILE/$FILE/g" R_code/04_"$FILE"_PCA_plots.R
R CMD BATCH R_code/04_"$FILE"_PCA_plots.R
rm R_code/04_"$FILE"_PCA_plots.R
rm vcf_analyses/admixture/sample.list*

########################################################Nuclear population statistics######################################################################################################
#Create lists to specify populations in vcftools
sed 's/_n_s.bam//' vcf_analyses/admixture/vcf.list.1 > vcf_analyses/admixture/vcf.list.2
cat vcf_analyses/admixture/pop-species.1.list | while read line
do
grep $line vcf_analyses/admixture/vcf.list.1 > vcf_analyses/admixture/"$line".indv
done
cat vcf_analyses/admixture/pop-species.1.list | while read line
module purge
do
#Population stats as calculated by vcftools, see manual http://vcftools.sourceforge.net/man_latest.html
module load gcc/8.3.0
module load vcftools/0.1.16
#allele frequencies
vcftools --gzvcf vcf/"$FILE"_noMAFapplied.recode.vcf.gz --keep vcf_analyses/admixture/$line.indv --freq2 --out vcf_analyses/popstats/allele_frequencies/$line
#Fst stat, not automated here, the user will have to specify individuals belonging to populations of interest for FST calculation
#mkdir vcf_analyses/popstats/fst
#vcftools --gzvcf vcf/"$FILE"_noMAFapplied.recode.vcf.gz --keep vcf_analyses/admixture/$line.indv --weir-fst-pop User_specified_list_of_pairwise_individuals --fst-window-size $window --fst-window-step 5000 --out vcf_analyses/popstats/fst/$line
#TajimasD stat; next commands tidy up datasets to merge prior to import into R
vcftools --gzvcf vcf/"$FILE"_noMAFapplied.recode.vcf.gz --keep vcf_analyses/admixture/$line.indv --TajimaD $window --out vcf_analyses/popstats/TajimaD/$line
sed '1s/$/\tloc/;2,$s/$/\t'"$line"/ vcf_analyses/popstats/TajimaD/"$line".Tajima.D > vcf_analyses/popstats/TajimaD/"$line".R1.Tajima.D
tail -n +2 vcf_analyses/popstats/TajimaD/"$line".R1.Tajima.D > vcf_analyses/popstats/TajimaD/"$line".R2.Tajima.D
head -1 vcf_analyses/popstats/TajimaD/"$line".R1.Tajima.D > vcf_analyses/popstats/TajimaD/AA.R2.Tajima.D
#Individual inbreeding coefficients/heterozygosity; next commands tidy up datasets to merge prior to import into R
vcftools --gzvcf vcf/"$FILE"_noMAFapplied.recode.vcf.gz --keep vcf_analyses/admixture/$line.indv --het --out vcf_analyses/popstats/heterozygosity/$line
tail -n +2 vcf_analyses/popstats/heterozygosity/"$line".het > vcf_analyses/popstats/heterozygosity/"$line".R.het
head -1 vcf_analyses/popstats/heterozygosity/"$line".het > vcf_analyses/popstats/heterozygosity/AA.R.het
#optional nucleotide diversity per site; note, a window approach would be better, but vcftools does not report values for genomic windows without SNPs, baising scores upwards in low diversity populations; the user can correct by adding 0s back into a windowed approach, but this must factor in genome size, fragmentation...details that cannot be easily automated
#vcftools --gzvcf vcf/"$FILE"_noMAFapplied.recode.vcf.gz --keep vcf_analyses/admixture/$line.indv --site-pi --out vcf_analyses/popstats/nucleotide_diversity/$line
#sed '1s/$/\tloc/;2,$s/$/\t'"$line"/ vcf_analyses/popstats/nucleotide_diversity/"$line".sites.pi > vcf_analyses/popstats/nucleotide_diversity/"$line".R1.sites.pi
#tail -n +2 vcf_analyses/popstats/nucleotide_diversity/"$line".R1.sites.pi > vcf_analyses/popstats/nucleotide_diversity/"$line".R2.sites.pi
#head -1 vcf_analyses/popstats/nucleotide_diversity/"$line".R1.sites.pi > vcf_analyses/popstats/nucleotide_diversity/AA.R2.sites.pi
#nucleotide diversity using a sliding window approach; note, vcftools does not output a value for windows lacking variant positions, upward biasing estimates (since 0 values are excluded from the table); set a large window to minimize this bias
vcftools --gzvcf vcf/"$FILE"_noMAFapplied.recode.vcf.gz --keep vcf_analyses/admixture/$line.indv --window-pi 50000 --out vcf_analyses/popstats/nucleotide_diversity/$line
sed '1s/$/\tloc/;2,$s/$/\t'"$line"/ vcf_analyses/popstats/nucleotide_diversity/"$line".windowed.pi > vcf_analyses/popstats/nucleotide_diversity/"$line".R1.windowed.pi
tail -n +2 vcf_analyses/popstats/nucleotide_diversity/"$line".R1.windowed.pi > vcf_analyses/popstats/nucleotide_diversity/"$line".R2.windowed.pi
head -1 vcf_analyses/popstats/nucleotide_diversity/"$line".R1.windowed.pi > vcf_analyses/popstats/nucleotide_diversity/AA.R2.windowed.pi
done
#relatedness 
module load gcc/8.3.0
module load vcftools/0.1.16
vcftools --gzvcf vcf/"$FILE"_noMAFapplied.recode.vcf.gz --relatedness --out vcf_analyses/popstats/relatedness/$FILE
vcftools --gzvcf vcf/"$FILE"_noMAFapplied.recode.vcf.gz --relatedness2 --out vcf_analyses/popstats/relatedness/$FILE

#merge datasets for R
cat vcf_analyses/popstats/TajimaD/*.R2.Tajima.D > vcf_analyses/popstats/TajimaD/"$FILE".Tajima.D
cat vcf_analyses/popstats/heterozygosity/*R.het > vcf_analyses/popstats/heterozygosity/"$FILE".het
cat vcf_analyses/popstats/nucleotide_diversity/*.R2.sites.pi > vcf_analyses/popstats/nucleotide_diversity/"$FILE".sites.pi
cat vcf_analyses/popstats/nucleotide_diversity/*.R2.windowed.pi > vcf_analyses/popstats/nucleotide_diversity/"$FILE".windowed.pi
#Prep R script and plot figures
module purge
module load foss/2020b
module load fosscuda/2020b
module load r/4.0.4
cp R_code/06_nuclear_stats.R R_code/06_"$FILE"_nuclear_stats.R
sed -i "s/FILE/$FILE/g" R_code/06_"$FILE"_nuclear_stats.R
sed -i "s/POPS/$pops/g" R_code/06_"$FILE"_nuclear_stats.R
R CMD BATCH R_code/06_"$FILE"_nuclear_stats.R

###########################################################################Alignment and phylogenetic analysis of genomes#########################################################################
##MAUVE does not appear to perform well. The mapped read genomes are aligned, thus the workflow continues with these data

##A block to compile de novo assembled mitochondrial genomes, if the user chooses to align those; otherwise the workflow will further align genomes based on map reads
#cat sample.list | while read line
#do
#cp NOVOPlasty/mito/$line/*.fasta results/organellar/mito/"$line"_mito.fasta
#sed -i "s/Contig_1/"$line"_mito/g" results/organellar/mito/"$line"_mito.fasta
#done
#cat results/organellar/mito/*.fasta > results/organellar/mito/"$FILE"_mito.fasta

##Use MAUVE alignment algorithm to align genomes
#Dependencies/mauve_snapshot_2015-02-13/linux-x64/progressiveMauve --output=results/organellar/mito/alignments/"$FILE"_mito.xmfa results/organellar/Popgenome/mito/"$FILE"_mito.fasta
#path/to/MAUVE/progressiveMauve --alignment-output-dir=results/organellar/mito/alignments --alignment-output-format=nexus --output=results/organellar/mito/alignments/"$FILE"_mito.xmfa results/organellar/mito/"$FILE"_mito.fasta
#convert xmfa to fasta, optional step, please download perl script from https://github.com/kjolley/seq_scripts
#perl Dependencies/xmfa2fasta.pl --align --file results/organellar/mito/alignments/"$FILE"_mito.xmfa > results/organellar/mito/alignments/"$FILE"_mito.MAUVE.fasta

##compile plastid genomes for alignment; note, the stream uses the mapped (and therefore already aligned) read concensus files rather than the de novo genomes, as the latter is prone to breaking the workflow, either because the genomes don't assemble or they don't align properly properly without user input
##A block to compile de novo assembled mitochondrial genomes, if the user chooses to align those; otherwise the workflow will further align genomes based on map reads
#cat sample.list | while read line
#do
#cp NOVOPlasty/plastid/$line/Option_1*.fasta results/organellar/plastid/"$line"_plastid.fasta
#sed -i "s/Contig_1/"$line"_plastid/g" results/organellar/plastid/"$line"_plastid.fasta
#done
#cat results/organellar/plastid/*.fasta > results/organellar/plastid/"$FILE"_plastid.fasta
##Use MAUVE alignment algorithm to align genomes
#Dependencies/mauve_snapshot_2015-02-13/linux-x64/progressiveMauve --output=results/organellar/plastid/alignments/"$FILE"_plastid.xmfa results/organellar/Popgenome/chloro/"$FILE"_plastid.fasta
#path/to/MAUVE/progressiveMauve --alignment-output-format=nexus --output=results/organellar/plastid/alignments/"$FILE"_plastid.xmfa results/organellar/plastid/"$FILE"_plastid.fasta
#convert xmfa to fasta, optional step, please download perl script from https://github.com/kjolley/seq_scripts
#perl Dependencies/xmfa2fasta.pl --align --file results/organellar/plastid/alignments/"$FILE"_plastid.xmfa > results/organellar/plastid/alignments/"$FILE"_plastid.MAUVE.fasta

#############checkpoint in workflow here, default is to analyze variant positions from read mapping, but if is switches to MAUVE alignments, they should ensure the alignments are correct before proceeding, as misaligned regions will have significant downstream impacts on stats########################
#setup specimen info to import into R; the added "NA" to the popgenome.indv file does not affect population statistics, and is done to ensure all sample IDs are included (R reads the final line space into the final specimen ID, preventing R from including in calculations)
cat vcf_analyses/admixture/pop-species.1.list | while read line
do
grep $line vcf_analyses/admixture/vcf.list.2 > vcf_analyses/admixture/"$line".popgenome.tmp1.indv
paste -d, -s vcf_analyses/admixture/"$line".popgenome.tmp1.indv > vcf_analyses/admixture/"$line".popgenome.tmp2.indv
sed ' 1 s/.*/&,NA/' vcf_analyses/admixture/"$line".popgenome.tmp2.indv > vcf_analyses/admixture/"$line".popgenome.indv
done
##get organellar diversity stats via R package PopGenome; note populations (POP1, POP2, POP3...) are in the order specified in pop-species.list
module load foss/2020b
module load fosscuda/2020b
module load r/4.0.4
cp R_code/05_organellar_stats.R R_code/05_"$FILE"_organellar_stats.R
sed -i "s/FILE/$FILE/g" R_code/05_"$FILE"_organellar_stats.R
sed -i "s/POPS/$pops/g" R_code/05_"$FILE"_organellar_stats.R
R CMD BATCH R_code/05_"$FILE"_organellar_stats.R

module purge
module load gcccore/10.3.0
module load python/3.9.5
#convert vcf to fasta formats for phylogenetic analysis
Dependencies/vcf2phylip.py -i vcf_analyses/plink/"$FILE".LD_pruned.vcf --fasta
Dependencies/vcf2phylip.py -i vcf_analyses/plink/"$FILE".LD_pruned.vcf --nexus
mv vcf_analyses/plink/"$FILE".LD_pruned.vcf vcf
mv vcf_analyses/plink/"$FILE".LD_pruned.min4.fasta vcf
Dependencies/vcf2phylip.py -i vcf/"$FILE"_noMAFapplied.recode.vcf.gz --fasta
Dependencies/vcf2phylip.py -i vcf/"$FILE"_noMAFapplied.recode.vcf.gz --nexus
Dependencies/vcf2phylip.py -i vcf/"$FILE"_MAFapplied.recode.vcf.gz --fasta
Dependencies/vcf2phylip.py -i vcf/"$FILE"_MAFapplied.recode.vcf.gz --nexus
#RAxML analysis of vcf files using GTAGAMMA substitution model, no invariant positions (since only SNPs are being analysed)
module purge
module load foss/2019b
module load raxml/8.2.12-hybrid-avx
raxmlHPC -T $threads -f a -m GTRGAMMA -p 12345 -x 12345 -# 1000 -s vcf/"$FILE".LD_pruned.min4.fasta -n results/raxml/$FILE.LD_pruned
raxmlHPC -T $threads -f a -m GTRGAMMA -p 12345 -x 12345 -# 1000 -s vcf/"$FILE"_noMAFapplied.recode.min4.fasta -n results/raxml/$FILE.noMAFapplied
raxmlHPC -T $threads -f a -m GTRGAMMA -p 12345 -x 12345 -# 1000 -s vcf/"$FILE"_MAFapplied.recode.min4.fasta -n results/raxml/$FILE.MAFapplied
#RAxML analysis using GTAGAMMA substitution model on organellar alignments, with invariant positions factored into substitution model (since full genomes are considered)
raxmlHPC -T $threads -f a -m GTRGAMMAI -p 12345 -x 12345 -# 1000 -s results/organellar/Popgenome/mito/"$FILE"_mito.fasta -n results/raxml/$FILE_mito
raxmlHPC -T $threads -f a -m GTRGAMMAI -p 12345 -x 12345 -# 1000 -s results/organellar/Popgenome/chloro/"$FILE"_plastid.fasta -n results/raxml/$FILE_plastid

#If running Splitstree as part of the workflow, run below commands, then see below instructions. Unfortunately, splitstree output cannot be automated easily, you must run as an interactive session
#create executable files for splitstree
cp Dependencies/to_exe_LD_pruned Dependencies/"$FILE"_to_exe_LD_pruned
sed -i "s/REPLACE/"$FILE"/g" Dependencies/"$FILE"_to_exe_LD_pruned
cp Dependencies/to_exe_MAFapplied Dependencies/"$FILE"_to_exe_MAFapplied
sed -i "s/REPLACE/"$FILE"/g" Dependencies/"$FILE"_to_exe_MAFapplied
cp Dependencies/to_exe_noMAFapplied Dependencies/"$FILE"_to_exe_noMAFapplied
sed -i "s/REPLACE/"$FILE"/g" Dependencies/"$FILE"_to_exe_noMAFapplied
#add taxa block and fix formatting to nexus files
cp Dependencies/TAXA_BLOCK Dependencies/"$FILE"_TAXA_BLOCK
sed -i "s/BEGIN DATA/BEGIN CHARACTERS/g" Dependencies/"$FILE"_TAXA_BLOCK
sed -i -e "s/NTAX=XXX/"$(grep -E -o NTAX=[[:alnum:]][[:alnum:]] "$FILE"_MAFapplied.recode.min4.nexus)"/g" Dependencies/"$FILE"_TAXA_BLOCK
cat sample_list | tr '\n' ' ' > sample_list_taxablock
sed -i -e s/LABELSYYY/"$(cat sample_list_taxablock)"/g Dependencies/"$FILE"_TAXA_BLOCK
sed "#NEXUS" vcf/"$FILE".LD_pruned.min4.nexus
sed "#NEXUS" vcf/"$FILE"_MAFapplied.recode.min4.nexus
sed "#NEXUS" vcf/"$FILE"_noMAFapplied.recode.min4.nexus
cat Dependencies/"$FILE"_TAXA_BLOCK 
cat Dependencies/"$FILE"_TAXA_BLOCK vcf/"$FILE".LD_pruned.min4.nexus > vcf/"$FILE".LD_pruned.min4.splitstree.nexus
cat Dependencies/"$FILE"_TAXA_BLOCK vcf/"$FILE"_MAFapplied.recode.min4.splitstree.nexus
cat Dependencies/"$FILE"_TAXA_BLOCK vcf/"$FILE"_noMAFapplied.recode.min4.splitstree.nexus

##need to run splitstree as an interactive session, see github page for more information on how to do this: https://github.com/tbringloe/WGS-NOVAC)
#There are two options: install and use splitstree interactively on windows or mac, or run interactively through remote server with x11 display variable set using mobaXterm (already x11 enabled, need to specify if using puTTY configuration)
#very little resources are needed, be sure to run with x11 specified. Copy below command into command line interface while logged into HPC server
#sinteractive --x11 --time=0:10:00 --ntasks 1
#run command for splitstree, swap in file path to installed splitstree, a png of network will appear in the results folder
#/file/path/to/splitstree4/SplitsTree -g -c Dependencies/"$FILE"_to_exe_LD_pruned
#/file/path/to/splitstree4/SplitsTree -g -c Dependencies/"$FILE"_to_exe_MAFapplied
#/file/path/to/splitstree4/SplitsTree -g -c Dependencies/"$FILE"_to_exe_noMAFapplied

########################################################################Citation information###################################################################################
#none of this is possible without the hard work of those creating bioinformatic programs, please cite the studies below if you used this workflow and are publishing results

#Please see useful pages here:
#Scripts reporduced and modified with permission by joanam: https://github.com/speciationgenomics/scripts
#Scripts for xmfa to fasta conversion by kjolley: https://github.com/kjolley/seq_scripts
#Very useful guide for filtering and visualizing SNPs (from which some of this workflow is derived): https://speciationgenomics.github.io/filtering_vcfs/
#Mermaid wofkflow diagram: https://mermaid.live/
#Clever one liner for VCF compilation step: https://www.ecseq.com/support/ngs-snippets/how-to-run-time-consuming-data-analysis-processes-in-parallel-on-unix-systems

#Citations for programs:
#Alexander, D. H., Lange, K. (2011). Enhancements to the ADMIXTURE algorithm for individual ancestry estimation. BMC Bioinformatics, 12, 1-6. https://dalexander.github.io/admixture/download.html
#Danacek, P., Auton, A., Goncalo, A., Albers, C. A., Banks, E., DePristo, M. A., Handsaker, R., Lunter, G., Marth, G., Sherry, S. T., McVean, G., Durbin, R. & 1000 Genomes Project Analysis Group. (2011). The variant call format and VCFtools. Bioinformatics 27: 2156-8. http://vcftools.sourceforge.net/
#Danacek, P., Bonfield, J. K., Liddle, J., Marshall, J., Ohan, V., Pollard, M. O., Whitwham. A., Keane, T., McCarthy, S. A., Davies, R. M., & Li, H. (2021). Twelve years of SAMtools and BCFtools. Gigascience, 10, giab008. https://github.com/samtools/bcftools.
#Dierckxsens, N., Mardulyn, P., & Smits, G. (2017). NOVOPlasty: de novo assembly of organelle genomes from whole genome data. Nucleic Acids Research, 45, e18. https://github.com/ndierckx/NOVOPlasty
#Huson, D. H., & Bryant, D. (2006). Application of phylogenetic networks in evolutionary studies. Molecular Biology and Evolution, 23, 254-67. https://software-ab.informatik.uni-tuebingen.de/download/splitstree4/welcome.html
#Langmead, B., & Salzberg, S. L. (2012). Fast gapped-read alignment with Bowtie 2. Nature Methods, 9, 357-9. http://bowtie-bio.sourceforge.net/bowtie2/index.shtml
#Li, H., Handsaker, B., Wysoker, A., Fennell, T., Ruan, J., Homer, N., Marth, G., Abecasis, G., Durbin, R., & 1000 Genome Project Data Processing Subgroup. (2009). The sequence alignment/map format and SAMtools. Bioinformatics 25: 2078-9. http://www.htslib.org/
#Pfeifer, B., Wittelsbürger, U., Ramos-Onsins, S. E., & Lercher, M. J. (2014). PopGenome: An efficient swiss army knife for population genomic analyses in R. Molecular Biology and Evolution, 31, 1929-1936.https://cran.r-project.org/web/packages/PopGenome/index.html
#Purcell, S., Neale, B., Todd-Brown, K., Thomas, L., Ferreira, M. A. R., Bender, D., Maller, J., Sklar, P., de Bakker, P. I. W., Daly, M. J., & Sham, P. C. (2007). PLINK: a tool set for whole-genome association and population based linkage analyses. American Journal of Human Genetics, 81, 559-575. https://www.cog-genomics.org/plink/
#Stamatakis, A. (2014). RAxML version 8: a tool for phylogenetic analysis and post-analysis of large phylogenies. Bioinformatics, 30, 1312-3. https://cme.h-its.org/exelixis/web/software/raxml/