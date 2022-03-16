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
# Created by the University of Melbourne job script generator for SLURM
# Mon Dec 09 2019 15:50:02 GMT+1100 (Australian Eastern Daylight Time)

# Partition for the job:
#SBATCH --partition=physical

# Multithreaded (SMP) job: must run on one node and the cloud partition
#SBATCH --nodes=1

# The name of the job:
#SBATCH --job-name="Variant_dataset_build_test"

# The project ID which this job should run under:
#SBATCH --account="punim0608"

# Maximum number of tasks/CPU cores used by the job:
#SBATCH --nodes=1
#SBATCH --ntasks=1
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
#SBATCH --mail-user=trevor.bringloe@unimelb.edu.au

# The maximum running time of the job in days-hours:mins:sec
#SBATCH --time=0-1:59:59

# check that the script is launched with sbatch
if [ "x$SLURM_JOB_ID" == "x" ]; then
   echo "You need to submit your job to the queuing system with sbatch"
   exit 1
fi

# Run the job from the directory where it was launched (default)

# The modules to load:
module load gcccore/11.2.0
module load gcc/8.3.0
module load gcc/10.2.0
module load gcc/10.3.0
module load foss/2019b
module load seqtk/1.3
module load perl/5.34.0
module load trimmomatic/0.39-java-11.0.2
module load bowtie2/2.4.2
module load samtools/1.13
module load bcftools/1.12
module load vcftools/0.1.16
module load plink/2.00-alpha2-x86_64_avx2
module load admixture/1.3.0
module load foss/2020b
module load fosscuda/2020b
module load r/4.0.4

################################################User inputs#############################################################
##Number of threads for analysis
threads=1 #be sure to specify this in the slurm script
EXAMPLE_SAMPLE=AT002 #put an example sample name here for compilation step (uses this sample to extract contig names)
FILE=AE_test #an informative prefix to carry through analyses

##User inputs following file settings
reference_genome=alaria_v3.fasta
repeat_regions=taxa_hardmask.bed # optional file to specify repeat regions to exclude from VCF file
mito_seed=Alaria_COI.fasta
plastid_seed=Alaria_rbcL.fasta

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
map_param=0.3

##User inputs following filtering settings for compiling VCF file
min_cov=15
max_cov=100
allelic_balance_low=0.2
allelic_balance_high=5
minor_allele_frequency=0.02
min_Q=30 #follows phred scale, 30=1/1000 chance of SNP calling error
max_missing=0.9 # value is proportion of present data needed to keep a SNP site, so 0.9=10% missingness

##User specifies following parameters for removal of linked SNPs using PLINK
plink_r2=0.15
plink_window_size=25 #in kb

##User specifies values of k for ADMIXTURE analysis
k_low=2
k_high=3
pops=Arctic,Atlantic,Greenland1,Greenland2,Faroe_Islands,Norway,Ireland #here, list populations in order they should appear in admixture plot, seperated by a comma
################################################PIPELINE####################################################################
##Step 1 is to assemble the organellar genomes using NOVOPlasty
##Mito
cat sample.list | while read line
do
#setup file paths for analysis
mkdir NOVOPlasty/mito/$line
#Extract subreads for assembly, since we don't typically need all the reads to assemble the organelles; seqtk outputs uncompressed reads, which we need anyways for NOVOPlasty; -s flag specifies a random seed, which we need to grab same forward and reverse reads
seqtk sample -s100 read_files/"$line"_1.fq.gz "$Number_reads_to_extract" > read_files/"$line"_sub_1.fq
seqtk sample -s100 read_files/"$line"_2.fq.gz "$Number_reads_to_extract" > read_files/"$line"_sub_2.fq
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
cp NOVOPlasty/NOVOPlasty4.2.pl NOVOPlasty/mito/$line
cp NOVOPlasty/README.md NOVOPlasty/mito/$line
#Run NOVOPlasty; be sure to check seed region, as the seed sequence can sometimes get incorporated into the final genome, presumably when coverage is low
NOVOPlasty/mito/$line/NOVOPlasty4.2.pl -c NOVOPlasty/mito/$line/config_mito_"$line".txt
mv *"$line"* NOVOPlasty/mito/$line
done

##Plastid## Can remove this portion of the script if not a photosynthesic organism
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
cp NOVOPlasty/NOVOPlasty4.2.pl NOVOPlasty/plastid/$line
cp NOVOPlasty/README.md NOVOPlasty/plastid/$line
#Run NOVOPlasty
NOVOPlasty/mito/$line/NOVOPlasty4.2.pl -c NOVOPlasty/plastid/$line/config_plastid_"$line".txt
#Output is in working directory, so move to relevant sample folder
mv *"$line"* NOVOPlasty/plastid/$line
done
#Remove subsampled reads
rm read_files/*sub_1.fq
rm read_files/*sub_2.fq

#Index reference genome for mapping
module load gcc/10.2.0
module load bowtie2/2.4.2
bowtie2-build reference_genome/$reference_genome IDX/reference_genome

#Loop will perform QC, trim, and map
cat sample.list | while read line
do
#Perform QC on raw files
mkdir read_files/raw/$line
fastqc -thread $threads read_files/"$line"_1.fq.gz read_files/"$line"_2.fq.gz -o read_files/raw
#trim raw read files based on specified criteria; trimmed files are deposited into tmp/trimming before moving into tmp as an annoying workaround, trimmomatic otherwise deletes and overwrites files in the loop
java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE -threads $threads read_files/"$line"_1.fq.gz read_files/"$line"_2.fq.gz tmp/trimming/"$line"_1_tp.fq.gz tmp/trimming/"$line"_1_tup.fq.gz tmp/trimming/"$line"_2_tp.fq.gz tmp/trimming/"$line"_2_tup.fq.gz TRAILING:$Trailing HEADCROP:$Headcrop AVGQUAL:$Average_quality MINLEN:$Min_length
mv tmp/trimming/"$line"_1_tp.fq.gz tmp/"$line"_trimmed_1.fq.gz
mv tmp/trimming/"$line"_2_tp.fq.gz tmp/"$line"_trimmed_2.fq.gz
mkdir read_files/trimmed/$line
fastqc -thread $threads tmp/"$line"_trimmed_1.fq.gz tmp/"$line"_trimmed_2.fq.gz -o read_files/trimmed
done

#Summarize QC, be sure to evaluate
module load foss/2019b
module load multiqc/1.9-python-3.7.4
multiqc read_files/raw -o read_files/raw
multiqc read_files/trimmed -o read_files/trimmed

#Map reads to reference genome, then convert to sorted bam file
module load gcc/10.2.0
module load bowtie2/2.4.2
cat sample.list | while read line
do
bowtie2 --score-min L,0,-"$map_param" --no-unal -p $threads -x IDX/reference_genome -1 tmp/"$line"_trimmed_1.fq.gz -2 tmp/"$line"_trimmed_2.fq.gz -S tmp/"$line"_n.sam
done

#Sort, convert, and index Sequence Alignment Map to binary format for compiling step; steps are broken into multiple loops because environment dependencies are not compatible
module load gcc/8.3.0
module load gcc/10.2.0
module load samtools/1.12
module load bcftools/1.12
cat sample.list | while read line
do
samtools view -S -@ "$threads" -b tmp/"$line"_n.sam > tmp/"$line"_n.bam
#File paths get incorporated into vcf file downstream, so sorted bams are deposited into working directory and moved later
samtools sort -@ "$threads" -o "$line"_n_s.bam tmp/"$line"_n.bam
#rm tmp/$line*
samtools index -@ "$threads" "$line"_n_s.bam
done

module load gcc/8.3.0
module load gcc/10.2.0
module load samtools/1.12
module load bcftools/1.12

#Compile BAM files into VCF file; this step takes a long time depending on the size of the dataset, the process is done in parallel here by dividing the task up by contig
ls *_n_s.bam > sample_list
#mpileup and call on each contig, specified parallel pileups by creating xargs list of each contig. Use samtools view to obtain a list of all contig names, and input each as a xargs command, with -P as number of parallel jobs
samtools view -@ "$threads" -H "$EXAMPLE_SAMPLE"_n_s.bam | grep "\@SQ" | sed 's/^.*SN://g' | cut -f 1 | xargs -I {} -n 1 -P $threads sh -c "bcftools mpileup --threads $threads -Ou -f reference_genome/$reference_genome -a AD,ADF,ADR,DP,DP4 -r {} -b sample_list | bcftools call --threads $threads -f GQ,GP -mv -Oz > tmp/tmp.{}.vcf.gz"
#list temporary vcf files
ls tmp/tmp* > list_tmp
#concat into single vcf
bcftools concat -f list_tmp -Oz > vcf/"$FILE"_raw.vcf.gz
rm tmp/tmp*
rm list_tmp
mv *_n_s.bam sorted_bam
mv *n_s.bam.bai sorted_bam

##Congrats, the computationally intensive part is over, consider running the remaining script with less threads if job hangs or you need to optmize parameters##

module purge
module load gcc/8.3.0
module load vcftools/0.1.16
##QC plots of the raw vcf dataset; aspects of the plots can and should be altered in the R script; please restart pipeline from here if initial thresholds can be optimized
vcftools --gzvcf vcf/"$FILE"_raw.vcf.gz --depth --out vcf_analyses/QC_RAW/$FILE
vcftools --gzvcf vcf/"$FILE"_raw.vcf.gz --site-mean-depth --out vcf_analyses/QC_RAW/$FILE
vcftools --gzvcf vcf/"$FILE"_raw.vcf.gz --freq2 --max-alleles 2 --out vcf_analyses/QC_RAW/$FILE
vcftools --gzvcf vcf/"$FILE"_raw.vcf.gz --site-quality --out vcf_analyses/QC_RAW/$FILE
vcftools --gzvcf vcf/"$FILE"_raw.vcf.gz --missing-indv --out vcf_analyses/QC_RAW/$FILE
vcftools --gzvcf vcf/"$FILE"_raw.vcf.gz --missing-site --out vcf_analyses/QC_RAW/$FILE
vcftools --gzvcf vcf/"$FILE"_raw.vcf.gz --het --out vcf_analyses/QC_RAW/$FILE
cp R_code/01_vcfQC_RAW.R R_code/01_vcfQC_"$FILE"_RAW.R
sed -i "s/FILE/$FILE/g" R_code/01_vcfQC_"$FILE"_RAW.R
module purge
module load foss/2020b
module load fosscuda/2020b
module load r/4.0.4
R CMD BATCH R_code/01_vcfQC_"$FILE"_RAW.R
rm R_code/01_vcfQC_"$FILE"_RAW.R

##Apply filtering criteria
#set heterozygous alleles with allelic imbalance biased towards reference allele (e.g. if set to 5, then 5 reference to each 1 alternate allele) to missing
module purge
module load foss/2019b
module load gcc/8.3.0
module load bcftools/1.9
bcftools +setGT vcf/"$FILE"_raw.vcf.gz -- -t q -i 'GT="het" & FMT/AD[1-13:0]/FMT/AD[1-13:1]>'"$allelic_balance_high" -n ./. > vcf/"$FILE"_tmp1.vcf
#set heterozygous alleles with allelic imbalance biased towards alternate allele (e.g. if set to <0.2, then 1 reference to each 5 alternate allele) to missing
bcftools +setGT vcf/"$FILE"_tmp1.vcf -- -t q -i 'GT="het" & FMT/AD[1-13:0]/FMT/AD[1-13:1]<'"$allelic_balance_low" -n ./. > vcf/"$FILE"_tmp2.vcf
#set alleles with low read depth to missing
bcftools +setGT vcf/"$FILE"_tmp2.vcf -- -t q -i 'FMT/DP<'"$min_cov" -n ./. > vcf/"$FILE"_tmp3.vcf
##set alleles with high read depth to missing; note raw QC plots can be used to inform parameter settings
bcftools +setGT vcf/"$FILE"_tmp3.vcf -- -t q -i 'FMT/DP>'"$max_cov" -n ./. > vcf/"$FILE"_tmp4.vcf
#optional step to remove samples with low overall coverage, i.e. sequencing effort is insufficient in these individuals and they should be removed before downstream analysis; replace XXX with samples names in vcf file (will be the sorted bam file name, e.g. XXX_n_s.bam), --remove-indv XXX can be specified multiple times for individual samples; if using this command, remove the hashtag prior to the command line and change $FILE_tmp4.vcf to $FILE_tmp4.1.recode.vcf in the next command step
#vcftools --gzvcf "$FILE"_tmp4.vcf --remove-indv XXX --recode --out "$FILE"_4.1.vcf
#optional step to remove sites from low complexity and repeat regions; bed file can be created from output of repeat analysis when annotating the reference genome; if using this command, remove the hashtag prior to the command line and change $FILE_tmp4.vcf or $FILE_tmp4.1.recode.vcf to $FILE_tmp4.1.recode.vcf or $FILE_tmp4.2.recode.vcf in the next command step
#vcftools --gzvcf "$FILE"_tmp4.1.vcf --exclude-bed vcf_analyses/$repeat_regions --recode --out "$FILE"_tmp4.2.vcf
##Filter sites for Q values [quality scores] and missing data above specified threshold
module purge
module load gcc/8.3.0
module load vcftools/0.1.16
vcftools --gzvcf vcf/"$FILE"_tmp4.vcf --max-missing "$max_missing" --minQ "$min_Q" --recode --out vcf/"$FILE"_noMAFapplied
##Use vcftools to remove sites with minor allele frequency of less than 0.02; note, some pop statistics will be performed on dataset not filtered for MAF
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
vcftools --gzvcf vcf/"$FILE"_MAFapplied.recode.vcf.gz --missing-indv --out vcf_analyses/QC_FILTERED/MAFapplied/"$FILE"_MA
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

##Admixture analysis to infer population structure at different values of k; always use unlinked (LD pruned) datasets for these analyses
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
##Get sample IDs
bcftools query -l vcf/"$FILE"_raw.vcf.gz > vcf_analyses/admixture/sample.list.1
sed 's/_/\ /' vcf_analyses/admixture/sample.list.1 > vcf_analyses/admixture/sample.list.2
sed 's/_n_s.bam//' vcf_analyses/admixture/sample.list.2 > vcf_analyses/admixture/sample.list.3
cut -f 2 -d ' ' vcf_analyses/admixture/sample.list.3 > vcf_analyses/admixture/sample.list.4
paste vcf_analyses/admixture/sample.list.1 vcf_analyses/admixture/sample.list.4 > vcf_analyses/admixture/ind.list
rm vcf_analyses/admixture/sample.list*
##The R script provided here was written by Joana Meier, September 2019, modified by Trevor Bringloe March 2022, https://github.com/speciationgenomics/scripts/blob/master/plotADMIXTURE.r
module load foss/2020b
module load fosscuda/2020b
module load r/4.0.4
Rscript --vanilla R_code/03_admixture_plots.R -p $FILE -i vcf_analyses/admixture/ind.list -k "$k_high" -l $pops

##plot PCA of results
##Get list of species/populations to plot
sort vcf_analyses/admixture/sample.list.4 | uniq > vcf_analyses/admixture/sample.list.5
paste -d, -s vcf_analyses/admixture/sample.list.5 > vcf_analyses/admixture/pop-species.list
#PCA plot through R using output through plink
cp R_code/04_PCA_plots.R R_code/04_"$FILE"_PCA_plots.R
sed -i "s/FILE/$FILE/g" R_code/04_"$FILE"_PCA_plots.R
R CMD BATCH R_code/04_"$FILE"_PCA_plots.R
rm R_code/04_"$FILE"_PCA_plots.R
