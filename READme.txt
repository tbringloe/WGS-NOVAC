########################WGS-NOVAC workflow documentation#########################3

Start by consulting https://github.com/tbringloe/WGS-NOVAC, see references at the end of the github page for papers and other github pages relevant to development of the workflow

In order to launch the workflow, start by put all the downloaded github files into a single directory in your workspace
run the following commands from the terminal (preferably in a new folder for your workspace), this will setup all the folders and move files used for the analysis
It is ok to create folders and move files around on shared HPC servers, but never submit workflow commands through the terminal, I know people who have crashed the whole server this way :///// Needless to say angry emails ensued
mkdir Dependencies
mkdir IDX
mkdir NOVOPlasty
mkdir NOVOPlasty/mito
mkdir NOVOPlasty/plastid
mkdir NOVOPlasty/seed_files
mkdir R_code
mkdir read_files
mkdir reference_genome
mkdir results
mkdir results/organellar
mkdir results/organellar/mito
mkdir results/organellar/mito/alignments
mkdir results/organellar/mito/mapping
mkdir results/organellar/plastid
mkdir results/organellar/plastid/alignments
mkdir results/organellar/plastid/mapping
mkdir results/organellar/Popgenome
mkdir results/organellar/Popgenome/chloro
mkdir results/organellar/Popgenome/mito
mkdir results/admixture
mkdir results/nuclear_stats
mkdir results/pca
mkdir results/raxml
mkdir results/read_QC
mkdir results/read_QC/raw
mkdir results/read_QC/trimmed
mkdir results/splitstree
mkdir results/vcf_QC
mkdir results/vcf_QC/raw
mkdir results/vcf_QC/MAFapplied
mkdir results/vcf_QC/noMAFapplied
mkdir sorted_bam
mkdir tmp
mkdir tmp/trimming
mkdir vcf
mkdir vcf_analyses
mkdir vcf_analyses/admixture
mkdir vcf_analyses/plink
mkdir vcf_analyses/popstats
mkdir vcf_analyses/popstats/allele_frequencies
mkdir vcf_analyses/popstats/heterozygosity
mkdir vcf_analyses/popstats/nucleotide_diversity
mkdir vcf_analyses/popstats/relatedness
mkdir vcf_analyses/popstats/TajimaD
mkdir vcf_analyses/QC_FILTERED
mkdir vcf_analyses/QC_FILTERED/MAFapplied
mkdir vcf_analyses/QC_FILTERED/noMAFapplied
mkdir vcf_analyses/QC_RAW
mv TAXA_BLOCK Dependencies
mv to_exe_LD_pruned Dependencies
mv to_exe_MAFapplied Dependencies
mv to_exe_noMAFapplied Dependencies
mv vcf2phylip.py Dependencies
mv config_mito.txt NOVOPlasty
mv config_plastid.txt NOVOPlasty
mv 01_vcfQC_RAW.R R_code
mv 02_vcfQC_MAF_FILTERED.R R_code
mv 02_vcfQC_noMAF_FILTERED.R R_code
mv 03_admixture_plots.R R_code
mv 04_PCA_plots.R R_code
mv 05_organellar_stats.R R_code
mv 06_nuclear_stats.R R_code

If using MAUVE, NOVOPlasty, or Splitstree, these likely need to be downloadeded and installed, with the former two dropped into respective folder locations. Depending on the date and version, you will have to go into the workflow script and modify the program filepath to ensure it is called.
The user will also need to check what modules are available if using a shared HPC environment, or otherwise download, install, and change pathways for executing all programs in the script (sorry, dependency management remains a major barrier in need of improvement)
To check if a module is available on your HPC cluster, type 'module spider <program of interest>' into the terminal. If available, options and versions numbers will be listed.
If available, also check if other modules are needed to load a program by typing 'module spider <program of interest with version number> into the terminal. This should list any other modules to list before loading the program.
Alter all the load module sections of the script to properly load dependencies. Note, the workflow has only been tested with the versions listed in the workflow script. As far as I know, the only critical version to use is bcftools1.9 for variant filtering (later versions have changed expression syntax)
Sometimes modules are incompatible with each other on an HPC, hence the 'module purge' command in the workflow to reset the module environment when switching to a new program
Also load R and install libraries needed to execute the R scripts. These include: Here, tidyverse, rprojroot, ggplot2, optparse, PopGenome, gplots, and rmarkdown. Again, dependency management remains a barrier, I was able to run all of these libraries using R.4.0.4

MAUVE linux download here: https://darlinglab.org/mauve/download.html
Drop uncompressed download into /Dependencies

NOVOPlasty download here: https://github.com/ndierckx/NOVOPlasty
Drop files, including filter_reads.pl, LICENSE, NOVOPlasty4.2.pl (or whatever version you are using), and README.md into /NOVOPlasty
Check permissions of the NOVOPlastyX.X.pl file with 'ls -l NOVOPlastyX.X.pl. If permissions are limited, change using 'chmod a+rwx NOVOPlastyX.X.pl'
Do not replace the config_mito.txt and config_plastid.txt files, these are edited to work with the NOVAC workflow

Splitstree linux download here: https://software-ab.informatik.uni-tuebingen.de/download/splitstree4/welcome.html

Make sure you put reference genomes in the reference_genome folder, and seed files if using NOVOPlasty in NOVOPlasty/seed_files
Make sure you specify user input variables in the NOVAC_master_30iii22.txt file
Make sure you configure the read files to match sample IDs in sample.list, following this format: 01_SampleID_Pop/Species_1/2.fq.gz; make sure specimens are in a logical order, grouped by populations, in an order that is meaninful for figures.
Make sure you provide a list of sample IDs in sample.list. A quick way to do this is to move to the read_files folder then type 'ls *_1.fq.gz > sample.list', then 'sed -i 's/_1.fq.gz//' sample.list'. That will remove suffixes from the list. Use 'nano sample.list' to ensure sample IDs are correct (##_SampleID_Pop/species), then move sample.list back to root working space folder where you will launch the job
Once the above steps are completed run the workflow from the root workspace directory by typing the following command in the terminal: sbatch NOVAC_master_30iii22.slurm

Remember other useful commands for monitoring your job

check_project_usage # provides information on project storage
check_scratch_usage # provides information on scratch storage
squeue -u <username> # provides status of submitted jobs
sshare -l -A <projectID> -u <username> # provides information on usage relative to other users both within and across projects, which informs you on que priority. LevelFS is the most relevant metric, values above 1 mean you have priority, values below 1 indicate less priority. The more resources you use, the lower the value, and thus the lower your priority relative to other jobs. Think of usage in terms of core-hours (number of cores*number of hours a job runs for). The half life on resource usage should be ~1 week, so LevelFS will climb over time. Project is punim0608 for algal lab
ls -l --block-size=M # provides file sizes of current folder in MB
nano <filename> # opens files and allows for modifications (e.g. altering slurm script parameters)
nohup <command> # for jobs that run through the command interface, this will ensure the job will continue if you quit your session
logout # quit session

An Rmarkdown script is provided. In order to run, download the results folder from the workspace, open the script in Rstudio (this has built in functions to knit an html of results) and set working directory to the results folder before running the script.

If you run into troubles or have suggestions for improvements, you can open an issue at https://github.com/tbringloe/WGS-NOVAC