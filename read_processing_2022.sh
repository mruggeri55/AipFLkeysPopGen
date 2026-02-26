##### read processing for Aiptasia pop gen (2bRAD) proj #####
# each 'sample' is actually  a pool of samples -- used 1/8th reduction scheme
# some slight modifications to OG script
# processing both 2022 runs separately then will concatenate clean files with 2021 run

#expand files - run as job (make sure copy of gzip is in each directory)
#module load gzip
#gzip.sh <- gzip -d *.gz

sbatch -J gzip gzip.sh #only took ~10min
#
#
#
#
##### count number of reads for each samples and quality using fastqc and multiqc #####
#first, count number of reads in each sequence file and record
countreads.pl > numRawReads.txt

#write a little bash script to iterate fastqc over all files
# make new output direct for fastqc
mkdir fastqc_out

# copy fastqc.sh into directory with reads
# modify the output directory to /scratch1/mruggeri/AipFLkeys2022/run1/fastqc_out OR run2

sbatch -J fastqc fastqc.sh #run this in run1 directory

#run multiqc to aggregate fastqc data
#needed to install multiqc locally
#module load python
#pip install multiqc
#run multiqc as job or doesn't work

module load multiqc 
multiqc /scratch1/mruggeri/AipFLkeys2022/run1/fastqc_out # or run 2

####### concatenate reads across lanes (i.e. for each pool of samples) #######
# first make list of IDs
ls -1 *.fastq | awk -F '_' '{print $1 "_" $2}' | sort | uniq > ID

#then, print list of job commands
for i in `cat ./ID`; do echo cat $i\_*.fastq \> $i\.fastq; done
# copy list of commands into conc.sh and run as job
# move all raw files to their own directory

# counted reads by pool
countreads.pl > reads_by_pool.txt
awk -F' ' '{sum += $2} END{print sum}' reads_by_pool.txt # to sum reads then divide by 33 pools to get avg
# ~5.086 million reads / pool on avg for run 1
# ~4.177 million reads / pool on avg for run 2

#### deduplicate reads and separate pools into samples
for i in `cat ./ID`; do echo trim2bRAD_2barcodes_dedup_N2.pl input=$i.fastq site=\"\.{12}CGA\.{6}TGC\.{12}\|\.{12}GCA\.{6}TCG\.{12}\" adaptor=\"AGATC?\" sampleID=100 deduplicate=1 bc=[ATGC]{4} minBCcount=10000; done 
# copy list of commands into job script and run
sbatch -J pool_adap pool_adap.sh
# will separate out samples and rename them into PoolName_AdaptorSeq.tr0 
# example MR_A_AGAC.tr0
# having issues de-multiplexing only 16 samples worked for run1 and 32 for run2 -- it was because the sequencer already trimmed the adaptor sequence
# works fine with un-trimmed reads

# Carly fixed trimming issue by changing trim2bRAD_2barcodes_dedup_N2.pl to not look for adaptor sequence following the internal BC > now trim2bRAD_2barcodes_dedup_N2_NoAdap.pl
# BUT looks like sequencing facility trimmed the adaptors off (which they shouldn't have)

# Stephanie gave us the non-trimmed files -- used rclone to copy onto HPC
# followed https://www.carc.usc.edu/user-information/user-guides/data-management/transferring-files-rclone with some modifications
# when it prompts for root_folder_id enter the team drive ID
# i.e. from https://drive.google.com/drive/folders/1HwfpwERMTl1CLiqjGebgMSLxt1yThJwI
# enter 1HwfpwERMTl1CLiqjGebgMSLxt1yThJwI
# then at the end when it asks if you want to configure team drive put yes
# now this folder is referred to as google-drive (or whatever you named it)
# made run1 directory 
# rclone copy "google-drive:WithoutAdaptorTrimming_Fastq" "run1"
# rclone copy "untrimmed-run:FASTQ_Generation_2022-06-13_17_47_24Z-572800640.zip" "ckenkel_26"
# yayy it worked -- start over from the top

# now rename files to actual pool names
# first need to replace all of the dashes with underscores bc steph fucked up
rename - _ *.tr0
# let's also remove the _S#_ because not in OG sample name
# honestly cannot figure out a good way to do this so did it an easy but silly way
# see noSnum.sh

# open sftp window
sftp mruggeri@discovery@usc.edu
# transfer pool_adapt_to_new_filename.csv to scratch
put /Users/maria/Desktop/Kenkel_lab/Aip_FLkeys/computational/2bRAD/pool_adapt_to_new_filename.csv /scratch1/mruggeri/AipFLkeys2022/run1
dos2unix pool_adapt_to_new_filename.csv # to convert to unix and remove any weird line endings

#now rename all samples using this file
awk -F ',' '{ print "mv " $1 " " $2 }' pool_adapt_to_new_filename.csv
#copied output into a job script rename.sh and ran

# make new ID list to organize samples into their own directories
ls -1 *.tr0 | awk -F '.' '{print $1}'  > ID

organize2.sh ID

#### now quality filter and adaptor filter ####
# copy this into proc.sh script

module load gcc
module load jdk

#quality filter minimum q20 over 100% of the read
cat $NAME'.tr0' | fastq_quality_filter -q 20 -p 100 > $NAME'.trim'

#adaptor filtering
bbduk.sh in=$NAME'.trim' ref=/project/ckenkel_26/RefSeqs/adaptors.fasta k=12 stats=stats.txt out=$NAME'.clean'

batch.sh proc.sh ID proc

#### count number of reads passing each filter using counts2.sh
# change sequence identifier to match your reads
# commented out mapping perc part bc will concatenate clean reads and then map
sbatch -J counts counts2.sh

##### now concatenate runs before mapping
# concatenating across all 3 runs -- old run + new run 1 & 2
# this is ridiculous but I think it will work
for i in `cat ./ID`; do echo cat run1/$i/$i\.clean run2/$i/$i\.clean oldrun/clean/$i\.clean \> $i\.cat ; done
# copied output into job script conc_runs.sh and ran as job

# count up reads and see if these sums match
countreads.pl \.cat > sum_cat_runs.txt

# need to rename a couple of samples not on pool_adap to new filename list
# see pool_adap_to_new_filename_clean.csv for all IDs
mv CM1_A_GTGA.cat LK_sub2_rt1_127.cat
mv CM1_D_TCAG.cat MI_sub2_rt3_D_.cat
mv CM2_F_TCAC.cat TK_sub1_rt3_75.cat
mv CM2_G_GTGA.cat MI_sub3_rt4_E_TR1.cat
mv CM2_H_CATC.cat MI_sub3_rt1_B_.cat
# MR_H_ATAC.cat actually cannot be found in list -- only 3782 reads in list though so maybe failed prep?

# going to leave all tech reps for now and will remove after downstream processing
# intentional tech reps labeled with *_TR1.cat or *_TR2.cat
# unintentional tech reps listed below
# OK_sub2_rt3_A1.tr0 and OK_sub2_rt3_A2.tr0 tech reps
# BP 257A and 257B are 257, added to plate twice so probs tech reps
# SA 632A and 632B are 632 added to plate twice so probs tech reps

# note: re-organize samples into sample specific directories in order to use batch script
# first make new ID list
ls -1 *.cat | awk -F '.' '{print $1}' | sort | uniq > ID
organize2.sh ID

############### MAPPING ###############
#run mapping script -- had to increase memory a lot, mapped very fast though so probs can cut back to an hour
#GENOME_FASTA=/project/ckenkel_26/mruggeri/Aip_sym_catGenome/Aiptasia_sym_catGenome.fasta
#bowtie2 --threads 16 --no-unal --score-min L,16,1 --local -L 16 -x $GENOME_FASTA -U $NAME".cat" -S $NAME"bt2_HostSym.sam"
batch.sh HostSymMap.sh ID HostSymMap

#get overall alignment rates
grep 'overall alignment rate' */HostSymMap*.err | sed 's/\/HostSymMap.*.err:/\'$'\t/g' > map_perc.txt
sed 's/% overall alignment rate//g' map_perc.txt > sample_map_perc.txt
# will merge this with sum_cat_sort.txt later

# filter out samples with less than 25% alignment
awk '$2>=25' sample_map_perc.txt | cut -f 1 -d " " | sort | uniq > goods
wc -l goods 
# 316 goods -- probably filtered out all SA's
# how many above 50% ? 292 so ~92.4% of goods actually above 50% -- cool
awk '$2>=50' sample_map_perc.txt | cut -f 1 -d " " | sort | uniq | wc -l

# let's look at bads and filter them out
awk '$2<=25' sample_map_perc.txt | cut -f 1 -d " " | sort | uniq > bads
# mostly SAs (I don't think these samples were aiptasia -- collected by UG)
# also BP_sub1_rt2_286, BP_sub3_rt1_285, SI_sub1_rt2_39, TK_sub2_rt2_84 bad
for i in `cat ./bads_ID`; do mv $i/ -t bad_samples/; done

# now count reads mapping to syms
# look up long contig in aiptasia file NW_018384101.1 1835033bps
# moved all sam files mapped to host and sym to their own directory HostSymSam
mkdir HostSymSam
mv */*.sam HostSymSam/
zooxType.pl host="NW_018384101.1" > zooxCounts.txt #run this as job in future, taking a minute

# now need to separate into host and sym sam files
# sym genomes labeled chr11, chr12, chr13, and chr14 so reverse grep host
awk '{print $1}' goods > ID_clean
for i in `cat ./ID`; do echo 'egrep -v' '"chr11|chr12|chr13|chr14"' $i'bt2_HostSym.sam' \> $i\_host.sam ; done
# put output into job script hostSam.sh and run

# move host files into their own directories
ls *.sam > ID
sed 's/\.sam//g' ID > ID_clean
organize2.sh ID_clean

#### bam.sh -- convert to bam files ###
module load gcc
module load intel
module load jdk
module load picard
module load samtools

# enter your job specific code below this line
#GENOME_FASTA=/project/ckenkel_26/mruggeri/Aip_sym_catGenome/Aiptasia_sym_catGenome.fasta
#samtools view -bt $GENOME_FASTA $NAME'.sam' > unsorted.bam && samtools sort -o sorted.bam unsorted.bam && java -Xmx5g -jar $PICARD AddOrReplaceReadGroups INPUT=sorted.bam OUTPUT=$NAME'.bam' RGID=group1 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=$NAME && samtools index $NAME'.bam'

batch.sh bam.sh ID_clean bam
