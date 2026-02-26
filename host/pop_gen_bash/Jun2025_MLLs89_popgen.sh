#### pop gen with 90 genets ####
# jun 27, 2025
# using high coverage MLLs from strict threshold in poppr (89 bams)

IN=/scratch1/mruggeri/AipFLkeys/bams_all
OUT=/scratch1/mruggeri/AipFLkeys/angsd/NoClones/MLLs89
prefix=MLLs89
N=72
maxDP=890

cd $IN
FILTERS="-uniqueOnly 1 -remove_bads 1 -minMapQ 20 -minQ 25 -dosnpstat 1 -doHWE 1 -sb_pval 1e-5 -hetbias_pval 1e-5 -skipTriallelic 1 -snp_pval 1e-5 -minMaf 0.05"
TODO="-doMajorMinor 1 -doMaf 1 -doCounts 1 -makeMatrix 1 -doIBS 1 -doCov 1 -doGeno 8 -doBcf 1 -doPost 1 -doGlf 2"
angsd -b /scratch1/mruggeri/AipFLkeys/lists/$prefix -GL 1 $FILTERS -minInd $N -maxDepth $maxDP $TODO -P 1 -out $OUT/$prefix
# 814 loci

SCRIPTS=/scratch1/mruggeri/AipFLkeys/scripts
# prune linked sites
nano $SCRIPTS/ngsLD_prefix.sh
# set prefix variable to dataset prefix
# set Nind to total samples
cd $OUT
sbatch $SCRIPTS/ngsLD_prefix.sh

# check how many loci after pruning
gunzip -c MLLs89_unlinked.beagle.gz | wc -l
# 478 loci left

# also changed NGSadmix script to include a prefix
nano $SCRIPTS/NGSadmix_v2.sh
mkdir NGSadmix
# ran NGSadmix_v2.sh as job -- see contents below
# for j in `seq 1 20` ;
# do
#     for i in `seq 1 7` ; 
#     do
#         NGSadmix -likes $prefix'_unlinked.beagle.gz' -minMaf 0.05 -P 20 -K ${i} -o NGSadmix/$prefix'_k'${i}'run'${j}
#     done
# done
sbatch $SCRIPTS/NGSadmix_v2.sh

# scp all log and qopt files onto computer


### EEEMS #### also need unlinked bcf for hierfstat #
# make bed file from linkage pruned loci
# first need to prune bcf file
# prefix_unlinked.pos contains loci to keep
prefix=MLLs89
module load bcftools 

sed 's/:/\t/' $prefix'_unlinked.pos' > $prefix'_unlinked.posBCFtools'
bcftools view $prefix'.bcf' -Oz -o $prefix'.bcf.gz'
bcftools index -t $prefix'.bcf.gz'
bcftools view -R $prefix'_unlinked.posBCFtools' $prefix'.bcf.gz' -Oz -o $prefix'.unlinked.bcf.gz'

bcftools view -H $prefix'.unlinked.bcf.gz' | wc -l
wc -l $prefix'_unlinked.pos'
# sftp for hierfstat
##### STOPPING HERE #####

# make diff file using plink
module load plink2 bcftools
plink2 --vcf AipClean90.unlinked.bcf.gz --allow-extra-chr --make-bed --out AipClean90.unlinked
# also need to fix chromosome codes
awk '{$1="0";print $0}' AipClean90.unlinked.bim > AipClean90.unlinked.bim.tmp
mv AipClean90.unlinked.bim.tmp AipClean90.unlinked.bim
# now run bed2diffs to make diff file
bed2diffs_v1 --bfile AipClean90.unlinked

# also need '.outer' file with coords of polygon to plot -- note did this in google earth then imported
# and '.coord' file with coords of sample locations
# both need same prefix as input file 'AipClean90.unlinked'

# make param file 
# note this need to be full or relative path (cannot store in different directory and give relative path from current directory bc that won't work)
datapath = ../AipClean90.unlinked
mcmcpath = Aip90_d300_s454
nIndiv = 90
nSites = 454
nDemes = 300
diploid = true
numMCMCIter = 5000000
numBurnIter = 2000000
numThinIter = 9999

# run EEMS 
sbatch -J runEEMS runEEMS.sh 

# then plot results in R
module load geos gdal proj udunits r
R

library(rEEMSplots)
setwd('/scratch1/mruggeri/AipFLkeys/angsd/NoClones/AipClean90/EEMS/')

mcmcpath = "Aip90_d400_s454/"
plotpath = "Aip90_d400_s454/Aip90_d400_s454"
projection_none <- "+proj=longlat +datum=WGS84"
projection_mercator <- "+proj=merc +datum=WGS84"

eems.plots(mcmcpath, plotpath, longlat = F, projection.in = projection_none, projection.out = projection_mercator, out.png=FALSE, add.grid=F, add.demes=T, add.outline=T, add.map=T, add.abline = T, add.r.squared = T)

# repeated for 200, 300 and 400 demes





