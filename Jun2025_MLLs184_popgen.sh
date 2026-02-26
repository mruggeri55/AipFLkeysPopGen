# 29 June 2025
# trying with strict MLL threshold (aka 0.022 predicted by farthest neighbor in poppr)
# sftp MLLs184 onto HPC (note removed tech reps but many were called as unique MLLs)

IN=/scratch1/mruggeri/AipFLkeys/bams_all
OUT=/scratch1/mruggeri/AipFLkeys/angsd/NoClones/MLLs184
prefix=MLLs184
# set N to 80% of total (0.8 x 184 = 147.2)
N=147
maxDP=1840

cd $IN
FILTERS="-uniqueOnly 1 -remove_bads 1 -minMapQ 20 -minQ 25 -dosnpstat 1 -doHWE 1 -sb_pval 1e-5 -skipTriallelic 1 -hetbias_pval 1e-5 -snp_pval 1e-5 -minMaf 0.05"
TODO="-doMajorMinor 1 -doMaf 1 -doCounts 1 -makeMatrix 1 -doIBS 1 -doCov 1 -doGeno 8 -doBcf 1 -doPost 1 -doGlf 2"
angsd -b /scratch1/mruggeri/AipFLkeys/lists/$prefix -GL 1 $FILTERS -minInd $N -maxDepth $maxDP $TODO -P 1 -out $OUT/$prefix
# 661 loci

SCRIPTS=/scratch1/mruggeri/AipFLkeys/scripts
# prune linked sites
nano $SCRIPTS/ngsLD_prefix.sh
# set prefix variable to dataset prefix
# set Nind to total samples
cd $OUT
sbatch $SCRIPTS/ngsLD_prefix.sh

# check how many loci after pruning
gunzip -c $prefix'_unlinked.beagle.gz' | wc -l
# 377 loci after pruning

# also change NGSadmix script to include prefix
nano $SCRIPTS/NGSadmix_v2.sh
mkdir NGSadmix
sbatch $SCRIPTS/NGSadmix_v2.sh

# scp all log and qopt files onto computer
