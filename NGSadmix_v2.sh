#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBACTH --mem=4GB
#SBATCH --time=01:00:00
#SBATCH --error=NGSadmix/%x_%j.err
#SBATCH --output=NGSadmix/%x_%j.out
#SBATCH --mail-type=END
#SBATCH --mail-user=mruggeri@usc.edu

prefix=MLLs89

for j in `seq 1 20` ;
do
    for i in `seq 1 7` ; 
    do
        NGSadmix -likes $prefix'_unlinked.beagle.gz' -minMaf 0.05 -P 20 -K ${i} -o NGSadmix/$prefix'_k'${i}'run'${j}
    done
done
