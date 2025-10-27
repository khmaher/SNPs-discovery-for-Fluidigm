#!/bin/bash

#SBATCH --job-name=07_clean_bam
#SBATCH --output=07_clean_bam.log
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH -A molecolb
#SBATCH -p molecolb
#SBATCH --mem-per-cpu=12GB
#SBATCH --time=48:00:00

source ~/.bash_profile
conda activate snps

src=$PWD

# make new directories for flagstat reports and cleaned bam files
mkdir $src/flagstat
mkdir $src/clean_aligned

# run flagstat to check mapping efficiency 
# use samtools to remove all unmapped reads from the bam (using the -F 4 flag option)
# rerun flagstat to check read numbers in clean bam files

for f in $src/aligned/*_all.bam;
do FBASE=$(basename $f)
	BASE=${FBASE%_all.bam}
	samtools flagstat $src/aligned/${BASE}_all.bam > $src/flagstat/${BASE}.flagstat
	samtools view -F 4 -b $src/aligned/${BASE}_all.bam > $src/clean_aligned/${BASE}_all_mapped.bam
	samtools flagstat $src/clean_aligned/${BASE}_all_mapped.bam > $src/flagstat/${BASE}_mapped.flagstat
done

