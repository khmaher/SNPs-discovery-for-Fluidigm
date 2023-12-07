#!/bin/bash

#SBATCH --job-name=05_fastqc2
#SBATCH --output=05_fastqc2.log
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH -A molecolb
#SBATCH -p molecolb
#SBATCH --mem-per-cpu=12GB
#SBATCH --time=07:00:00

source ~/.bash_profile
conda activate multiqc

src=$PWD

mkdir $src/fastqc2
mkdir $src/fastqc2/fastqc_F

for f in $src/trim/*_trimmed_paired_R1.fastq.gz;
do fastqc $f -o $src/fastqc2/fastqc_F
done

multiqc $src/fastqc2/fastqc_F -o $src/fastqc2/fastqc_F/multiqc


mkdir $src/fastqc2/fastqc_R

for f in $src/trim/*_trimmed_paired_R2.fastq.gz;
do fastqc $f -o $src/fastqc2/fastqc_R
done

multiqc $src/fastqc2/fastqc_R -o $src/fastqc2/fastqc_R/multiqc

