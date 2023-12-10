#!/bin/bash

#SBATCH --job-name=06_alignment
#SBATCH --output=06_alignment.log
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH -A molecolb
#SBATCH -p molecolb
#SBATCH --mem-per-cpu=4GB
#SBATCH --time=72:00:00

source ~/.bash_profile

helpFunction()
{
   echo ""
   echo "Usage: $0 -g parameterG"
   echo -e "\t-f the name of your reference genome"
   exit 1 # Exit script after printing help
}

while getopts "g:" opt
do
   case "$opt" in
      g ) parameterG="$OPTARG" ;;
      ? ) helpFunction ;; # Print helpFunction in case parameter is non-existent
   esac
done

# Print helpFunction in case parameters are empty
if [ -z "$parameterG" ]
then
   echo "Some or all of the parameters are empty";
   helpFunction
fi


src=$PWD

# make directories
# only need to run once
#mkdir $src/genome
mkdir $src/aligned

# align paired reads using bwa mem and output as bam file using samtools
for f in $src/trim/*_trimmed_paired_R1.fastq.gz;
do 		FBASE=$(basename $f)
        BASE=${FBASE%_trimmed_paired_R1.fastq.gz}
        bwa mem -t 20 \
		$src/genome/$parameterG \
		$src/trim/${BASE}_trimmed_paired_R1.fastq.gz \
		$src/trim/${BASE}_trimmed_paired_R2.fastq.gz| \
		samtools sort -o $src/aligned/${BASE}_paired.bam
done

# combine single end reads into one file
# align unpaired reads using bwa mem and output as bam file using samtools
# merge paired and unpaired alignments


for f in $src/trim/*_trimmed_unpaired_R1.fastq.gz;
do 		FBASE=$(basename $f)
        BASE=${FBASE%_trimmed_unpaired_R1.fastq.gz}
        zcat $src/trim/${BASE}_trimmed_unpaired_R1.fastq.gz $src/trim/${BASE}_trimmed_unpaired_R1.fastq.gz > $src/trim/${BASE}_trimmed_unpaired_both.fastq.gz
        bwa mem -t 20 \
		$src/genome/$parameterG \
		$src/trim/${BASE}_trimmed_unpaired_both.fastq.gz| \
		samtools sort -o $src/aligned/${BASE}_unpaired.bam
		samtools merge -@ 4 $src/aligned/${BASE}_all.bam $src/aligned/${BASE}_paired.bam $src/aligned/${BASE}_unpaired.bam
done

