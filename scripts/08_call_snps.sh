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
   echo "Usage: $0 -g parameterG -o parameterO"
   echo -e "\t-g Description of what is parameterG"
   echo -e "\t-o Description of what is parameterO"
   echo -e "\t-o Description of what is parameterA"
   echo -e "\t-o Description of what is parameterB"
   exit 1 # Exit script after printing help
}

while getopts "g:o:a:b:" opt
do
   case "$opt" in
      g ) parameterG="$OPTARG" ;;
      o ) parameterO="$OPTARG" ;;
      a ) parameterA="$OPTARG" ;;
      b ) parameterB="$OPTARG" ;;
      ? ) helpFunction ;; # Print helpFunction in case parameter is non-existent
   esac
done

# Print helpFunction in case parameters are empty
if [ -z "$parameterG" ]|| [ -z "$parameterO" ]|| [ -z "$parameterA" ]|| [ -z "$parameterB" ]
then
   echo "Some or all of the parameters are empty";
   helpFunction
fi

### Make A and B parameters optional? Any other parameters we want to include?


src=$PWD

# make new directories for vcf file and cleaned file
mkdir $src/vcf

# index genome file for bcftools to use
samtools faidx $src/genome/$parameterG

# Make file of list of bam files
ls $src/clean_aligned/*.bam > $src/bamFiles.txt


# run bcftools mpileup
##-Ou: ouput an uncompressed bam file. This is the option to use when piping the output to another command for optimum performance/speed.
##–max-depth 10000: the maximum number of sequences considered per position
##-q 20: filter out alignments with mapping quality <20
##-Q 20: filter out bases with QS < 20
##-P ILLUMINA: use Illumina platform for indels
##-a FORMAT/DP,FORMAT/AD: output depth and allelic depth
##-f specify the genome reference file, which must be faidx-indexed
##-b list of input bam alignment files
# pipe output and call snps
##-m: use the multiallelic caller
##-v: output variants only
##-f GQ: output genotype quality
##-O z: output in compressed VCF format

bcftools mpileup -Ou \
--max-depth 10000 -q $parameterA -Q $parameterB \
-P ILLUMINA -a FORMAT/DP,FORMAT/AD \
-f $src/genome/$parameterG \
-b $src/bamFiles.txt | \
bcftools call -mv -f GQ \
-O z -o $src/vcf/$parameterO.vcf.gz
