#!/bin/bash

#SBATCH --job-name=04_trimmomatic
#SBATCH --output=04_trimmomatic.log
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH -A molecolb
#SBATCH -p molecolb
#SBATCH --mem-per-cpu=16GB
#SBATCH --time=72:00:00

source ~/.bash_profile
conda activate multiqc

helpFunction()
{
   echo ""
   echo "Usage: $0 -f parameterF -r parameterR -k parameterK -s parameterS -l parameterL -t parameterT -c parameterC -h parameterH -m parameterM"
   echo -e "\t-f the file extension for your forward reads"
   echo -e "\t-r the file extension for your reverse reads"
   echo -e "\t-k parameters for ILLUMINACLIP"
   echo -e "\t-s parameters for SLIDINGWINDOW"
   echo -e "\t-l parameters for LEADING"
   echo -e "\t-t parameters for TRAILING"
   echo -e "\t-c parameters for CROP"
   echo -e "\t-h parameters for HEADCROP"
   echo -e "\t-m parameters for MINLEN"
   exit 1 # Exit script after printing help
}

# ILLUMINACLIP: Cut adapter and other illumina-specific sequences from the read.
#SLIDINGWINDOW: Perform a sliding window trimming, cutting once the average quality within the window falls below a threshold.
#LEADING: Cut bases off the start of a read, if below a threshold quality
#TRAILING: Cut bases off the end of a read, if below a threshold quality
#CROP: Cut the read to a specified length
#HEADCROP: Cut the specified number of bases from the start of the read
#MINLEN: Drop the read if it is below a specified length


while getopts "f:r:k:s:l:t:c:h:m:" opt
do
   case "$opt" in
      f ) parameterF="$OPTARG" ;;
      r ) parameterR="$OPTARG" ;;
      k ) parameterK="$OPTARG" ;;
      s ) parameterS="$OPTARG" ;;
      l ) parameterL="$OPTARG" ;;
      t ) parameterT="$OPTARG" ;;
      c ) parameterC="$OPTARG" ;;
      h ) parameterH="$OPTARG" ;;
      m ) parameterM="$OPTARG" ;;
      #? ) helpFunction ;; # Print helpFunction in case parameter is non-existent
   esac
done

# Print helpFunction in case parameters are empty
if [ -z "$parameterF" ] || [ -z "$parameterR" ]
then
   echo "Some or all of the parameters are empty";
   helpFunction
fi

# run command

src=$PWD

mkdir $src/trim


for f in $src/raw_data/*$parameterF;
do FBASE=$(basename $f)
	BASE=${FBASE%$parameterF}
	trimmomatic PE -threads 4 -phred33 \
	$src/raw_data/${BASE}$parameterF \
	$src/raw_data/${BASE}$parameterR \
	$src/trim/${BASE}_trimmed_paired_R1.fastq.gz \
	$src/trim/${BASE}_trimmed_unpaired_R1.fastq.gz \
	$src/trim/${BASE}_trimmed_paired_R2.fastq.gz \
	$src/trim/${BASE}_trimmed_unpaired_R2.fastq.gz \
	$parameterK $parameterS $parameterL $arameterT $parameterC $parameterH $parameterM
done

#ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:12 SLIDINGWINDOW:4:30 MINLEN:80
	#$parameterK $parameterS $parameterL $arameterT $parameterC $parameterH $parameterM
