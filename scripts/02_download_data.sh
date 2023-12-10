#!/bin/bash

#SBATCH --job-name=01_download_genome
#SBATCH --output=02_download_data.log
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH -A molecolb
#SBATCH -p molecolb
#SBATCH --mem-per-cpu=12GB
#SBATCH --time=48:00:00

source ~/.bash_profile

src=$PWD

helpFunction()
{
   echo ""
   echo "Usage: $0 -f parameterF"
   echo -e "\t-f a file containing a list of SRR names/numbers for the samples you want to download which should be located in your raw_data folder"
   exit 1 # Exit script after printing help
}

while getopts "f:" opt
do
   case "$opt" in
      f ) parameterF="$OPTARG" ;;
      ? ) helpFunction ;; # Print helpFunction in case parameter is non-existent
   esac
done

# Print helpFunction in case parameters are empty
if [ -z "$parameterF" ]
then
   echo "Some or all of the parameters are empty";
   helpFunction
fi

# Begin script in case all parameters are correct

### download 

cd $src/raw_data

# illumina short read WGR data


SRRname=$parameterF

SRR=$(cat $SRRname)

for file in $SRR; do
    echo $file
    fastq-dump $file --split-files --origfmt --gzip
done