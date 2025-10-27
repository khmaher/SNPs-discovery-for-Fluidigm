#!/bin/bash

#SBATCH --job-name=03_fastqc
#SBATCH --output=03_fastqc.log
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
##SBATCH -A molecolb
##SBATCH -p molecolb
#SBATCH --mem-per-cpu=12GB
#SBATCH --time=07:00:00

source ~/.bash_profile
conda activate multiqc

helpFunction()
{
   echo ""
   echo "Usage: $0 -f parameterF -r parameterR"
   echo -e "\t-f the file extension for your forward reads"
   echo -e "\t-r the file extension for your reverse reads"
   exit 1 # Exit script after printing help
}

while getopts "f:r:" opt
do
   case "$opt" in
      f ) parameterF="$OPTARG" ;;
      r ) parameterR="$OPTARG" ;;
      ? ) helpFunction ;; # Print helpFunction in case parameter is non-existent
   esac
done

# Print helpFunction in case parameters are empty
if [ -z "$parameterF" ] || [ -z "$parameterR" ]
then
   echo "Some or all of the parameters are empty";
   helpFunction
fi

src=$PWD

mkdir $src/fastqc
mkdir $src/fastqc/fastqc_F

for f in $src/raw_data/*$parameterF;
do fastqc $f -o $src/fastqc/fastqc_F
done

multiqc $src/fastqc/fastqc_F -o $src/fastqc/fastqc_F/multiqc


mkdir $src/fastqc/fastqc_R

for f in $src/raw_data/*$parameterR;
do fastqc $f -o $src/fastqc/fastqc_R
done

multiqc $src/fastqc/fastqc_R -o $src/fastqc/fastqc_R/multiqc

