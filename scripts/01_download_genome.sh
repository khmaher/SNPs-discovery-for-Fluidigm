#!/bin/bash

#SBATCH --job-name=01_download_genome
#SBATCH --output=01_download_genome.log
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH -A molecolb
#SBATCH -p molecolb
#SBATCH --mem-per-cpu=12GB
#SBATCH --time=07:00:00

source ~/.bash_profile

src=$PWD

helpFunction()
{
   echo ""
   echo "Usage: $0 -w parameterW -g parameterG"
   echo -e "\t-w Description of what is parameterW"
   echo -e "\t-g Description of what is parameterG"
   exit 1 # Exit script after printing help
}

while getopts "w:g:" opt
do
   case "$opt" in
      w ) parameterW="$OPTARG" ;;
      g ) parameterG="$OPTARG" ;;
      ? ) helpFunction ;; # Print helpFunction in case parameter is non-existent
   esac
done

# Print helpFunction in case parameters are empty
if [ -z "$parameterW" ] || [ -z "$parameterG" ]
then
   echo "Some or all of the parameters are empty";
   helpFunction
fi

# Begin script in case all parameters are correct

cd $src/genome

### download genome
wget $parameterW

### index genome
bwa index $parameterG