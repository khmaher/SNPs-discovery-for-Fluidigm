#!/bin/bash

#SBATCH --job-name=09_filter_vcf
#SBATCH --output=09_filter_vcf.log
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH -A molecolb
#SBATCH -p molecolb
#SBATCH --mem-per-cpu=16GB
#SBATCH --time=07:00:00

source ~/.bash_profile

helpFunction()
{
   echo ""
   echo "Usage: $0 -o parameterO -r parameterR -q parameterQ -i parameterI -m parameterM -a parameterA -g parameterN -g parameterG"
   echo -e "\t-g Description of what is parameterO"
   echo -e "\t-o Description of what is parameterR"
   echo -e "\t-o Description of what is parameterQ"
   echo -e "\t-o Description of what is parameterI"
   echo -e "\t-o Description of what is parameterM"
   echo -e "\t-o Description of what is parameterA"
   echo -e "\t-o Description of what is parameterN"
   echo -e "\t-o Description of what is parameterG"
   exit 1 # Exit script after printing help
}

while getopts "o:r:q:i:m:a:n:g:" opt
do
   case "$opt" in
      o ) parameterO="$OPTARG" ;;
      r ) parameterR="$OPTARG" ;;
      q ) parameterQ="$OPTARG" ;;
      i ) parameterI="$OPTARG" ;;
      m ) parameterM="$OPTARG" ;;
      a ) parameterA="$OPTARG" ;;
      n ) parameterN="$OPTARG" ;;
      g ) parameterG="$OPTARG" ;;
      ? ) helpFunction ;; # Print helpFunction in case parameter is non-existent
   esac
done

# Print helpFunction in case parameters are empty
if [ -z "$parameterO" ]|| [ -z "$parameterR" ]|| [ -z "$parameterQ" ]|| [ -z "$parameterI" ]|| [ -z "$parameterM" ]|| [ -z "$parameterA" ]|| [ -z "$parameterN" ] || [ -z "$parameterG" ]
then
   echo "Some or all of the parameters are empty";
   helpFunction
fi


#minimum number of reads $parameterR
#quality threshold $parameterQ
#number of individuals $parameterI
#MAF $parameterM
#remove sites with an average genotype depth higher than X $parameterA
#Number of SNPs to pull out
# genome name $parameterG



src=$PWD

#mkdir $src/primer_design

cd $src/vcf

# index vcf

bcftools index $parameterO.vcf.gz

# keep only biallelic SNPs

bcftools view -v snps -m 2 -M 2 $parameterO.vcf.gz -O z > $parameterO.bi_snps.vcf.gz

# index new file

bcftools index $parameterO.bi_snps.vcf.gz


# filter vcf file

# exclude SNP calls informed by less than X reads
# exclude SNPs below a quality threshold of X
# exclude SNPs genotyped for less than X individuals
bcftools filter -S . -e "FMT/DP<$parameterR" $parameterO.bi_snps.vcf.gz | \
bcftools view -e "QUAL<$parameterQ || AN/2<$parameterI" -O z > \
$parameterO.bi_snps.NOGTDP$parameterR.Q$parameterQ.SAMP$parameterI.vcf.gz

# index new file

bcftools index $parameterO.bi_snps.NOGTDP$parameterR.Q$parameterQ.SAMP$parameterI.vcf.gz 

# filter to remove SNPs with a MAF lower than X

bcftools view -e "MAF<$parameterM" -O z $parameterO.bi_snps.NOGTDP$parameterR.Q$parameterQ.SAMP$parameterI.vcf.gz > $parameterO.bi_snps.NOGTDP$parameterR.Q$parameterQ.SAMP$parameterI.$parameterM.vcf.gz

# filter to remove sites with an average genotype depth higher than X
bcftools view -e "AVG(FMT/DP)>$parameterA" -O z $parameterO.bi_snps.NOGTDP$parameterR.Q$parameterQ.SAMP$parameterI.vcf.gz > $parameterO.bi_snps.NOGTDP$parameterR.Q$parameterQ.SAMP$parameterI.$parameterA.vcf.gz
bcftools index $parameterO.bi_snps.NOGTDP$parameterR.Q$parameterQ.SAMP$parameterI.$parameterA.vcf.gz


# find only SNPs including at least one hom for ref, one hom for alt and one het
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' $parameterO.bi_snps.NOGTDP$parameterR.Q$parameterQ.SAMP$parameterI.$parameterA.vcf.gz | less -S \
| awk '/0\/0/ && /1\/1/ && /0\/1/ && /0\/1/' > homref_homalt_het.txt

# find only sites including at least one hom for ref, one hom for alt and one het - just chrom and site info
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' $parameterO.bi_snps.NOGTDP$parameterR.Q$parameterQ.SAMP$parameterI.$parameterA.vcf.gz | less -S \
| awk '/0\/0/ && /1\/1/ && /0\/1/ && /0\/1/' |awk '{ printf "%5s\t%s\n", $1, $2 }' > sites_homref_homalt_het.txt 

# filter to keep sites based on list created above
bcftools view -O z -R sites_homref_homalt_het.txt $parameterO.bi_snps.NOGTDP$parameterR.Q$parameterQ.SAMP$parameterI.$parameterA.vcf.gz > $parameterO.bi_snps.NOGTDP$parameterR.Q$parameterQ.SAMP$parameterI.$parameterA.homref.alt.het.vcf.gz
bcftools index $parameterO.bi_snps.NOGTDP$parameterR.Q$parameterQ.SAMP$parameterI.$parameterA.homref.alt.het.vcf.gz


# randomly sample X SNP sites and assign to file
bcftools view -H $parameterO.bi_snps.NOGTDP$parameterR.Q$parameterQ.SAMP$parameterI.$parameterA.homref.alt.het.vcf.gz | shuf -n $parameterN | awk '{ printf "%5s\t%s\n", $1, $2 }' > $parameterN.SNPs.txt

# minus 250bp and add 250bp for range to make bed file
awk -v s=250 '{print $1, $2-s, $2+s}' $parameterN.SNPs.txt | sed 's/ /\t/g' > $parameterN.SNPs.bed

# now extract these sequences from genome file
bedtools getfasta -fi $src/genome/$parameterG -bed $parameterN.SNPs.bed > $src/primer_design/$parameterN.SNPs.fasta

# filter to keep sites based on list created above
bcftools view -O z -R $parameterN.SNPs.txt $parameterO.bi_snps.NOGTDP$parameterR.Q$parameterQ.SAMP$parameterI.$parameterA.homref.alt.het.vcf.gz > $parameterO.bi_snps.NOGTDP$parameterR.Q$parameterQ.SAMP$parameterI.$parameterA.homref.alt.het.$parameterN.vcf.gz
bcftools index $parameterO.bi_snps.NOGTDP$parameterR.Q$parameterQ.SAMP$parameterI.$parameterA.homref.alt.het.$parameterN.vcf.gz

# make vcf of X SNP sites and extract info
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' $parameterO.bi_snps.NOGTDP$parameterR.Q$parameterQ.SAMP$parameterI.$parameterA.homref.alt.het.$parameterN.vcf.gz > $src/primer_design/$parameterN.SNPs_info.txt


