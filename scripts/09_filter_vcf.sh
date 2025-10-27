#!/bin/bash

#SBATCH --job-name=09_filter_vcf
#SBATCH --output=09_filter_vcf.out.log
#SBATCH --error=09_filter_vcf.err.log
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
##SBATCH -A molecolb
##SBATCH -p molecolb
#SBATCH --mem-per-cpu=16GB
#SBATCH --time=07:00:00

source ~/.bash_profile
conda activate snps

helpFunction()
{
   echo ""
   echo "Usage: $0 -o parameterO -g parameterG -r parameterR -q parameterQ -i parameterI -m parameterM -a parameterA -c parameterC -w parameterW -n parameterN"
   echo -e "\t-o the name you want to call your VCF, this should match the name you specified in the previous step"
   echo -e "\t-g the name of the genome which was used to align the data"
   echo -e "\t-r minimum depth needed to retain a SNP site"
   echo -e "\t-q the minimum quality threshold for a SNP to be retained (all SNPs with a lower quality score will be excluded)"
   echo -e "\t-i the minimum number of individuals typed to retain a SNP"
   echo -e "\t-m the MAF"
   echo -e "\t-a exclude sites where the average genotype depth is below this threshold"
   echo -e "\t-c the maximum correlation allowed between neighbouring SNPs (together with -w is used to avoid LD)"
   echo -e "\t-w the windown within which to consider SNPs as neighbours (together with -c is used to avoid LD)"
   echo -e "\t-n number of SNPs to retain for primer design"
   exit 1 # Exit script after printing help
}

while getopts "o:r:q:i:m:a:n:g:c:w:" opt
do
   case "$opt" in
      o ) parO="$OPTARG" ;;
      g ) parG="$OPTARG" ;;
      r ) parR="$OPTARG" ;;
      q ) parQ="$OPTARG" ;;
      i ) parI="$OPTARG" ;;
      m ) parM="$OPTARG" ;;
      a ) parA="$OPTARG" ;;
      c ) parC="$OPTARG" ;;
      w ) parW="$OPTARG" ;;
      n ) parN="$OPTARG" ;;
      ? ) helpFunction ;; # Print helpFunction in case parameter is non-existent
   esac
done

# Print helpFunction in case parameters are empty
if [ -z "$parO" ] || [ -z "$parG" ] || [ -z "$parR" ] || [ -z "$parQ" ] || [ -z "$parI" ] || [ -z "$parM" ] || [ -z "$parA" ] || [ -z "$parC" ] || [ -z "$parW" ] || [ -z "$parN" ]
then
   echo "Some or all of the parameters are empty";
   helpFunction
fi

# vcf file name $parO
# genome name $parG
# minimum number of reads $parR
# quality threshold $parQ
# number of individuals $parI
# MAF $parM
# remove sites with an average genotype depth higher than $parA
# LD pruning correlation $parC
# LD pruning window $parW
# Number of SNPs to pull out $parN

src=$PWD
mkdir $src/primer_design
cd $src/vcf

# PREP STEP: index vcf
printf "VCF filtering\n\nPreparation: indexing ${parO}.vcf.gz\n"
bcftools index $parO.vcf.gz
printf "Number of SNPs present in initial file, ${parO}.vcf.gz :\n"
zgrep -vc '^#' $parO.vcf.gz

# STEP 1: keep only biallelic SNPs
printf "\nSTEP 1: Retaining only biallelic SNPs.\n"
bcftools view -v snps -m 2 -M 2 $parO.vcf.gz -O z > $parO.bi_snps.vcf.gz
# index file and print output
bcftools index $parO.bi_snps.vcf.gz
printf "Number of SNPs in ${parO}.bi_snps.vcf.gz :\n"
zgrep -vc '^#' $parO.bi_snps.vcf.gz

# STEP 2: remove SNPs closer than 250 bp to chromosome/contig ends
printf "\nSTEP 2: Retaining SNPs at least 250 bp from the start or end of each chromosome or contig.\n"
samtools faidx $src/genome/$parG
# Create a bedfile of positions not within 250 of contig starts/ends
awk 'BEGIN{FS="\t"; OFS=FS} {print $1,"250",$2-250}' $src/genome/$parG.fai | awk '{if ($3>=250) print $0}' > $src/genome/$parG.250flank.bed
# Only include SNPs in bedfile
bcftools view -R $src/genome/$parG.250flank.bed $parO.bi_snps.vcf.gz -O z > $parO.bi_snps.not_ends.vcf.gz 
# index file and print output
bcftools index $parO.bi_snps.not_ends.vcf.gz
printf "Number of SNPs in $parO.bi_snps.not_ends.vcf.gz :\n"
zgrep -vc '^#' $parO.bi_snps.not_ends.vcf.gz

# STEP 3: Filter by overall depth, quality and sample number. Exclude SNP calls informed by less than parR reads, those below a quality threshold of parQ, and those SNPs genotyped for less than parI individuals
printf "\nSTEP 3: Retaining SNPs with a minimum genotype depth of ${parR} reads, a quality score of at least ${parQ}, and genotyped in at least ${parI} indivdiauls.\n"
bcftools filter -S . -e "FMT/DP<$parR" $parO.bi_snps.not_ends.vcf.gz | bcftools view -e "QUAL<$parQ || AN/2<$parI" -O z > $parO.bi_snps.not_ends.NOGTDP$parR.Q$parQ.SAMP$parI.vcf.gz
# index file and print output
bcftools index $parO.bi_snps.not_ends.NOGTDP$parR.Q$parQ.SAMP$parI.vcf.gz
printf "Number of SNPs in ${parO}.bi_snps.not_ends.NOGTDP${parR}.Q${parQ}.SAMP${parI}.vcf.gz :\n"
zgrep -vc '^#' $parO.bi_snps.not_ends.NOGTDP$parR.Q$parQ.SAMP$parI.vcf.gz

# STEP 4: Filter by minimum allele frequency. Exclude SNPs with a MAF lower than parM
printf "\nSTEP 4: Retaining SNPs with a minimum allele frequency (MAF) of ${parM}.\n"
bcftools view -e "MAF<$parM" -O z $parO.bi_snps.not_ends.NOGTDP$parR.Q$parQ.SAMP$parI.vcf.gz > $parO.bi_snps.not_ends.NOGTDP$parR.Q$parQ.SAMP$parI.MAF$parM.vcf.gz
# index file and print output
bcftools index $parO.bi_snps.not_ends.NOGTDP$parR.Q$parQ.SAMP$parI.MAF$parM.vcf.gz
printf "Number of SNPs in ${parO}.bi_snps.not_ends.NOGTDP${parR}.Q${parQ}.SAMP${parI}.MAF${parM}.vcf.gz :\n"
zgrep -vc '^#' $parO.bi_snps.not_ends.NOGTDP$parR.Q$parQ.SAMP$parI.MAF$parM.vcf.gz

# STEP 5: Filter by average genotype depth
printf "\nSTEP 5: Considering all sites, the average genotype read depth is :\n"
# calculate the average depth (as an integer) across individuals across all sites. Save it as an environment variable and print it
AVDEPTH="$(bcftools query -f '%CHROM\t%POS\t%AN\t%DP\n' $parO.bi_snps.not_ends.NOGTDP$parR.Q$parQ.SAMP$parI.MAF$parM.vcf.gz | awk '{if ($3 > 0) sum_depth += 2 * ($4 / $3); count++} END {if (count > 0) print int(sum_depth/count)}')"
echo ${AVDEPTH}
AvFilt=$((parA * AVDEPTH))
# exclude SNPs if they have an average depth (across all individuals) that is greater than parA times the average across all sites
printf "Retaining SNPs if the average genotype depth for a given site is less than $AvFilt (calculated as $parA times $AVDEPTH).\n"
bcftools view -e "AVG(FMT/DP) > $AvFilt" -O z $parO.bi_snps.not_ends.NOGTDP$parR.Q$parQ.SAMP$parI.MAF$parM.vcf.gz > $parO.bi_snps.not_ends.NOGTDP$parR.Q$parQ.SAMP$parI.MAF$parM.AVMDP$AvFilt.vcf.gz
# index file and print output
bcftools index $parO.bi_snps.not_ends.NOGTDP$parR.Q$parQ.SAMP$parI.MAF$parM.AVMDP$AvFilt.vcf.gz
printf "Number of SNPs in ${parO}.bi_snps.not_ends.NOGTDP${parR}.Q${parQ}.SAMP${parI}.MAF${parM}.AVMDP${AvFilt}.vcf.gz :\n"
zgrep -vc '^#' $parO.bi_snps.not_ends.NOGTDP$parR.Q$parQ.SAMP$parI.MAF$parM.AVMDP$AvFilt.vcf.gz

# STEP 6: Filter by diversity of allelic combinations
printf "\nSTEP 6: Ensuring diversity of allelic combinations. Retaining SNPs for which there is at least one hom for ref, one hom for alt and one het.\n"
# find only SNPs including at least one hom for ref, one hom for alt and one het
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' $parO.bi_snps.not_ends.NOGTDP$parR.Q$parQ.SAMP$parI.MAF$parM.AVMDP$AvFilt.vcf.gz | less -S \
| awk '/0\/0/ && /1\/1/ && /0\/1/ && /0\/1/' > homref_homalt_het.txt
# find only sites including at least one hom for ref, one hom for alt and one het - just chrom and site info
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' $parO.bi_snps.not_ends.NOGTDP$parR.Q$parQ.SAMP$parI.MAF$parM.AVMDP$AvFilt.vcf.gz | less -S \
| awk '/0\/0/ && /1\/1/ && /0\/1/ && /0\/1/' |awk '{ printf "%5s\t%s\n", $1, $2 }' > sites_homref_homalt_het.txt 
# filter to keep sites based on list created above
bcftools view -O z -R sites_homref_homalt_het.txt $parO.bi_snps.not_ends.NOGTDP$parR.Q$parQ.SAMP$parI.MAF$parM.AVMDP$AvFilt.vcf.gz > $parO.bi_snps.not_ends.NOGTDP$parR.Q$parQ.SAMP$parI.MAF$parM.AVMDP$AvFilt.homref.alt.het.vcf.gz
# index file and print output
bcftools index $parO.bi_snps.not_ends.NOGTDP$parR.Q$parQ.SAMP$parI.MAF$parM.AVMDP$AvFilt.homref.alt.het.vcf.gz
printf "Number of SNPs in ${parO}.bi_snps.not_ends.NOGTDP${parR}.Q${parQ}.SAMP${parI}.MAF${parM}.AVMDP${AvFilt}.homref.alt.het.vcf.gz :\n"
zgrep -vc '^#' $parO.bi_snps.not_ends.NOGTDP$parR.Q$parQ.SAMP$parI.MAF$parM.AVMDP$AvFilt.homref.alt.het.vcf.gz

# STEP 7: prune based on LD. Exlude SNPs with an r2 value larger than parC within a parW window.
printf "\nSTEP 7: Retaining SNPs with a correlation (R2) of less than ${parC} with other SNPs within a ${parW} window:\n"
bcftools +prune -m $parC -w $parW $parO.bi_snps.not_ends.NOGTDP$parR.Q$parQ.SAMP$parI.MAF$parM.AVMDP$AvFilt.homref.alt.het.vcf.gz -Oz -o $parO.bi_snps.not_ends.NOGTDP$parR.Q$parQ.SAMP$parI.MAF$parM.AVMDP$AvFilt.homref.alt.het.cor$parC.win$parW.vcf.gz
# index file and print output
bcftools index $parO.bi_snps.not_ends.NOGTDP$parR.Q$parQ.SAMP$parI.MAF$parM.AVMDP$AvFilt.homref.alt.het.cor$parC.win$parW.vcf.gz
printf "Number of SNPs in ${parO}.bi_snps.not_ends.NOGTDP${parR}.Q${parQ}.SAMP${parI}.MAF${parM}.AVMDP${AvFilt}.homref.alt.het.cor${parC}.win${parW}.vcf.gz :\n"
PRUNE="$(zgrep -vc '^#' $parO.bi_snps.not_ends.NOGTDP$parR.Q$parQ.SAMP$parI.MAF$parM.AVMDP$AvFilt.homref.alt.het.cor$parC.win$parW.vcf.gz)"
echo ${PRUNE}

# STEP 8: Random sampling. Check that there are sufficient SNPs available and carry out random sampling. Print an informative message about relationship between final pruned number and parN
printf "\nSTEP 8: Selection of a random subset of SNPs to be taken forward for primer design.\n"
if [ $PRUNE -lt $parN ]
then
    printf "ERROR: Following pruning only ${PRUNE} SNPs were detected. This is less than the the number of SNPs specified to be \n retained for primer design (${parN}). You will have to rerun this filtration script with some filters relaxed. Refer to the \noutput of this log file to see where drop out of SNPs has occured.\n"
elif [ $PRUNE -lt $((2 * $parN)) ]
then
    printf "WARNING: Following pruning only ${PRUNE} SNPs were detected. This is less than two times greater than the the number of \n SNPs specified to be retained for primer design(${parN}). Consider rerunning this filtration script with some filters relaxed. \nRefer to the output of this log file to see where drop out of SNPs has occured.\n"
elif [ $PRUNE -ge $((2 * $parN)) ]
then
    printf "Following pruning ${PRUNE} SNPs were detected. ${parN} will be retained for primer design. Proceeding with final step.\n"
else
    printf "ERROR: Unable to determine how many SNPs have been filtered.\n"
fi

# Carry out subsampling if there are enough SNPs
if [ $PRUNE -ge $parN ]
then
# randomly sample X SNP sites and assign to file
    bcftools view -H $parO.bi_snps.not_ends.NOGTDP$parR.Q$parQ.SAMP$parI.MAF$parM.AVMDP$AvFilt.homref.alt.het.cor$parC.win$parW.vcf.gz | shuf -n $parN | awk '{ printf "%5s\t%s\n", $1, $2 }' > $parN.SNPs.txt
# minus 250bp and add 250bp for range to make bed file
    awk -v s=250 '{print $1, $2-s, $2+s}' $parN.SNPs.txt | sed 's/ /\t/g' > $parN.SNPs.bed
# now extract these sequences from genome file
    bedtools getfasta -fi $src/genome/$parG -bed $parN.SNPs.bed > $src/primer_design/$parN.SNPs.fasta
# filter to keep sites based on list created above. Make vcf of $parN.SNPs.txt sites
    bcftools view -O z -R $parN.SNPs.txt $parO.bi_snps.not_ends.NOGTDP$parR.Q$parQ.SAMP$parI.MAF$parM.AVMDP$AvFilt.homref.alt.het.cor$parC.win$parW.vcf.gz > $parO.bi_snps.not_ends.NOGTDP$parR.Q$parQ.SAMP$parI.MAF$parM.AVMDP$AvFilt.homref.alt.het.cor$parC.win$parW.$parN.vcf.gz
# index new file
    bcftools index $parO.bi_snps.not_ends.NOGTDP$parR.Q$parQ.SAMP$parI.MAF$parM.AVMDP$AvFilt.homref.alt.het.cor$parC.win$parW.$parN.vcf.gz
# extract sequence info
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' $parO.bi_snps.not_ends.NOGTDP$parR.Q$parQ.SAMP$parI.MAF$parM.AVMDP$AvFilt.homref.alt.het.cor$parC.win$parW.$parN.vcf.gz > $src/primer_design/$parN.SNPs_info.txt
    printf "A subset of ${parN} SNPs have been selected for primer design. Please refer to the files primer_design/${parN}.SNPs.fasta and primer_design/${parN}.SNPs_info.txt. Filtration and subsampling complete.\n"
elif [ $PRUNE -lt $parN ]
then
    printf "\nERROR: Unable to subsample SNPs for primer design because there were too few left after filtration.\n"
else
    printf "\nERROR: Unable to complete filtration and subsampling. Please inspect above results and file 09_filter_vcf.err.log"
fi
