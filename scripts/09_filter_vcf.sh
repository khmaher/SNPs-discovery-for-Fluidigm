#!/bin/bash

#Settings for the Sun Grid Engine

# run time for job in hours:mins:sec (max 168:0:0, jobs with h_rt < 8:0:0 have priority)

#$ -l h_rt=07:59:59

#$ -pe openmp 1

# request memory for job (default 2G, max 256G)

#$ -l rmem=24G

#$ -q evolgen.q

#$ -P evolgen

# give the job a name (optional):

#$-N filter_snps



cd /fastdata/bo4kma/monkparakeet/vcf

# index vcf

bcftools index monkparakeet.vcf.gz

# keep only biallelic SNPs

bcftools view -v snps -m 2 -M 2 monkparakeet.vcf.gz -O z > monkparakeet_bi_snps.vcf.gz

# index new file

bcftools index monkparakeet_bi_snps.vcf.gz


# filter vcf file

# exclude SNP calls informed by less than 3 reads
# exclude SNPs below a quality threshold of 20
# exclude SNPs genotyped for less than 4 individuals
bcftools filter -S . -e 'FMT/DP<3' monkparakeet_bi_snps.vcf.gz | \
bcftools view -e 'QUAL<20 || AN/2<4' -O z > \
monkparakeet_bi_snps.NOGTDP3.Q20.SAMP4.vcf.gz

# index new file

bcftools index monkparakeet_bi_snps.NOGTDP3.Q20.SAMP4.vcf.gz 

# filter to remove SNPs with a MAF lower than 0.25

bcftools view -e 'MAF<0.25' -O z monkparakeet_bi_snps.NOGTDP3.Q20.SAMP4.vcf.gz > monkparakeet_bi_snps.NOGTDP3.Q20.SAMP4.MAF0.25.vcf.gz

# keep only SNPs with a MAF=0.5
bcftools view -i 'MAF=0.5' -O z monkparakeet_bi_snps.NOGTDP3.Q20.SAMP4.vcf.gz > monkparakeet_bi_snps.NOGTDP3.Q20.SAMP4.MAF0.5.vcf.gz
bcftools index monkparakeet_bi_snps.NOGTDP3.Q20.SAMP4.MAF0.5.vcf.gz

# find only SNPs including at least one hom for ref, one hom for alt and one het
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' monkparakeet_bi_snps.NOGTDP3.Q20.SAMP4.MAF0.5.vcf.gz | less -S \
| awk '/0\/0/ && /1\/1/ && /0\/1/ && /0\/1/' > homref_homalt_het.txt

# find only sites including at least one hom for ref, one hom for alt and one het - just chrom and site info
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' monkparakeet_bi_snps.NOGTDP3.Q20.SAMP4.MAF0.5.vcf.gz | less -S \
| awk '/0\/0/ && /1\/1/ && /0\/1/ && /0\/1/' |awk '{ printf "%5s\t%s\n", $1, $2 }' > sites_homref_homalt_het.txt 

# filter to keep sites based on list created above
bcftools view -O z -R sites_homref_homalt_het.txt monkparakeet_bi_snps.NOGTDP3.Q20.SAMP4.MAF0.5.vcf.gz > monkparakeet_bi_snps.NOGTDP3.Q20.SAMP4.MAF0.5.homref.alt.het.vcf.gz
bcftools index monkparakeet_bi_snps.NOGTDP3.Q20.SAMP4.MAF0.5.homref.alt.het.vcf.gz

# keep only autosomal sites
bcftools view -O z -r CM029917.1,CM029918.1,CM029919.1,CM029920.1,CM029921.1,CM029922.1,CM029923.1,CM029924.1,CM029925.1,CM029926.1,CM029927.1,CM029928.1,CM029929.1,CM029930.1,CM029931.1,CM029932.1,CM029933.1,CM029934.1,CM029935.1,CM029936.1,CM029937.1,CM029938.1,CM029939.1  monkparakeet_bi_snps.NOGTDP3.Q20.SAMP4.MAF0.5.homref.alt.het.vcf.gz > monkparakeet_bi_snps.NOGTDP3.Q20.SAMP4.MAF0.5.homref.alt.het.chrom.vcf.gz
bcftools index monkparakeet_bi_snps.NOGTDP3.Q20.SAMP4.MAF0.5.homref.alt.het.chrom.vcf.gz

# randomly sample 108 SNP sites and assign to file
bcftools view -H monkparakeet_bi_snps.NOGTDP3.Q20.SAMP4.MAF0.5.homref.alt.het.chrom.vcf.gz | shuf -n 108 | awk '{ printf "%5s\t%s\n", $1, $2 }' > 108_SNPs.txt

# minus 250bp and add 250bp for range to make bed file
awk -v s=250 '{print $1, $2-s, $2+s}' 108_SNPs.txt | sed 's/ /\t/g' > 108_SNPs.bed

# now extract these sequences from genome file
bedtools getfasta -fi ../genome/GCA_017639245.1_MMon_1.0_genomic.fna -bed 108_SNPs.bed > ../primer_design/108_SNPs.fasta
