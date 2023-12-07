#!/bin/bash

#Settings for the Sun Grid Engine

# run time for job in hours:mins:sec (max 168:0:0, jobs with h_rt < 8:0:0 have priority)

#$ -l h_rt=47:59:59

#$ -pe openmp 1

# request memory for job (default 2G, max 256G)

#$ -l rmem=24G

#$ -q evolgen.q

#$ -P evolgen

# give the job a name (optional):

#$-N call_snps

source /usr/local/extras/Genomics/.bashrc

./08_call_snps.sh /fastdata/bo4kma/monkparakeet
