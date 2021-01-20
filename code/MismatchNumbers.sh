#!/bin/bash

#SBATCH --job-name=MismatchNumbers
#SBATCH --time=24:00:00
#SBATCH --partition=broadwl
#SBATCH --mem=32G
#SBATCH --mail-type=END
#SBATCH --output=MismatchNumbers.out
#SBATCH --error=MismatchNumbers.err

source ~/activate_anaconda.sh
conda activate comp_threeprime_env


for i in $(ls ../data/TestMM2_SeondaryRead/*.sort.bam)
do
sample=$(echo ${i} | cut -d "/" -f 4  | cut -f 1 -d ".")
samtools view $i | cut -f15 | sort -n | uniq -c > ../data/TestMM2_mismatch/${sample}_mismatch.txt
done

for i in $(ls ../data/TestMM2_PrimaryRead/*.sort.bam)
do
sample=$(echo ${i} | cut -d "/" -f 4  | cut -f 1 -d ".")
samtools view $i | cut -f15 | sort -n | uniq -c > ../data/TestMM2_mismatch/${sample}_mismatch.txt
done
