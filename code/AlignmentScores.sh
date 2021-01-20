#!/bin/bash

#SBATCH --job-name=AlignmentScores
#SBATCH --time=24:00:00
#SBATCH --partition=broadwl
#SBATCH --mem=32G
#SBATCH --mail-type=END
#SBATCH --output=AlignmentScores.out
#SBATCH --error=AlignmentScores.err

source ~/activate_anaconda.sh
conda activate comp_threeprime_env


for i in $(ls ../data/TestMM2/*.255reads.sort.bam )
do
sample=$(echo ${i} | cut -d "/" -f 4  | cut -f 1 -d ".")
samtools view $i | cut -f14 | sort -n | uniq -c > ../data/TestMM2_AS/${sample}_255_AS.txt
samtools view ../data/TestMM2_SeondaryRead/${sample}_secondary.sort.bam | cut -f14 | sort -n | uniq -c > ../data/TestMM2_AS/${sample}_Secondary.txt
samtools view ../data/TestMM2_PrimaryRead/${sample}_primary.sort.bam  | cut -f14 | sort -n | uniq -c > ../data/TestMM2_AS/${sample}_Primary.txt
done
