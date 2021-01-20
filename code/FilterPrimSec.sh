#!/bin/bash

#SBATCH --job-name=FilterPrimSec
#SBATCH --time=24:00:00
#SBATCH --partition=broadwl
#SBATCH --mem=32G
#SBATCH --mail-type=END
#SBATCH --output=FilterPrimSec.out
#SBATCH --error=FilterPrimSec.err

source ~/activate_anaconda.sh
conda activate comp_threeprime_env



for i in $(ls ../data/TestMM2/*.sort.bam)
do
sample=$(echo ${i} | cut -d "/" -f 4  | cut -f 1 -d ".")
python filterSecondaryread.py $i  ../data/TestMM2_SeondaryRead/${sample}_secondary.bam
samtools sort -o ../data/TestMM2_SeondaryRead/${sample}_secondary.sort.bam -O bam  ../data/TestMM2_SeondaryRead/${sample}_secondary.bam
samtools index ../data/TestMM2_SeondaryRead/${sample}_secondary.sort.bam
python filterPrimaryread.py $i  ../data/TestMM2_PrimaryRead/${sample}_primary.bam
samtools sort -o  ../data/TestMM2_PrimaryRead/${sample}_primary.sort.bam -O bam  ../data/TestMM2_PrimaryRead/${sample}_primary.bam
samtools index ../data/TestMM2_PrimaryRead/${sample}_primary.sort.bam
done
