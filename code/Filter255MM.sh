#!/bin/bash

#SBATCH --job-name=Filter255
#SBATCH --time=24:00:00
#SBATCH --partition=broadwl
#SBATCH --mem=32G
#SBATCH --mail-type=END
#SBATCH --output=Filter255.out
#SBATCH --error=Filter255.err

source ~/activate_anaconda.sh
conda activate comp_threeprime_env



for i in $(ls ../data/TestMM2/*.sort.bam)
do
sample=$(echo ${i} | cut -d "/" -f 4  | cut -f 1 -d ".")
samtools view -q 255 -b $i >  ../data/TestMM2/${sample}.255reads.bam
samtools sort -o ../data/TestMM2/${sample}.255reads.sort.bam -O bam  ../data/TestMM2/${sample}.255reads.bam
samtools index  ../data/TestMM2/${sample}.255reads.sort.bam
samtools view ../data/TestMM2/${sample}.255reads.sort.bam | cut -f15 | sort -n | uniq -c > ../data/TestMM2_mismatch/${sample}_255mismatch.txt
done
