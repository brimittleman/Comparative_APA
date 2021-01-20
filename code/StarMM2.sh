#!/bin/bash

#SBATCH --job-name=StarMM2
#SBATCH --time=24:00:00
#SBATCH --partition=bigmem2
#SBATCH --mem=100G
#SBATCH --mail-type=END
#SBATCH --output=StarMM2.out
#SBATCH --error=StarMM2.err



source ~/activate_anaconda.sh
conda activate comp_threeprime_env

#"STAR --runThreadN 4 --genomeDir {params.genome} --readFilesIn {input} --outFilterMultimapNmax 10 --outSAMmultNmax 1 --sjdbGTFfile {params.gtffile} --outSAMtype BAM SortedByCoordinate --outStd BAM_SortedByCoordinate --outFileNamePrefix {output} > {output}"
#/project2/gilad/briana/genome_anotation_data/hg38_try2/hg38.gtf
# /project2/gilad/kenneth/References/human/genome/hg38.fa

for i in $(ls ../Human/data/fastq/*)
do
sample=$(echo ${i} | cut -d "/" -f 5  | cut -f 1 -d ".")
STAR --runThreadN 4 --genomeDir  /project2/gilad/briana/genome_anotation_data/hg38_try2/ --readFilesIn $i --outFilterMultimapNmax 10 --outSAMmultNmax 2 --sjdbGTFfile /project2/gilad/briana/genome_anotation_data/hg38_try2/hg38.gtf --outSAMtype BAM SortedByCoordinate --outStd BAM_SortedByCoordinate --outFileNamePrefix ../data/TestMM2/${sample}.MM2.bam > ../data/TestMM2/${sample}.MM2.bam
samtools sort -o ../data/TestMM2/${sample}.MM2.sort.bam -O bam ../data/TestMM2/${sample}.MM2.bam
samtools index ../data/TestMM2/${sample}.MM2.sort.bam
done
