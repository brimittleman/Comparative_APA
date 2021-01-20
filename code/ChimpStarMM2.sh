#!/bin/bash

#SBATCH --job-name=ChimpStarMM2
#SBATCH --time=24:00:00
#SBATCH --partition=bigmem2
#SBATCH --mem=100G
#SBATCH --mail-type=END
#SBATCH --output=ChimpStarMM2.out
#SBATCH --error=ChimpStarMM2.err



source ~/activate_anaconda.sh
conda activate comp_threeprime_env



for i in $(ls ../Chimp/data/fastq/*)
do
sample=$(echo ${i} | cut -d "/" -f 5  | cut -f 1 -d ".")
STAR --runThreadN 4 --genomeDir  /project2/gilad/briana/genome_anotation_data/Chimp_genome/ --readFilesIn $i --outFilterMultimapNmax 10 --outSAMmultNmax 2 --sjdbGTFfile /project2/gilad/briana/genome_anotation_data/Chimp_genome/pantro6_ncbiRefseq.gtf --outSAMtype BAM SortedByCoordinate --outStd BAM_SortedByCoordinate --outFileNamePrefix ../data/TestMM2/${sample}.MM2.bam > ../data/TestMM2/${sample}.MM2.bam
samtools sort -o ../data/TestMM2/${sample}.MM2.sort.bam -O bam ../data/TestMM2/${sample}.MM2.bam
samtools index ../data/TestMM2/${sample}.MM2.sort.bam
done
