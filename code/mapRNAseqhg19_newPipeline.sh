#!/bin/bash


#SBATCH --job-name=maphg19_new
#SBATCH --time=24:00:00
#SBATCH --partition=bigmem2
#SBATCH --mem=200G
#SBATCH --tasks-per-node=4
#SBATCH --mail-type=END
#SBATCH --output=maphg19_new.out
#SBATCH --error=maphg19_new.err



module load STAR
module load samtools



STAR --runThreadN 4  --twopassMode Basic --genomeDir /project2/gilad/briana/genome_anotation_data/star_genome/ --readFilesIn ../Human/data/RNAseq/fastq/YG-BM-S10-18504H-Total_S10_R1_001.fastq   --outFilterMultimapNmax 10 --outSAMmultNmax 1 --outSAMtype BAM SortedByCoordinate --outStd BAM_SortedByCoordinate --outFileNamePrefix ../Human/data/RNAseq/bam_hg19/NA18504_hg19_new > ../Human/data/RNAseq/bam_hg19/NA18504_hg19_new.bam

samtools sort ../Human/data/RNAseq/bam_hg19/NA18504_hg19_new.bam > ../Human/data/RNAseq/bam_hg19/NA18504_hg19_new-sort.bam
samtools index ../Human/data/RNAseq/bam_hg19/NA18504_hg19_new-sort.bam

STAR --runThreadN 4  --twopassMode Basic --genomeDir /project2/gilad/briana/genome_anotation_data/star_genome/ --readFilesIn ../Human/data/RNAseq/fastq/YG-BM-S11-18510H-Total_S11_R1_001.fastq  --outFilterMultimapNmax 10 --outSAMmultNmax 1 --outSAMtype BAM SortedByCoordinate --outStd BAM_SortedByCoordinate --outFileNamePrefix ../Human/data/RNAseq/bam_hg19/NA18510_hg19_new > ../Human/data/RNAseq/bam_hg19/NA18510_hg19_new.bam

samtools sort ../Human/data/RNAseq/bam_hg19/NA18510_hg19_new.bam > ../Human/data/RNAseq/bam_hg19/NA18510_hg19_new-sort.bam
samtools index ../Human/data/RNAseq/bam_hg19/NA18510_hg19_new-sort.bam
