#!/bin/bash




#SBATCH --job-name=maphg19
#SBATCH --time=24:00:00
#SBATCH --partition=bigmem2
#SBATCH --mem=200G
#SBATCH --tasks-per-node=4
#SBATCH --mail-type=END
#SBATCH --output=maphg19.out
#SBATCH --error=maphg19.err



module load STAR



STAR --runThreadN 4 --genomeDir /project2/gilad/briana/genome_anotation_data/star_genome/ --readFilesIn ../Human/data/fastq/human_combined_18498_N.fastq --outFilterMultimapNmax 1 --outSAMtype BAM SortedByCoordinate --outStd BAM_SortedByCoordinate --outFileNamePrefix ../Human/data/bam_hg19/human_combined_18498_N.bam > ../Human/data/bam_hg19/human_combined_18498_N.bam

samtools sort ../Human/data/bam_hg19/human_combined_18498_N.bam > ../Human/data/sort_hg19/human_combined_18498_N.sort.bam
samtools index ../Human/data/sort_hg19/human_combined_18498_N.sort.bam

STAR --runThreadN 4 --genomeDir /project2/gilad/briana/genome_anotation_data/star_genome/ --readFilesIn ../Human/data/fastq/human_combined_18498_T.fastq --outFilterMultimapNmax 1 --outSAMtype BAM SortedByCoordinate --outStd BAM_SortedByCoordinate --outFileNamePrefix ../Human/data/bam_hg19/human_combined_18498_T.bam > ../Human/data/bam_hg19/human_combined_18498_T.bam

samtools sort ../Human/data/bam_hg19/human_combined_18498_T.bam > ../Human/data/sort_hg19/human_combined_18498_T.sort.bam
samtools index ../Human/data/sort_hg19/human_combined_18498_T.sort.bam

STAR --runThreadN 4 --genomeDir /project2/gilad/briana/genome_anotation_data/star_genome/ --readFilesIn ../Human/data/fastq/human_combined_18499_N.fastq --outFilterMultimapNmax 1 --outSAMtype BAM SortedByCoordinate --outStd BAM_SortedByCoordinate --outFileNamePrefix ../Human/data/bam_hg19/human_combined_18499_N.bam > ../Human/data/bam_hg19/human_combined_18499_N.bam

samtools sort ../Human/data/bam_hg19/human_combined_18499_N.bam > ../Human/data/sort_hg19/human_combined_18499_N.sort.bam
samtools index ../Human/data/sort_hg19/human_combined_18499_N.sort.bam

STAR --runThreadN 4 --genomeDir /project2/gilad/briana/genome_anotation_data/star_genome/ --readFilesIn ../Human/data/fastq/human_combined_18499_T.fastq --outFilterMultimapNmax 1 --outSAMtype BAM SortedByCoordinate --outStd BAM_SortedByCoordinate --outFileNamePrefix ../Human/data/bam_hg19/human_combined_18499_T.bam > ../Human/data/bam_hg19/human_combined_18499_T.bam


samtools sort ../Human/data/bam_hg19/human_combined_18499_T.bam > ../Human/data/sort_hg19/human_combined_18499_T.sort.bam
samtools index ../Human/data/sort_hg19/human_combined_18499_T.sort.bam


STAR --runThreadN 4 --genomeDir /project2/gilad/briana/genome_anotation_data/star_genome/ --readFilesIn ../Human/data/fastq/human_combined_18502_N.fastq --outFilterMultimapNmax 1 --outSAMtype BAM SortedByCoordinate --outStd BAM_SortedByCoordinate --outFileNamePrefix ../Human/data/bam_hg19/human_combined_18502_N.bam > ../Human/data/bam_hg19/human_combined_18502_N.bam

samtools sort ../Human/data/bam_hg19/human_combined_18502_N.bam > ../Human/data/sort_hg19/human_combined_18502_N.sort.bam
samtools index ../Human/data/sort_hg19/human_combined_18502_N.sort.bam

STAR --runThreadN 4 --genomeDir /project2/gilad/briana/genome_anotation_data/star_genome/ --readFilesIn ../Human/data/fastq/human_combined_18502_T.fastq --outFilterMultimapNmax 1 --outSAMtype BAM SortedByCoordinate --outStd BAM_SortedByCoordinate --outFileNamePrefix ../Human/data/bam_hg19/human_combined_18502_T.bam > ../Human/data/bam_hg19/human_combined_18502_T.bam

samtools sort ../Human/data/bam_hg19/human_combined_18502_T.bam > ../Human/data/sort_hg19/human_combined_18502_T.sort.bam
samtools index ../Human/data/sort_hg19/human_combined_18502_T.sort.bam
