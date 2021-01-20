#!/bin/bash




#SBATCH --job-name=maphg19_sub
#SBATCH --time=24:00:00
#SBATCH --partition=bigmem2
#SBATCH --mem=200G
#SBATCH --tasks-per-node=4
#SBATCH --mail-type=END
#SBATCH --output=maphg19_sub.out
#SBATCH --error=maphg19_sub.err


source ~/activate_anaconda.sh
conda activate comp_threeprime_env



subjunc -i /project2/gilad/briana/genome_anotation_data/hg38/hg38_subread -r ../Human/data/fastq/human_combined_18498_N.fastq -T 8 > ../Human/data/bam_hg19_sub/human_combined_18498_N.bam

samtools sort ../Human/data/bam_hg19_sub/human_combined_18498_N.bam > ../Human/data/sort_hg19_sub/human_combined_18498_N.sort.bam
samtools index ../Human/data/sort_hg19_sub/human_combined_18498_N.sort.bam

subjunc -i /project2/gilad/briana/genome_anotation_data/hg38/hg38_subread -r ../Human/data/fastq/human_combined_18498_T.fastq -T 8 > ../Human/data/bam_hg19_sub/human_combined_18498_T.bam

samtools sort ../Human/data/bam_hg19_sub/human_combined_18498_T.bam > ../Human/data/sort_hg19_sub/human_combined_18498_T.sort.bam
samtools index ../Human/data/sort_hg19_sub/human_combined_18498_T.sort.bam


subjunc -i /project2/gilad/briana/genome_anotation_data/hg38/hg38_subread -r ../Human/data/fastq/human_combined_18499_N.fastq -T 8 > ../Human/data/bam_hg19_sub/human_combined_18499_N.bam

samtools sort ../Human/data/bam_hg19_sub/human_combined_18499_N.bam > ../Human/data/sort_hg19_sub/human_combined_18499_N.sort.bam
samtools index ../Human/data/sort_hg19_sub/human_combined_18499_N.sort.bam


subjunc -i /project2/gilad/briana/genome_anotation_data/hg38/hg38_subread -r ../Human/data/fastq/human_combined_18499_T.fastq -T 8 > ../Human/data/bam_hg19_sub/human_combined_18499_T.bam

samtools sort ../Human/data/bam_hg19_sub/human_combined_18499_T.bam > ../Human/data/sort_hg19_sub/human_combined_18499_T.sort.bam
samtools index ../Human/data/sort_hg19_sub/human_combined_18499_T.sort.bam


subjunc -i /project2/gilad/briana/genome_anotation_data/hg38/hg38_subread -r ../Human/data/fastq/human_combined_18502_N.fastq -T 8 > ../Human/data/bam_hg19_sub/human_combined_18502_N.bam

samtools sort ../Human/data/bam_hg19_sub/human_combined_18502_N.bam > ../Human/data/sort_hg19_sub/human_combined_18502_N.sort.bam
samtools index ../Human/data/sort_hg19_sub/human_combined_18502_N.sort.bam



subjunc -i /project2/gilad/briana/genome_anotation_data/hg38/hg38_subread -r ../Human/data/fastq/human_combined_18502_T.fastq -T 8 > ../Human/data/bam_hg19_sub/human_combined_18502_T.bam

samtools sort ../Human/data/bam_hg19_sub/human_combined_18502_T.bam > ../Human/data/sort_hg19_sub/human_combined_18502_T.sort.bam
samtools index ../Human/data/sort_hg19_sub/human_combined_18502_T.sort.bam
