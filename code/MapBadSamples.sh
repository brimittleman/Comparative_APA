#!/bin/bash

#SBATCH --job-name=MapBadSamples
#SBATCH --account=pi-yangili1
#SBATCH --time=24:00:00
#SBATCH --output=MapBadSamples.out
#SBATCH --error=MapBadSamples.err
#SBATCH --partition=bigmem2
#SBATCH --mem=100G
#SBATCH --mail-type=END

source ~/activate_anaconda.sh
conda activate comp_threeprime_env


#map 4973 human

STAR --runThreadN 4 --twopassMode Basic --genomeDir /project2/gilad/briana/genome_anotation_data/hg38_try2/ --readFilesIn ../Chimp/data/RNAseq/fastq/YG-BM-S5-4973C-Total_S5_R1_001.fastq --outFilterMultimapNmax 10 --outSAMmultNmax 1 --sjdbGTFfile  /project2/gilad/briana/genome_anotation_data/hg38_try2/hg38.gtf --outSAMtype BAM SortedByCoordinate --outStd BAM_SortedByCoordinate --outFileNamePrefix ../data/TwoBadSampleAnalysis/NA4973inHuman > ../data/TwoBadSampleAnalysis/NA4973inHuman


#18493 to chimp

STAR --runThreadN 4 --twopassMode Basic --genomeDir /project2/gilad/briana/genome_anotation_data/Chimp_genome/ --readFilesIn ../Human/data/RNAseq/fastq/YG-BM-S7-18498H-Total_S7_R1_001.fastq --outFilterMultimapNmax 10 --outSAMmultNmax 1 --sjdbGTFfile  /project2/gilad/briana/genome_anotation_data/Chimp_genome/pantro5.gtf --outSAMtype BAM SortedByCoordinate --outStd BAM_SortedByCoordinate --outFileNamePrefix ../data/TwoBadSampleAnalysis/NA18498inChimp > ../data/TwoBadSampleAnalysis/NA18498inChimp
