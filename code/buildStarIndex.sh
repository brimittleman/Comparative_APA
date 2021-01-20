#!/bin/bash

#SBATCH --job-name=generateStarIndex
#SBATCH --account=pi-yangili1
#SBATCH --time=24:00:00
#SBATCH --output=generateStarIndexHuman.out
#SBATCH --error=generateStarIndexHuman.err
#SBATCH --partition=bigmem2
#SBATCH --mem=200G
#SBATCH --mail-type=END

source ~/activate_anaconda.sh
conda activate comp_threeprime_env



STAR --runThreadN 4 --runMode genomeGenerate --genomeDir /project2/gilad/briana/genome_anotation_data/hg38_try2 --genomeFastaFiles /project2/gilad/kenneth/References/human/genome/hg38.fa --sjdbGTFfile /project2/gilad/briana/genome_anotation_data/hg38/hg38_ncbiRefseq.gtf
