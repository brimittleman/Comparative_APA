#!/bin/bash

#SBATCH --job-name=generateStarIndex
#SBATCH --account=pi-yangili1
#SBATCH --time=24:00:00
#SBATCH --output=generateStarIndex.out
#SBATCH --error=generateStarIndex.err
#SBATCH --partition=bigmem2
#SBATCH --mem=200G
#SBATCH --mail-type=END

source ~/activate_anaconda.sh
conda activate comp_threeprime_env


STAR --runThreadN 4 --runMode genomeGenerate --genomeChrBinNbits = 15  --genomeDir /project2/gilad/briana/genome_anotation_data/Pantro5/ --genomeFastaFiles /project2/gilad/kenneth/References/chimp/genome/panTro5.fa --sjdbGTFfile /project2/gilad/briana/genome_anotation_data/Pantro5/Pan_troglodytes.Pan_tro_3.0.93.PANTRO5Version.gtf
