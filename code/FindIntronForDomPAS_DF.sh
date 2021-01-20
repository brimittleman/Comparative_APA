#!/bin/bash

#SBATCH --job-name=FindIntronForDomPAS
#SBATCH --account=pi-yangili1
#SBATCH --time=24:00:00
#SBATCH --output=FindIntronForDomPAS.out
#SBATCH --error=FindIntronForDomPAS.err
#SBATCH --partition=broadwl
#SBATCH --mem=36G
#SBATCH --mail-type=END

source ~/activate_anaconda.sh
conda activate comp_threeprime_env



bedtools intersect -S -sorted -loj -a ../data/DominantPAS_DF/SameDominantPAS_intronic.bed -b /project2/gilad/briana/genome_anotation_data/hg38_refseq_anno/hg38_refSeq_intron_sorted.bed  > ../data/DominantPAS_DF/SameDominantPAS_intronic_mapped2Intron.txt


bedtools intersect -S -sorted -loj -a ../data/DominantPAS_DF/DifferentDominantPAS_intronic.bed -b /project2/gilad/briana/genome_anotation_data/hg38_refseq_anno/hg38_refSeq_intron_sorted.bed  > ../data/DominantPAS_DF/DifferentDominantPAS_intronic_mapped2Intron.txt
