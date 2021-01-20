#!/bin/bash

#SBATCH --job-name=threeprimeOrthoFC
#SBATCH --account=pi-yangili1
#SBATCH --time=24:00:00
#SBATCH --output=threeprimeOrthoFC.out
#SBATCH --error=threeprimeOrthoFCcd.err
#SBATCH --partition=broadwl
#SBATCH --mem=36G
#SBATCH --mail-type=END




source ~/activate_anaconda.sh
conda activate comp_threeprime_env


#human
featureCounts -s 2 -a /project2/gilad/kenneth/OrthoExonPartialMapping/human.withM.gtf -t 'exon' -g 'gene_id' -o ../data/Threeprime2Ortho/Human_threeprime_orthoexon ../Human/data/sort_clean/*N*sort.bam


#human
featureCounts -s 2 -a /project2/gilad/kenneth/OrthoExonPartialMapping/chimp.withM.gtf -t 'exon' -g 'gene_id' -o ../data/Threeprime2Ortho/Chimp_threeprime_orthoexon ../Chimp/data/sort_clean/*N*sort.bam
