#!/bin/bash

#SBATCH --job-name=ConvertJunc2Bed
#SBATCH --account=pi-yangili1
#SBATCH --time=24:00:00
#SBATCH --output=ConvertJunc2Bed.out
#SBATCH --error=ConvertJunc2Bed.err
#SBATCH --partition=broadwl
#SBATCH --mem=36G
#SBATCH --mail-type=END



for i in $(ls /project2/gilad/briana/Comparative_APA/Chimp/data/RNAseq/bam/*.bamSJ.out.tab)
do
/project2/yangili1/bjf79/201911_DHX38/rna-seq-DHX38/code/scripts/SJ_to_junctions.sh $i $i.junc.bed
done

for i in $(ls /project2/gilad/briana/Comparative_APA/Human/data/RNAseq/bam/*.bamSJ.out.tab)
do
/project2/yangili1/bjf79/201911_DHX38/rna-seq-DHX38/code/scripts/SJ_to_junctions.sh $i $i.junc.bed
done
