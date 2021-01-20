#!/bin/bash

#SBATCH --job-name=bam2junc_remove
#SBATCH --account=pi-yangili1
#SBATCH --time=24:00:00
#SBATCH --output=bam2junc_remove.out
#SBATCH --error=bam2junc_remove.err
#SBATCH --partition=broadwl
#SBATCH --mem=16G
#SBATCH --mail-type=END

#source deactivate
source deactivate

#module load Anaconda3/5.3.0
source  activate leafcutter
module load samtools


for bamfile in `ls /project2/gilad/briana/Comparative_APA/Human/data/RNAseq/sort_removebad/*.bam`
do
    echo Converting $bamfile to $bamfile.junc
    sh /project2/gilad/briana/leafcutter/scripts/bam2junc.sh $bamfile $bamfile.junc
    echo $bamfile.junc >> /project2/gilad/briana/Comparative_APA/Human/data/RNAseq/DiffSplice_removeBad/human_juncfiles.txt
done


for bamfile in `ls /project2/gilad/briana/Comparative_APA/Chimp/data/RNAseq/sort_removebad/*.bam`
do
    echo Converting $bamfile to $bamfile.junc
    sh /project2/gilad/briana/leafcutter/scripts/bam2junc.sh $bamfile $bamfile.junc
    echo $bamfile.junc >> /project2/gilad/briana/Comparative_APA/Chimp/data/RNAseq/DiffSplice_removeBad/chimp_juncfiles.txt
done
