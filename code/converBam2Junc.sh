#!/bin/bash

#SBATCH --job-name=bam2junc
#SBATCH --account=pi-yangili1
#SBATCH --time=24:00:00
#SBATCH --output=bam2junc.out
#SBATCH --error=bam2junc.err
#SBATCH --partition=broadwl
#SBATCH --mem=16G
#SBATCH --mail-type=END

source ~/activate_anaconda.sh
conda activate comp_threeprime_env

for bamfile in `ls /project2/gilad/briana/Comparative_APA/Human/data/RNAseq/sort/*.bam`
do
  echo Converting $bamfile to $bamfile.junc
  samtools index $bamfile
  regtools junctions extract -a 8 -m 10 -M 500000 -s 0 $bamfile -o $bamfile.junc
  echo $bamfile.junc >> /project2/gilad/briana/Comparative_APA/Human/data/RNAseq/DiffSplice/human_juncfiles.txt
done


for bamfile in `ls /project2/gilad/briana/Comparative_APA/Chimp/data/RNAseq/sort/*.bam`
do
    echo Converting $bamfile to $bamfile.junc
    samtools index $bamfile
    regtools junctions extract -a 8 -m 10 -M 500000 -s 0 $bamfile -o $bamfile.junc
    echo $bamfile.junc >> /project2/gilad/briana/Comparative_APA/Chimp/data/RNAseq/DiffSplice/chimp_juncfiles.txt
done
