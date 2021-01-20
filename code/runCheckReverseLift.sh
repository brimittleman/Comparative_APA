#!/bin/bash

#SBATCH --job-name=FilterReverseLift
#SBATCH --account=pi-yangili1
#SBATCH --time=24:00:00
#SBATCH --output=FilterReverseLift.out
#SBATCH --error=FilterReverseLift.err
#SBATCH --partition=broadwl
#SBATCH --mem=36G
#SBATCH --mail-type=END

source ~/activate_anaconda.sh
conda activate comp_threeprime_env


for i in $(ls /project2/gilad/briana/Comparative_APA/Chimp/data/RNAseq/sort/*sort.bam.junc)
do
describer=$(echo ${i} |  cut -d - -f 1-5 )
Rscript ReverseLiftFilter.R -J ${describer}-sort.bam.junc -L ${describer}-sort.bam.junc.2Human.backChimp -O ${describer}-sort.bam.junc.2Human.backChimp.SamePlace
done

for i in $(ls /project2/gilad/briana/Comparative_APA/Human/data/RNAseq/sort/*sort.bam.junc)
do
describer=$(echo ${i} |  cut -d - -f 1-5 )
Rscript ReverseLiftFilter.R -J ${describer}-sort.bam.junc -L ${describer}-sort.bam.junc.2Chimp.backHuman.junc -O ${describer}-sort.bam.junc.2Chimp.backHuman.junc.SamePlace
done
