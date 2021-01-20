#!/bin/bash

#SBATCH --job-name=mergedbamRNAand2bw
#SBATCH --account=pi-yangili1
#SBATCH --time=36:00:00
#SBATCH --output=mergedbamRNAand2bw.out
#SBATCH --error=mergedbamRNAand2bw.err
#SBATCH --cpus-per-task=2
#SBATCH --partition=bigmem2
#SBATCH --mem=216G
#SBATCH --mail-type=END

source ~/activate_anaconda.sh
conda activate comp_threeprime_env

#human
samtools merge ../Human/data/RNAseq/mergeBam/mergedHumanRNAseq.bam ../Human/data/RNAseq/sort/*sort.bam

samtools sort ../Human/data/RNAseq/mergeBam/mergedHumanRNAseq.bam -o ../Human/data/RNAseq/mergeBam/mergedHumanRNAseq.sort.bam

samtools index ../Human/data/RNAseq/mergeBam/mergedHumanRNAseq.sort.bam

bamCoverage -b ../Human/data/RNAseq/mergeBam/mergedHumanRNAseq.sort.bam -o  ../Human/data/RNAseq/mergeBW/mergedHumanRNAseq.bw


#chimp
samtools merge ../Chimp/data/RNAseq/mergeBam/mergedChimpRNAseq.bam ../Chimp/data/RNAseq/sort/*sort.bam

samtools sort ../Chimp/data/RNAseq/mergeBam/mergedChimpRNAseq.bam -o ../Chimp/data/RNAseq/mergeBam/mergedChimpRNAseq.sort.bam

samtools index ../Chimp/data/RNAseq/mergeBam/mergedChimpRNAseq.sort.bam

bamCoverage -b ../Chimp/data/RNAseq/mergeBam/mergedChimpRNAseq.sort.bam -o  ../Chimp/data/RNAseq/mergeBW/mergedChimpRNAseq.bw
