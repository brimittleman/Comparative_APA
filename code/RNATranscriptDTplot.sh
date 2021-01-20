#!/bin/bash

#SBATCH --job-name=RNATranscriptDTplot.sh
#SBATCH --account=pi-yangili1
#SBATCH --time=24:00:00
#SBATCH --output=RNATranscriptDTplot.out
#SBATCH --error=RNATranscriptDTplot.err
#SBATCH --partition=bigmem2
#SBATCH --mem=100G
#SBATCH --mail-type=END


source ~/activate_anaconda.sh
conda activate comp_threeprime_env


computeMatrix scale-regions -S ../Human/data/RNAseq/mergeBW/mergedHumanRNAseq.bw  -R /project2/gilad/briana/genome_anotation_data/hg38_refseq_anno/hg38_ncbiRefseq_transcripts.sort.bed -b 500  -a 500 --skipZeros --transcript_id_designator 4 -out ../data/DTmatrix/HumanRNA_transcripts.gz

plotHeatmap --sortRegions descend -m ../data/DTmatrix/HumanRNA_transcripts.gz --plotTitle "Human RNA" --heatmapHeight 4 --colorMap magma  --averageTypeSummaryPlot "mean"  -out ../output/dtPlots/HumanRNA_transcripts.png


computeMatrix scale-regions -S ../Chimp/data/RNAseq/mergeBW/mergedChimpRNAseq.bw  -R /project2/gilad/briana/genome_anotation_data/Chimp_refseqAnno/pantro6_ncbiRefseq_transcripts.sort.bed -b 500  -a 500 --skipZeros --transcript_id_designator 4 -out ../data/DTmatrix/ChimpRNA_transcripts.gz

plotHeatmap --sortRegions descend -m ../data/DTmatrix/ChimpRNA_transcripts.gz --plotTitle "Chimp RNA" --heatmapHeight 4 --colorMap magma  --averageTypeSummaryPlot "mean"  -out ../output/dtPlots/ChimpRNA_transcripts.png
