#!/bin/bash

#SBATCH --job-name=H3K36me3DTplot_transcript.sh
#SBATCH --account=pi-yangili1
#SBATCH --time=24:00:00
#SBATCH --output=H3K36me3DTplot_transcript.out
#SBATCH --error=H3K36me3DTplot_transcript.err
#SBATCH --partition=broadwl
#SBATCH --mem=36G
#SBATCH --mail-type=END


source ~/activate_anaconda.sh
conda activate comp_threeprime_env


computeMatrix scale-regions -S ../data/H3K36me3/MergedIndiv.H3k36me3.sort.bw   -R /project2/gilad/briana/genome_anotation_data/hg38_refseq_anno/hg38_ncbiRefseq_transcripts.sort.bed -b 10000  -a 10000 --skipZeros --transcript_id_designator 4 -out ../data/H3K36me3/H3k36me3PASTranscript.gz

plotHeatmap --sortRegions descend -m ../data/H3K36me3/H3k36me3PASTranscript.gz --plotTitle "Transcript H3L36me3" --heatmapHeight 4 --colorMap magma  --averageTypeSummaryPlot "mean" --startLabel "TSS"  --endLabel "TES" -out ../data/H3K36me3/H3k36me3PASTranscript.png
