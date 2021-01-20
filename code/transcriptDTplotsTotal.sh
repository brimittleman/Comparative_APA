#!/bin/bash

#SBATCH --job-name=TotalTranscriptDTplot.sh
#SBATCH --account=pi-yangili1
#SBATCH --time=24:00:00
#SBATCH --output=TotalTranscriptDTplot.out
#SBATCH --error=TotalTranscriptDTplot.err
#SBATCH --partition=bigmem2
#SBATCH --mem=100G
#SBATCH --mail-type=END


source ~/activate_anaconda.sh
conda activate comp_threeprime_env




computeMatrix scale-regions -S ../Human/data/mergedbw_byFrac/human_Total.SamplesMerged.sort.bw  -R /project2/gilad/briana/genome_anotation_data/hg38_refseq_anno/hg38_ncbiRefseq_transcripts.sort.bed -b 500  -a 500 --skipZeros --transcript_id_designator 4 -out ../data/DTmatrix/HumanTotal_transcripts.gz

plotHeatmap --sortRegions descend -m ../data/DTmatrix/HumanTotal_transcripts.gz --plotTitle "Human Total" --heatmapHeight 7 --colorMap YlGnBu  --averageTypeSummaryPlot "mean"  -out ../output/dtPlots/HumanTotal_transcripts.png


computeMatrix scale-regions -S ../Chimp/data/mergedbw_byFrac/chimp_Total.SamplesMerged.sort.bw  -R /project2/gilad/briana/genome_anotation_data/Chimp_refseqAnno/pantro6_ncbiRefseq_transcripts.sort.bed -b 500  -a 500 --skipZeros --transcript_id_designator 4 -out ../data/DTmatrix/ChimpTotal_transcripts.gz

plotHeatmap --sortRegions descend -m ../data/DTmatrix/ChimpTotal_transcripts.gz --plotTitle "Chimp Total" --heatmapHeight 7 --colorMap YlGnBu  --averageTypeSummaryPlot "mean"  -out ../output/dtPlots/ChimpTotal_transcripts.png
