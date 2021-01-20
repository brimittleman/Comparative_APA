#!/bin/bash

#SBATCH --job-name=nuclearTranscriptDTplot.sh
#SBATCH --account=pi-yangili1
#SBATCH --time=24:00:00
#SBATCH --output=nuclearTranscriptDTplot.out
#SBATCH --error=nuclearTranscriptDTplot.err
#SBATCH --partition=bigmem2
#SBATCH --mem=100G
#SBATCH --mail-type=END


source ~/activate_anaconda.sh
conda activate comp_threeprime_env




computeMatrix scale-regions -S ../Human/data/mergedbw_byFrac/human_Nuclear.SamplesMerged.sort.bw  -R /project2/gilad/briana/genome_anotation_data/hg38_refseq_anno/hg38_ncbiRefseq_transcripts.sort.bed -b 500  -a 500 --skipZeros --transcript_id_designator 4 -out ../data/DTmatrix/HumanNuclear_transcripts.gz

plotHeatmap --sortRegions descend -m ../data/DTmatrix/HumanNuclear_transcripts.gz --plotTitle "Human 3' Seq" --heatmapHeight 4 --colorMap magma  --averageTypeSummaryPlot "mean"  -out ../output/dtPlots/HumanNuclear_transcripts.png


computeMatrix scale-regions -S ../Chimp/data/mergedbw_byFrac/chimp_Nuclear.SamplesMerged.sort.bw  -R /project2/gilad/briana/genome_anotation_data/Chimp_refseqAnno/pantro6_ncbiRefseq_transcripts.sort.bed -b 500  -a 500 --skipZeros --transcript_id_designator 4 -out ../data/DTmatrix/ChimpNuclear_transcripts.gz

plotHeatmap --sortRegions descend -m ../data/DTmatrix/ChimpNuclear_transcripts.gz --plotTitle "Chimp 3' Seq" --heatmapHeight 4 --colorMap magma  --averageTypeSummaryPlot "mean"  -out ../output/dtPlots/ChimpNuclear_transcripts.png
