#!/bin/bash

#SBATCH --job-name=NuclearDTUTR
#SBATCH --account=pi-yangili1
#SBATCH --time=24:00:00
#SBATCH --output=NuclearDTUTRt.out
#SBATCH --error=NuclearDTUTR.err
#SBATCH --partition=broadwl
#SBATCH --mem=50G
#SBATCH --mail-type=END


source ~/activate_anaconda.sh
conda activate comp_threeprime_env

computeMatrix scale-regions -S ../Human/data/mergedbw_byFrac/human_Nuclear.SamplesMerged.sort.bw  -R ../data/orthoUTR/HumanDistal3UTR.sort.bed -b 500  -a 500 --skipZeros --transcript_id_designator 4 -out ../data/orthoUTR/HumanNuclear_UTR.gz

plotHeatmap --sortRegions descend -m ../data/orthoUTR/HumanNuclear_UTR.gz --plotTitle "Human" --heatmapHeight 4 --colorMap magma  --averageTypeSummaryPlot "mean" --startLabel "5'"  --endLabel "3'" -out ../data/orthoUTR/Human_UTR.png


computeMatrix scale-regions -S ../Chimp/data/mergedbw_byFrac/chimp_Nuclear.SamplesMerged.sort.bw  -R ../data/orthoUTR/ChimpDistal3UTR.sort.bed -b 500  -a 500 --skipZeros --transcript_id_designator 4 -out ../data/orthoUTR/ChimpNuclear_UTR.gz

plotHeatmap --sortRegions descend -m ../data/orthoUTR/ChimpNuclear_UTR.gz --plotTitle "Chimp" --heatmapHeight 4 --colorMap magma  --averageTypeSummaryPlot "mean"  --startLabel "5'"  --endLabel "3'"  -out ../data/orthoUTR/Chimp_UTR.png
