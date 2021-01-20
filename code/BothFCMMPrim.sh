#!/bin/bash

#SBATCH --job-name=BothFCMMPrim
#SBATCH --time=24:00:00
#SBATCH --partition=broadwl
#SBATCH --mem=32G
#SBATCH --mail-type=END
#SBATCH --output=BothFCMMPrim.out
#SBATCH --error=BothFCMMPrim.err



source ~/activate_anaconda.sh
conda activate comp_threeprime_env


featureCounts -O --primary -M -a ../data/cleanPeaks_anno/AllPAS_postLift.sort_LocAnnoPARSED_chimpLoc.SAF -F SAF -o ../data/testQuant/ALLPAS_postLift_LocParsed_Chimp_PrimaryandM ../Chimp/data/sort_clean/*_N-clean.sort.bam  -s 1

featureCounts -O --primary -M  -a ../data/cleanPeaks_anno/AllPAS_postLift.sort_LocAnnoPARSED.SAF -F SAF -o ../data/testQuant/ALLPAS_postLift_LocParsed_Human_PrimaryandM ../Human/data/sort_clean/*_N-clean.sort.bam -s 1
