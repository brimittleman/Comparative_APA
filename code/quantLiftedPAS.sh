#!/bin/bash

#SBATCH --job-name=quantLiftedPAS
#SBATCH --account=pi-yangili1
#SBATCH --time=24:00:00
#SBATCH --output=quantLiftedPAS.out
#SBATCH --error=quantLiftedPAS.err
#SBATCH --partition=broadwl
#SBATCH --mem=36G
#SBATCH --mail-type=END

source ~/activate_anaconda.sh
conda activate comp_threeprime_env




featureCounts -O -a ../data/cleanPeaks_anno/AllPAS_postLift.sort_LocAnnoPARSED_chimpLoc.SAF -F SAF -o ../Chimp/data/CleanLiftedPeaks_FC/ALLPAS_postLift_LocParsed_Chimp ../Chimp/data/sort_clean/*.bam -s 1


featureCounts -O -a ../data/cleanPeaks_anno/AllPAS_postLift.sort_LocAnnoPARSED.SAF -F SAF -o ../Human/data/CleanLiftedPeaks_FC/ALLPAS_postLift_LocParsed_Human ../Human/data/sort_clean/*.bam -s 1
