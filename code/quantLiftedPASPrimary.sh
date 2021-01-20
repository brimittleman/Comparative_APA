#!/bin/bash

#SBATCH --job-name=quantLiftedPASPrimary
#SBATCH --account=pi-yangili1
#SBATCH --time=24:00:00
#SBATCH --output=quantLiftedPASPrimary.out
#SBATCH --error=quantLiftedPASPrimary.err
#SBATCH --partition=broadwl
#SBATCH --mem=36G
#SBATCH --mail-type=END

source ~/activate_anaconda.sh
conda activate comp_threeprime_env




featureCounts -O -M --primary -a ../data/cleanPeaks_anno/AllPAS_postLift.sort_LocAnnoPARSED_chimpLoc.SAF -F SAF -o ../data/CleanLiftedPeaks_FC_primary/ALLPAS_postLift_LocParsed_Primary_Chimp ../Chimp/data/sort_clean/*.bam -s 1


featureCounts -O -M --primary  -a ../data/cleanPeaks_anno/AllPAS_postLift.sort_LocAnnoPARSED.SAF -F SAF -o ../data/CleanLiftedPeaks_FC_primary/ALLPAS_postLift_LocParsed_Primary_Human ../Human/data/sort_clean/*.bam -s 1
