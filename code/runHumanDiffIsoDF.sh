#!/bin/bash

#SBATCH --job-name=run_Humanleafcutter_dF
#SBATCH --account=pi-yangili1
#SBATCH --time=24:00:00
#SBATCH --output=run_Humanleafcutter_dF.out
#SBATCH --error=run_Humanleafcutter_dF.err
#SBATCH --partition=broadwl
#SBATCH --mem=50G
#SBATCH --mail-type=END


module load Anaconda3/5.3.0
source  activate leafcutter


for i in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
do
Rscript /project2/gilad/briana/davidaknowles-leafcutter-c3d9474/scripts/leafcutter_ds.R --num_threads 4  ../Human/data/DiffIso_Human_DF/ALLPAS_postLift_LocParsed_Human_fixed4LC_chr${i}.fc ../Human/data/DiffIso_Human/sample_groups.txt -o ../Human/data/DiffIso_Human_DF/TN_diff_isoform_chr${i}.txt
done
