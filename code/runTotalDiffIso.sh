#!/bin/bash

#SBATCH --job-name=run_Totalleafcutter_ds
#SBATCH --account=pi-yangili1
#SBATCH --time=24:00:00
#SBATCH --output=run_Totalleafcutter_ds.out
#SBATCH --error=run_Totalleafcutter_ds.err
#SBATCH --partition=broadwl
#SBATCH --mem=50G
#SBATCH --mail-type=END


module load Anaconda3/5.3.0
source  activate leafcutter


for i in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
do
Rscript /project2/gilad/briana/davidaknowles-leafcutter-c3d9474/scripts/leafcutter_ds.R --num_threads 4  ../data/DiffIso_Total/ALLPAS_postLift_LocParsed_HvC_Total_fixed4LC_chr${i}.fc ../data/TotalHvC/sample_goups.txt -o ../data/DiffIso_Total/TN_diff_isoform_chr${i}.txt
done
