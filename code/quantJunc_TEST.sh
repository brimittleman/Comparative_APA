#!/bin/bash

#SBATCH --job-name=quatJunc
#SBATCH --account=pi-yangili1
#SBATCH --time=24:00:00
#SBATCH --output=quatJunc.out
#SBATCH --error=quatJunc.err
#SBATCH --partition=broadwl
#SBATCH --mem=16G
#SBATCH --mail-type=END


source ~/activate_anaconda_python2.sh
#conda activate comp_threeprime_env
#source  activate leafcutter


touch ../data/DiffSplice_liftedJunc/Test_juncfiles.txt
for juncfile in `ls /project2/gilad/briana/Comparative_APA/Human/data/RNAseq/sort/*.junc`
do
    echo $juncfile
    python filterNumChroms.py $juncfile
    echo $juncfile.fixed.junc >> ../data/DiffSplice_liftedJunc/Test_juncfiles.txt

done

python /project2/gilad/briana/leafcutter/clustering/leafcutter_cluster_regtools.py -j /project2/gilad/briana/Comparative_APA/data/DiffSplice_liftedJunc/Test_juncfiles.txt -m 10 -o humanJunc -l 500000
