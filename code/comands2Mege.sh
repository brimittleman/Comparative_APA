#!/bin/bash


#SBATCH --job-name=mergethreeprime
#SBATCH --mail-type=END
#SBATCH --error=merge.err
#SBATCH --partition=broadwl
#SBATCH --mem=36G


cat /project2/gilad/briana/191121_NB501189_0663_AHFY73BGXC-YG-BM-24S-ln1-Adjusted/FastQ/YG-BM-24S-F2-18498-N_S14_R1_001.fastq /project2/gilad/briana/191123_NB501189_0665_AHG2NWBGXC-YG-BM-24S-ln2-Adjusted/FastQ/YG-BM-24S-F2-18498-N_S14_R1_001.fastq /project2/gilad/briana/191127_NB501189_0667_AHFY5TBGXC-YG-BM-24S-ln3/FastQ/*18498-N* > /project2/gilad/briana/Comparative_APA/Chimp/data/fastq/chimp_combined_4973_N.fastq
cat /project2/gilad/briana/191121_NB501189_0663_AHFY73BGXC-YG-BM-24S-ln1-Adjusted/FastQ/YG-BM-24S-E2-18498-T_S13_R1_001.fastq /project2/gilad/briana/191123_NB501189_0665_AHG2NWBGXC-YG-BM-24S-ln2-Adjusted/FastQ/YG-BM-24S-E2-18498-T_S13_R1_001.fastq /project2/gilad/briana/191127_NB501189_0667_AHFY5TBGXC-YG-BM-24S-ln3/FastQ/*18498-T* > /project2/gilad/briana/Comparative_APA/Chimp/data/fastq/chimp_combined_4973_T.fastq
cat /project2/gilad/briana/191121_NB501189_0663_AHFY73BGXC-YG-BM-24S-ln1-Adjusted/FastQ/YG-BM-24S-H2-18499-N_S16_R1_001.fastq /project2/gilad/briana/191123_NB501189_0665_AHG2NWBGXC-YG-BM-24S-ln2-Adjusted/FastQ/YG-BM-24S-H2-18499-N_S16_R1_001.fastq /project2/gilad/briana/191127_NB501189_0667_AHFY5TBGXC-YG-BM-24S-ln3/FastQ/*18499-N* > /project2/gilad/briana/Comparative_APA/Human/data/fastq/human_combined_18499_N.fastq
cat /project2/gilad/briana/191121_NB501189_0663_AHFY73BGXC-YG-BM-24S-ln1-Adjusted/FastQ/YG-BM-24S-G2-18499-T_S15_R1_001.fastq /project2/gilad/briana/191123_NB501189_0665_AHG2NWBGXC-YG-BM-24S-ln2-Adjusted/FastQ/YG-BM-24S-G2-18499-T_S15_R1_001.fastq /project2/gilad/briana/191127_NB501189_0667_AHFY5TBGXC-YG-BM-24S-ln3/FastQ/*18499-T* > /project2/gilad/briana/Comparative_APA/Human/data/fastq/human_combined_18499_T.fastq
cat /project2/gilad/briana/191121_NB501189_0663_AHFY73BGXC-YG-BM-24S-ln1-Adjusted/FastQ/YG-BM-24S-B3-18502-N_S18_R1_001.fastq /project2/gilad/briana/191123_NB501189_0665_AHG2NWBGXC-YG-BM-24S-ln2-Adjusted/FastQ/YG-BM-24S-B3-18502-N_S18_R1_001.fastq /project2/gilad/briana/191127_NB501189_0667_AHFY5TBGXC-YG-BM-24S-ln3/FastQ/*18502-N* > /project2/gilad/briana/Comparative_APA/Human/data/fastq/human_combined_18502_N.fastq
cat /project2/gilad/briana/191121_NB501189_0663_AHFY73BGXC-YG-BM-24S-ln1-Adjusted/FastQ/YG-BM-24S-A3-18502-T_S17_R1_001.fastq /project2/gilad/briana/191123_NB501189_0665_AHG2NWBGXC-YG-BM-24S-ln2-Adjusted/FastQ/*18502-T* /project2/gilad/briana/191127_NB501189_0667_AHFY5TBGXC-YG-BM-24S-ln3/FastQ/*18502-T* > /project2/gilad/briana/Comparative_APA/Human/data/fastq/human_combined_18502_T.fastq
cat /project2/gilad/briana/191121_NB501189_0663_AHFY73BGXC-YG-BM-24S-ln1-Adjusted/FastQ/YG-BM-24S-D3-18504-N_S20_R1_001.fastq /project2/gilad/briana/191123_NB501189_0665_AHG2NWBGXC-YG-BM-24S-ln2-Adjusted/FastQ/*18504-N* /project2/gilad/briana/191127_NB501189_0667_AHFY5TBGXC-YG-BM-24S-ln3/FastQ/*18504-N* > /project2/gilad/briana/Comparative_APA/Human/data/fastq/human_combined_18504_N.fastq
cat /project2/gilad/briana/191121_NB501189_0663_AHFY73BGXC-YG-BM-24S-ln1-Adjusted/FastQ/YG-BM-24S-C3-18504-T_S19_R1_001.fastq /project2/gilad/briana/191123_NB501189_0665_AHG2NWBGXC-YG-BM-24S-ln2-Adjusted/FastQ/*18504-T* /project2/gilad/briana/191127_NB501189_0667_AHFY5TBGXC-YG-BM-24S-ln3/FastQ/*18504-T* > /project2/gilad/briana/Comparative_APA/Human/data/fastq/human_combined_18504_T.fastq
cat /project2/gilad/briana/191121_NB501189_0663_AHFY73BGXC-YG-BM-24S-ln1-Adjusted/FastQ/YG-BM-24S-F3-18510-N_S22_R1_001.fastq /project2/gilad/briana/191123_NB501189_0665_AHG2NWBGXC-YG-BM-24S-ln2-Adjusted/FastQ/*18510-N* /project2/gilad/briana/191127_NB501189_0667_AHFY5TBGXC-YG-BM-24S-ln3/FastQ/*18510-N* > /project2/gilad/briana/Comparative_APA/Human/data/fastq/human_combined_18510_N.fastq
cat /project2/gilad/briana/191121_NB501189_0663_AHFY73BGXC-YG-BM-24S-ln1-Adjusted/FastQ/YG-BM-24S-E3-18510-T_S21_R1_001.fastq /project2/gilad/briana/191123_NB501189_0665_AHG2NWBGXC-YG-BM-24S-ln2-Adjusted/FastQ/*18510-T* /project2/gilad/briana/191127_NB501189_0667_AHFY5TBGXC-YG-BM-24S-ln3/FastQ/*18504-T* > /project2/gilad/briana/Comparative_APA/Human/data/fastq/human_combined_18510_T.fastq
cat /project2/gilad/briana/191121_NB501189_0663_AHFY73BGXC-YG-BM-24S-ln1-Adjusted/FastQ/YG-BM-24S-H3-18523-N_S24_R1_001.fastq /project2/gilad/briana/191123_NB501189_0665_AHG2NWBGXC-YG-BM-24S-ln2-Adjusted/FastQ/YG-BM-24S-H3-18523-N_S24_R1_001.fastq /project2/gilad/briana/191127_NB501189_0667_AHFY5TBGXC-YG-BM-24S-ln3/FastQ/*18523-N* > /project2/gilad/briana/Comparative_APA/Human/data/fastq/human_combined_18523_N.fastq
cat /project2/gilad/briana/191121_NB501189_0663_AHFY73BGXC-YG-BM-24S-ln1-Adjusted/FastQ/YG-BM-24S-G3-18523-T_S23_R1_001.fastq /project2/gilad/briana/191123_NB501189_0665_AHG2NWBGXC-YG-BM-24S-ln2-Adjusted/FastQ/YG-BM-24S-G3-18523-T_S23_R1_001.fastq /project2/gilad/briana/191127_NB501189_0667_AHFY5TBGXC-YG-BM-24S-ln3/FastQ/*18523-T* > /project2/gilad/briana/Comparative_APA/Human/data/fastq/human_combined_18523_T.fastq
cat /project2/gilad/briana/191121_NB501189_0663_AHFY73BGXC-YG-BM-24S-ln1-Adjusted/FastQ/YG-BM-24S-D2-18358-N_S12_R1_001.fastq /project2/gilad/briana/191123_NB501189_0665_AHG2NWBGXC-YG-BM-24S-ln2-Adjusted/FastQ/YG-BM-24S-D2-18358-N_S12_R1_001.fastq /project2/gilad/briana/191127_NB501189_0667_AHFY5TBGXC-YG-BM-24S-ln3/FastQ/*18358-N* > /project2/gilad/briana/Comparative_APA/Chimp/data/fastq/chimp_combined_18358_N.fastq
cat /project2/gilad/briana/191121_NB501189_0663_AHFY73BGXC-YG-BM-24S-ln1-Adjusted/FastQ/YG-BM-24S-C2-18358-T_S11_R1_001.fastq /project2/gilad/briana/191123_NB501189_0665_AHG2NWBGXC-YG-BM-24S-ln2-Adjusted/FastQ/YG-BM-24S-C2-18358-T_S11_R1_001.fastq /project2/gilad/briana/191127_NB501189_0667_AHFY5TBGXC-YG-BM-24S-ln3/FastQ/*18358-T* > /project2/gilad/briana/Comparative_APA/Chimp/data/fastq/chimp_combined_18358_T.fastq
cat /project2/gilad/briana/191121_NB501189_0663_AHFY73BGXC-YG-BM-24S-ln1-Adjusted/FastQ/YG-BM-24S-F1-3622-N_S6_R1_001.fastq /project2/gilad/briana/191123_NB501189_0665_AHG2NWBGXC-YG-BM-24S-ln2-Adjusted/FastQ/YG-BM-24S-F1-3622-N_S6_R1_001.fastq /project2/gilad/briana/191127_NB501189_0667_AHFY5TBGXC-YG-BM-24S-ln3/FastQ/*3622-N* > /project2/gilad/briana/Comparative_APA/Chimp/data/fastq/chimp_combined_3622_N.fastq
cat /project2/gilad/briana/191121_NB501189_0663_AHFY73BGXC-YG-BM-24S-ln1-Adjusted/FastQ/YG-BM-24S-E1-3622-T_S5_R1_001.fastq /project2/gilad/briana/191123_NB501189_0665_AHG2NWBGXC-YG-BM-24S-ln2-Adjusted/FastQ/YG-BM-24S-E1-3622-T_S5_R1_001.fastq /project2/gilad/briana/191127_NB501189_0667_AHFY5TBGXC-YG-BM-24S-ln3/FastQ/*3622-T* > /project2/gilad/briana/Comparative_APA/Chimp/data/fastq/chimp_combined_3622_T.fastq
cat /project2/gilad/briana/191121_NB501189_0663_AHFY73BGXC-YG-BM-24S-ln1-Adjusted/FastQ/YG-BM-24S-H1-3659-N_S8_R1_001.fastq /project2/gilad/briana/191123_NB501189_0665_AHG2NWBGXC-YG-BM-24S-ln2-Adjusted/FastQ/YG-BM-24S-H1-3659-N_S8_R1_001.fastq /project2/gilad/briana/191127_NB501189_0667_AHFY5TBGXC-YG-BM-24S-ln3/FastQ/*3659-N* > /project2/gilad/briana/Comparative_APA/Chimp/data/fastq/chimp_combined_3659_N.fastq
cat /project2/gilad/briana/191121_NB501189_0663_AHFY73BGXC-YG-BM-24S-ln1-Adjusted/FastQ/YG-BM-24S-G1-3659-T_S7_R1_001.fastq /project2/gilad/briana/191123_NB501189_0665_AHG2NWBGXC-YG-BM-24S-ln2-Adjusted/FastQ/YG-BM-24S-G1-3659-T_S7_R1_001.fastq /project2/gilad/briana/191127_NB501189_0667_AHFY5TBGXC-YG-BM-24S-ln3/FastQ/*3659-T* > /project2/gilad/briana/Comparative_APA/Chimp/data/fastq/chimp_combined_3659_T.fastq
cat /project2/gilad/briana/191121_NB501189_0663_AHFY73BGXC-YG-BM-24S-ln1-Adjusted/FastQ/YG-BM-24S-B2-4973-N_S10_R1_001.fastq /project2/gilad/briana/191123_NB501189_0665_AHG2NWBGXC-YG-BM-24S-ln2-Adjusted/FastQ/YG-BM-24S-B2-4973-N_S10_R1_001.fastq /project2/gilad/briana/191127_NB501189_0667_AHFY5TBGXC-YG-BM-24S-ln3/FastQ/*4973-N* > /project2/gilad/briana/Comparative_APA/Human/data/fastq/human_combined_18498_N.fastq
cat /project2/gilad/briana/191121_NB501189_0663_AHFY73BGXC-YG-BM-24S-ln1-Adjusted/FastQ/YG-BM-24S-A2-4973-T_S9_R1_001.fastq /project2/gilad/briana/191123_NB501189_0665_AHG2NWBGXC-YG-BM-24S-ln2-Adjusted/FastQ/YG-BM-24S-A2-4973-T_S9_R1_001.fastq /project2/gilad/briana/191127_NB501189_0667_AHFY5TBGXC-YG-BM-24S-ln3/FastQ/*4973-T* > /project2/gilad/briana/Comparative_APA/Human/data/fastq/human_combined_18498_T.fastq
cat /project2/gilad/briana/191121_NB501189_0663_AHFY73BGXC-YG-BM-24S-ln1-Adjusted/FastQ/YG-BM-24S-B1-pt30-N_S2_R1_001.fastq /project2/gilad/briana/191123_NB501189_0665_AHG2NWBGXC-YG-BM-24S-ln2-Adjusted/FastQ/YG-BM-24S-B1-pt30-N_S2_R1_001.fastq /project2/gilad/briana/191127_NB501189_0667_AHFY5TBGXC-YG-BM-24S-ln3/FastQ/*pt30-N* > /project2/gilad/briana/Comparative_APA/Chimp/data/fastq/chimp_combined_pt30_N.fastq
cat /project2/gilad/briana/191121_NB501189_0663_AHFY73BGXC-YG-BM-24S-ln1-Adjusted/FastQ/YG-BM-24S-A1-pt30-T_S1_R1_001.fastq /project2/gilad/briana/191123_NB501189_0665_AHG2NWBGXC-YG-BM-24S-ln2-Adjusted/FastQ/YG-BM-24S-A1-pt30-T_S1_R1_001.fastq /project2/gilad/briana/191127_NB501189_0667_AHFY5TBGXC-YG-BM-24S-ln3/FastQ/*pt30-T* > /project2/gilad/briana/Comparative_APA/Chimp/data/fastq/chimp_combined_pt30_T.fastq
cat /project2/gilad/briana/191121_NB501189_0663_AHFY73BGXC-YG-BM-24S-ln1-Adjusted/FastQ/YG-BM-24S-B1-pt30-N_S2_R1_001.fastq /project2/gilad/briana/191123_NB501189_0665_AHG2NWBGXC-YG-BM-24S-ln2-Adjusted/FastQ/YG-BM-24S-D1-pt91-N_S4_R1_001.fastq /project2/gilad/briana/191127_NB501189_0667_AHFY5TBGXC-YG-BM-24S-ln3/FastQ/*pt91-N* > /project2/gilad/briana/Comparative_APA/Chimp/data/fastq/chimp_combined_pt91_N.fastq
cat /project2/gilad/briana/191121_NB501189_0663_AHFY73BGXC-YG-BM-24S-ln1-Adjusted/FastQ/YG-BM-24S-C1-pt91-T_S3_R1_001.fastq /project2/gilad/briana/191123_NB501189_0665_AHG2NWBGXC-YG-BM-24S-ln2-Adjusted/FastQ/YG-BM-24S-C1-pt91-T_S3_R1_001.fastq /project2/gilad/briana/191127_NB501189_0667_AHFY5TBGXC-YG-BM-24S-ln3/FastQ/*pt91-T* > /project2/gilad/briana/Comparative_APA/Chimp/data/fastq/chimp_combined_pt91_T.fastq