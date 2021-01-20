#!/bin/bash

source ~/activate_anaconda.sh
conda activate comp_threeprime_env

snakemake \
    -kp \
    --ri \
    -j 100 \
    -s /project2/gilad/briana/Comparative_APA/code/SnakefilePAS \
    --rerun-incomplete \
    --max-jobs-per-second 5 \
    --latency-wait 60 \
    --configfile Config_human.yaml \
    --cluster-config clusterPAS.json \
    -c "sbatch \
        --mem={cluster.mem} \
        --nodes={cluster.n} \
        --tasks-per-node={cluster.tasks} \
        --partition=broadwl \
        --job-name={cluster.name} \
        --output={cluster.logfile}" \
    $*
