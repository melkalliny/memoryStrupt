#!/bin/bash
tar -C /lscratch/$SLURM_JOB_ID -xf /usr/local/matlab-compiler/v94.tar.gz;
cd /home/elkallinymm/consolidationProject/biowulf_packaging/calculate; ./run_consolidation_calculateSimilarity_LowAndHigh_biowulf_swarm.sh /lscratch/$SLURM_JOB_ID/v94  NIH043 6 PROBE_START 1 &> /home/elkallinymm/consolidationProject/biowulf_packaging/calculate/scripts/12.log
