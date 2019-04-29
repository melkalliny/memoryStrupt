#!/bin/bash
tar -C /lscratch/$SLURM_JOB_ID -xf /usr/local/matlab-compiler/v94.tar.gz;
cd /home/elkallinymm/consolidationProject/biowulf_packaging/calculate; ./run_consolidation_calculateSimilarity_LowAndHigh_biowulf_swarm.sh /lscratch/$SLURM_JOB_ID/v94  NIH050 4 RESPONSE 3 &> /home/elkallinymm/consolidationProject/biowulf_packaging/calculate/scripts/28.log
