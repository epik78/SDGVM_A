#!/bin/bash
#$ -l h_rt=03:00:00
#$ -l rmem=12G

sites=114

module load apps/matlab/2017b
matlab -nosplash -nodesktop -r "run_flux_sites($sites)"

for l in $(seq 1 $sites)
do
    
    qsub site_$l.sge

done