#!/bin/bash
#$ -l h_rt=0:30:00
#$ -l rmem=12G
#$ -cwd

sites=(17 31 46 52 57 109)

noruns=3

module load apps/matlab/2017b
matlab -nosplash -nodesktop -r "emu_runs([${sites[*]}])"

#for l in $(seq 1 $noruns)
#do
#    qsub disrun_$l.sge
#done

