#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=20G,h_vmem=20G,h_fsize=100G
#$ -pe local 4
#$ -N all_samples
#$ -o logs/all_samples.$TASK_ID.txt
#$ -e logs/all_samples.$TASK_ID.txt
#$ -m e
#$ -t 1-10
#$ -tc 10

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SGE_TASK_ID}"

## Load the R module (absent since the JHPCE upgrade to CentOS v7)
module load conda_R

## List current modules for reproducibility
module list

## Edit with your job command
#Rscript -e "options(width = 120); sessioninfo::session_info()"
Rscript all_samples.R --input ${SGE_TASK_ID}
echo "**** Job ends ****"
date

## This script was made using sgejobs version 0.99.1
## available from http://research.libd.org/sgejobs/