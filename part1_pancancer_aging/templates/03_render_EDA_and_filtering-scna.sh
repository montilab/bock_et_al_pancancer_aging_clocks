#!/bin/bash -l
#$ -N GISTIC_EDA
#$ -l h_rt=24:00:00
#$ -l mem_total=16G
#$ -P agedisease
#$ -l buyin
#$ -cwd
#$ -o logs/03_EDA_and_filtering-scna.stdout
#$ -e logs/03_EDA_and_filtering-scna.stderr
#$ -m bae
#$ -V

module load R/4.4.3
module load pandoc

echo "Job: $JOB_NAME | ID: $JOB_ID | User: $USER | Host: $HOSTNAME"
echo "Start: $(date) | Dir: $(pwd)"

Rscript /restricted/projectnb/agedisease/projects/pancancer_aging/scripts/scna_pn/03_render_EDA_and_filtering-scna.r

echo "Completed: $(date)"
