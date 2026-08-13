#!/bin/bash -l
#$ -N GISTIC_regressions
#$ -l h_rt=24:00:00
#$ -l mem_total=16G
#$ -P agedisease
#$ -l buyin
#$ -cwd
#$ -o logs/04_regressions-scna.stdout
#$ -e logs/04_regressions-scna.stderr
#$ -m bae
#$ -V

module load R/4.4.3
module load pandoc

echo "Job: $JOB_NAME | ID: $JOB_ID | User: $USER | Host: $HOSTNAME"
echo "Start: $(date) | Dir: $(pwd)"

Rscript /restricted/projectnb/agedisease/projects/pancancer_aging/scripts/scna_pn/04_render_regressions-scna.r

echo "Completed: $(date)"
