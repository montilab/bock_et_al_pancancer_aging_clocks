#!/bin/bash -l

#$ -N tcga_mut_EDA
#$ -l h_rt=24:00:00
#$ -l mem_total=16G
#$ -P agedisease
#$ -l buyin
#$ -cwd
#$ -o logs/03_EDA_and_filtering-mutation.stdout
#$ -e logs/03_EDA_and_filtering-mutation.stderr
#$ -m bae
#$ -V

# Submit from: /restricted/projectnb/agedisease/projects/pancancer_aging/scripts/mutation/

# Load Modules
module load R/4.4.3
module load pandoc

# Log Job Info
echo "Job: $JOB_NAME | ID: $JOB_ID | User: $USER | Host: $HOSTNAME"
echo "Start: $(date) | Dir: $(pwd)"

# Run R Script
Rscript /restricted/projectnb/agedisease/projects/pancancer_aging/scripts/mutation/03_render_EDA_and_filtering-mutation.r

echo "Completed: $(date)"
