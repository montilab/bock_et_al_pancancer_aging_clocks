#!/bin/bash -l
#$ -P agedisease
#$ -N download_tcga
#$ -j y

module load R/4.3.1

PROJECT=$1

Rscript 01_DownloadData-methylation.R ${PROJECT}



