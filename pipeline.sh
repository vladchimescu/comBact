#!/bin/bash

#SBATCH -N 1                                                                    
#SBATCH -p htc                                                                  
#SBATCH -t 00:45:00                                                             
#SBATCH -n 1                                                                    
#SBATCH --mem=20G                                                               
#SBATCH --mail-user=vkim@embl.de                                                
#SBATCH --mail-type=FAIL,END 

module load R
# module load UDUNITS/2.2.25-foss-2017b

DATADIR=$1
BATCH=$2
STRAIN=$3

echo "Data directory: " $DATADIR
echo "Processing batch: " $BATCH
echo "Strain: " $STRAIN

Rscript dataprep/annotate.R $DATADIR $BATCH $STRAIN
echo "Done: preprocessing data"

Rscript qualitycheck/qc_check.R $DATADIR $BATCH $STRAIN
echo "Done: quality control plots"

Rscript analysis/code_inter_scores_new.R $DATADIR $BATCH $STRAIN
echo "Done: interaction score computation"




