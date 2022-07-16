#!/bin/bash

#SBATCH -N 1                                                                    
#SBATCH -p htc                                                                  
#SBATCH -t 00:20:00                                                             
#SBATCH -n 1                                                                    
#SBATCH --mem=8G                                                               
#SBATCH --mail-user=vkim@embl.de                                                
#SBATCH --mail-type=FAIL,END 

module load R
module load UDUNITS/2.2.25-foss-2017b

DATADIR=$1
STRAIN=$2

echo "Data directory: " $DATADIR
echo "Strain: " $STRAIN

#Rscript qualitycheck/batcheffects.R $DATADIR $STRAIN
#echo "Done: batch effect plots"

Rscript analysis/permTest.R $DATADIR $STRAIN
echo "Done: permutation tests"



