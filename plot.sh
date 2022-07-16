#!/bin/bash

#SBATCH -N 1                                                                    
#SBATCH -p htc                                                                  
#SBATCH -t 01:00:00                                                             
#SBATCH -n 1                                                                    
#SBATCH --mem=16G                                                               
#SBATCH --mail-user=vkim@embl.de                                                
#SBATCH --mail-type=FAIL,END 

module load R
module load UDUNITS/2.2.25-foss-2017b

DATADIR=$1
STRAIN=$2

echo "Data directory: " $DATADIR
echo "Strain: " $STRAIN

Rscript analysis/combinePvals.R $DATADIR $STRAIN
echo "Done: combined the data from all batches"

Rscript analysis/multipletesting.R $DATADIR $STRAIN
echo "Done: multiple testing"

#Rscript analysis/binomial_test.R $DATADIR $STRAIN
#echo "Done: binomial test"

Rscript analysis/generate_comb_OD.R $DATADIR $STRAIN
echo "Done: Generated combination OD data for the Shiny app"

#Rscript analysis/lfdr-test.R $DATADIR $STRAIN
#echo "Done: Local FDR"

Rscript analysis/checkerboard.R $DATADIR $STRAIN
echo "Done: checkerboard plots"

#Rscript analysis/checkerboard-binomial.R $DATADIR $STRAIN
#echo "Done: checkerboard plots for binomial test results"

Rscript qualitycheck/detectOutliers.R $DATADIR $STRAIN
echo "Done: plotted Bliss scores in a single matrix"

#Rscript qualitycheck/threshold-checkerboard.R $DATADIR $STRAIN
#echo "Done: plotted checkerboards for interactions close to the threshold"
#
#Rscript qualitycheck/thresh-blissmat.R $DATADIR $STRAIN
#echo "Done: plotted Bliss score matrix for combs close to the threshold"

Rscript shiny/generate_single_drug_OD.R $DATADIR $STRAIN
echo "Done: generated single-drug OD data for the Shiny app"

Rscript shiny/generateAllbliss.R $DATADIR $STRAIN
echo "Done: aggregated all Bliss scores in one .rds file"

# clean up the Wilcoxon directory
#find $DATADIR/wilcoxon/ -name 'pvals-[0-9]*' -delete

