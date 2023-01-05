# Combinatorial Drug Screen

Code for processing and analysing the combinatorial screen data presented in the manuscript "High-throughput profiling of drug interactions in Gram-positive bacteria” (https://www.biorxiv.org/content/10.1101/2022.12.23.521747v1). 

All drug-drug interactions (with their effect size and adjusted p-value) detected with this pipeline are listed in Supplementary Table 2 (antibiotic screen conducted in *S. aureus* Newman and DSM 20231, *B. subtilis* 168 and *S. pneumoniae* D39V) and Supplementary Table 6 (non-antibiotic drug screen in *S. aureus* DSM 20231). https://www.biorxiv.org/content/10.1101/2022.12.23.521747v1.supplementary-material

The repository tree has the following structure:
+ `data`: unprocessed optical density measurements
+ `dataprep` folder: R scripts for pre-processing of the raw screen data
+ `qualitycheck` directory: scripts for quality control
+ `analysis`: high-level analysis of the screen data, computation of interaction scores, p-values and network analysis
+ `phylogeny`: code and data for the phylogeny analysis

Note that here we use the following codes for the 4 Gram-positive bacterial strains:
Code | Strain
---|---|
K| *S. aureus* Newman
M| *S. aureus* DSM 20231
D| *B. subtilis* 168
P| *S. pneumoniae* D39V

The main branch contains the data analysis pipeline for the drug combinatorial screen in Gram-positive bacterial species. In order to run the pipeline in the SLURM cluster environment:

```sbatch pipeline.sh datadir batchnumber strain```

For example, provided that the data is in  and you wannt to process the data from batch 2 (strain K), submit the job as follows:

```sbatch pipeline.sh data/screen_K 2 K```

Furthermore, the data can be analyzed at the strain level, so that the significant interactions are identified separately for each strain. First run permutation test in parallel:

```sbatch permtest.sh datadir strain```

For example if one wants to process all of the data for strain K:

```sbatch permtest.sh data/screen_K K```

Next run `plot.sh` script to generate visualizations of individual drug combinations:

```
sbatch plot.sh datadir strain
# strain K
sbatch plot.sh data/screen_K K
```
