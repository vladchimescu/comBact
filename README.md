# Combinatorial Drug Screen
*Elisabetta Cacace, Vladislav Kim (collaborator)*

Code for processsing and analysing the combinatorial screen data. The repository tree has the following structure: 
+ `dataprep` folder: R scripts for pre-processing of the raw screen data
+ `qualitycheck` directory: scripts for quality control
+ `analysis`: high-level analysis of the screen data, computation of interaction scores, p-values and network analysis


On this branch the scripts form the data analysis pipeline for the drug combinatorial screen. In order to run the pipeline in the SLURM cluster environment:

```sbatch pipeline.sh datadir batchnumber strain```

For example, provided that the data is in  and you wannt to process the data from batch 2 (strain K), submit the job as follows:

```sbatch pipeline.sh data/screen_K 2 K```

Furthermore, the data can be analyzed at the strain level, so that the significant interactions are identified separately for each strain. First run permutation test in parallel:

```sbatch permtest.sh datadir strain```

For example if one wants to process all of the data for strain K:

```sbatch permtest.sh data/screen_K K```

Next run `plot.sh` script:

```
sbatch plot.sh datadir strain
# strain K
sbatch plot.sh data/screen_K K
```

For pair-wise strain comparison run the following script:
```
sbatch compare_strains.sh data/screen_K data/screen_M
```
