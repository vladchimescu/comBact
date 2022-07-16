# run Wilcoxon test on all batches
argv = commandArgs(trailingOnly = TRUE)
# argv[1] = data/ directory
# argv[2] = strain

require(graphics)
require(stats)
require(Hmisc)
require(zoo)
require(gplots)
require(LSD)
require(corrplot)
library(igraph)
require(plotrix)
require (purrr)

strain = paste0("screen", argv[2])

script_dir = "~/Documents/embl/gitlab/drug_comb_screen/dataprep/"
if(!dir.exists(script_dir)) script_dir = "dataprep/"

# directory for writing epsilon's (interaction scores)
eps_dir = "~/Documents/embl/screen_K/interactions/"
if(!dir.exists(eps_dir)) {
  eps_dir = file.path(argv[1], "interactions/")
}


wilc_dir = "~/Documents/embl/screen_K/wilcoxon/"
if(!dir.exists(wilc_dir)) {
  wilc_dir = file.path(argv[1], "wilcoxon/")
}

# source the script with core functions
source(paste0(script_dir, "code_screen_core_fun.R"))


prep.interactions <- function(for_wilc) {
  #remove self interactions
  rows = apply(for_wilc[, c("Drug", "donor")], 1, function(i) length(unique(i)) > 1)
  for_wilc = for_wilc[rows,]
  
  for_wilc$recdon = for_wilc$combid_broad
  for_wilc$donrec = paste0(for_wilc$donor, "_", for_wilc$Drug)
  
  
  for_wilc$Drug <- as.character(for_wilc$Drug)
  for_wilc$donor <- as.character(for_wilc$donor)
  for_wilc$unique <- map2_chr(for_wilc$Drug, for_wilc$donor, ~ paste(sort(c(.x, .y)), collapse = "_"))
  
  for_wilc = subset(for_wilc, !Drug %in% c("ONE", "BUG", "NOT"))
  for_wilc = subset(for_wilc, !donor %in% c("ONE", "BUG", "NOT"))
  for_wilc
}

get.eps.list <- function(for_wilc) {
  tmp = dplyr::group_by(for_wilc, unique) %>%
    dplyr::summarise(epsilons = list(eps))
  
  eps_list = tmp$epsilons
  names (eps_list) = tmp$unique
  eps_list
}

get.exp.fitness <- function(for_wilc) {
  tmp = dplyr::group_by(for_wilc, unique) %>%
    dplyr::summarise(exp_fit = list(multipl))
  
  exp_fit_list = tmp$exp_fit
  names (exp_fit_list) = tmp$unique
  exp_fit_list
}

if(argv[2] == 'M') {
  drugleg = read.delim2(file.path(argv[1], 'legend_figs_V4.txt'))
  drugs_main = dplyr::filter(drugleg, main_or_HTD_screen == 'main')$code
  drugs_human = dplyr::filter(drugleg, main_or_HTD_screen == 'htd')$code
}

batches = grep("batch", list.dirs(eps_dir, recursive = F, 
                                  full.names = F), value = T)
eps = lapply(batches, function(b) {
  eps = read.table(paste0(eps_dir, b, "/interactions.tsv"),
                   sep = "\t", header = T)
  eps = prep.interactions(eps)
  if(argv[2] == 'M') {
    eps = dplyr::filter(eps, Drug %in% drugs_main & donor %in% drugs_main)
  }
  eps
})

eps_batches = lapply(eps, get.eps.list)
exp_fit_batches = lapply(eps, get.exp.fitness)


# merge eps lists by (unique) drug combination
keys <- unique(unlist(lapply(eps_batches, names)))
eps_list = setNames(do.call(mapply, c(FUN=c, lapply(eps_batches, `[`, keys))), keys)

# merge expected fitness lists by (unique) drug combination
keys <- unique(unlist(lapply(exp_fit_batches, names)))
exp_fit_list = setNames(do.call(mapply, c(FUN=c, lapply(exp_fit_batches, `[`, keys))), keys)

stopifnot(length(exp_fit_list) == length(eps_list))


# submit 200 jobs 
njobs = 200
inds = split(1:length(eps_list), sort(rep_len(1:njobs, length(eps_list))))

for (i in 1:njobs) {
fname = paste0(strain, "-job", i, ".sh")
cat("#!/bin/bash
#SBATCH -N 1                                                                    
#SBATCH -p htc                                                                  
#SBATCH -t 00:15:00                                                             
#SBATCH -n 1                                                                    
#SBATCH --mem=8G                                                               
#SBATCH --mail-user=vkim@embl.de                                                
#SBATCH --mail-type=FAIL

module load R
module load UDUNITS/2.2.25-foss-2017b

Rscript analysis/combinedWilcoxon.R", argv[1], argv[2], i, min(inds[[i]]), max(inds[[i]]), file = fname)
  
system(paste("sbatch", fname))
}

# for DSM20231 also run permutation test separately for human drugs
if(argv[2] == 'M') {
  
  batches = grep("batch", list.dirs(eps_dir, recursive = F, 
                                    full.names = F), value = T)
  eps = lapply(batches, function(b) {
    eps = read.table(paste0(eps_dir, b, "/interactions.tsv"),
                     sep = "\t", header = T)
    eps = prep.interactions(eps)
    if(argv[2] == 'M') {
      eps = dplyr::filter(eps, Drug %in% drugs_human | donor %in% drugs_human)
    }
    eps
  })
  
  eps_batches = lapply(eps, get.eps.list)
  exp_fit_batches = lapply(eps, get.exp.fitness)
  
  
  # merge eps lists by (unique) drug combination
  keys <- unique(unlist(lapply(eps_batches, names)))
  eps_list = setNames(do.call(mapply, c(FUN=c, lapply(eps_batches, `[`, keys))), keys)
  
  # merge expected fitness lists by (unique) drug combination
  keys <- unique(unlist(lapply(exp_fit_batches, names)))
  exp_fit_list = setNames(do.call(mapply, c(FUN=c, lapply(exp_fit_batches, `[`, keys))), keys)
  
  stopifnot(length(exp_fit_list) == length(eps_list))
  
  
  # submit 200 jobs 
  njobs = 200
  inds = split(1:length(eps_list), sort(rep_len(1:njobs, length(eps_list))))
  
  for (i in 1:njobs) {
    fname = paste0("htd-",strain, "-job", i, ".sh")
    cat("#!/bin/bash
#SBATCH -N 1                                                                    
#SBATCH -p htc                                                                  
#SBATCH -t 00:15:00                                                             
#SBATCH -n 1                                                                    
#SBATCH --mem=8G                                                               
#SBATCH --mail-user=vkim@embl.de                                                
#SBATCH --mail-type=FAIL

module load R
module load UDUNITS/2.2.25-foss-2017b

Rscript analysis/combinedWilcoxon_human.R", argv[1], argv[2], i, min(inds[[i]]), max(inds[[i]]), file = fname)
    
    system(paste("sbatch", fname))
  }
  
  
}


