# run Wilcoxon test on all batches
argv = commandArgs(trailingOnly = TRUE)
# argv[1] = data/ directory
# argv[2] = strain
# argv[3] = chunk number
# argv[4] = first index of eps_list of this job array
# argv[5] = last index of eps_list of this job array

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

set.seed(as.numeric(argv[3]))
strain = paste0("screen", argv[2])

params = rjson::fromJSON(file = "params.json")
eps_th = 0.1
# strongsyn = params$strongsyn
# strongant = params$strongant


script_dir = "~/Documents/embl/gitlab/drug_comb_screen/dataprep/"
if(!dir.exists(script_dir)) script_dir = "dataprep/"

figdir = "~/Documents/embl/gitlab/drug_comb_screen/figures/wilcoxon/"
if(!dir.exists(figdir)) {
  figdir = paste0("figures/", strain, "/wilcoxon/")
  dir.create(figdir)
}

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


batches = grep("batch", list.dirs(eps_dir, recursive = F, 
          full.names = F), value = T)
eps = lapply(batches, function(b) {
  eps = read.table(paste0(eps_dir, b, "/interactions.tsv"),
                    sep = "\t", header = T)
  eps = prep.interactions(eps)
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


# *******= 10000 permutations per hit *************
#this is the df from which to make the sampling later
all_bug_data = as.numeric(unlist(eps_list))  

sign_cutoff = 0.1

do_perm_test <- function(c) {
  comb = names(eps_list)[c]
  message(paste("Testing:", comb))
  
  comb_data = eps_list[[c]]
  comb_data = comb_data[!is.na(comb_data)]
  
  #Complete range
  comb_data_n = length(comb_data)
  
  #here minimum number of wells (values) to consider
  if(comb_data_n>3)
  {
    random_samples = as.data.frame(replicate(10000,sample(all_bug_data,comb_data_n)))
    permut_test = as.data.frame(sapply(random_samples,FUN=wilcox.test,y=comb_data))
    p_vals = as.numeric(as.matrix(permut_test[3,]))
  } else
  {p_vals = c(NA)}
  
  pval_wilc = (sum(p_vals > sign_cutoff) + 1) / (length(p_vals) + 1 )
  
  data.frame(comb = comb, pval = pval_wilc)
}

pvals = lapply(argv[4]:argv[5], do_perm_test)

pval.df = do.call(rbind, pvals)

write.table(pval.df, file = paste0(wilc_dir, "pvals-", argv[3], ".txt"),
            col.names = T, quote = F,
            row.names = F)



# #====================OBTAIN THE RESTRICTED EPSILON DISTRIBUTION=================
syn_th = 0.2
# antagonism cutoff
ant_th = 0.8

# restricted epsilon distribution for synergies
eps_syn = eps_list
eps_ant = eps_list
for(c in 1:length(eps_list))
{
  comb_data = eps_list[[c]]
  if(any(exp_fit_list[[c]] >= syn_th)){
    eps_syn[[c]] = comb_data[exp_fit_list[[c]] > syn_th]
  } else eps_syn[[c]] = NA

  if(any(exp_fit_list[[c]] < ant_th)) {
    eps_ant[[c]] = comb_data[exp_fit_list[[c]] <= ant_th]
  } else eps_ant[[c]] = NA
}

# now find restricted synergies
eps_list = eps_syn
syn_rest = lapply(argv[4]:argv[5], do_perm_test)
syn_rest_df = plyr::ldply(syn_rest)
write.table(syn_rest_df, file = paste0(wilc_dir, "synrest-pvals-", argv[3], ".txt"),
            col.names = T, quote = F,
            row.names = F)

# save restricted synergies

# find restricted antagonisms
eps_list = eps_ant
ant_rest = lapply(argv[4]:argv[5], do_perm_test)
ant_rest_df = plyr::ldply(ant_rest)
write.table(ant_rest_df, file = paste0(wilc_dir, "antrest-pvals-", argv[3], ".txt"),
            col.names = T, quote = F,
            row.names = F)



