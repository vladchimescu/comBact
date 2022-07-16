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

argv = commandArgs(trailingOnly = TRUE)

# argv[1] = data/ directory
# argv[2] = batch num
# argv[3] = strain

strain = paste0("screen", argv[3])

script_dir = "~/Documents/embl/gitlab/drug_comb_screen/dataprep/"
if(!dir.exists(script_dir)) script_dir = "dataprep/"

figdir = "~/Documents/embl/gitlab/drug_comb_screen/figures/wilcoxon/"
if(!dir.exists(figdir)) {
  figdir = paste0("figures/", strain, "/batch", argv[2], "/wilcoxon/")
  dir.create(figdir)
}

# directory for writing epsilon's (interaction scores)
eps_dir = "~/Documents/embl/screen_K/interactions/batch2/"
if(!dir.exists(eps_dir)) {
  eps_dir = file.path(argv[1], "interactions",
                      paste0("batch", argv[2], "/"))
}

wilc_dir = "~/Documents/embl/screen_K/wilcoxon/"
if(!dir.exists(wilc_dir)) {
  wilc_dir = file.path(argv[1], "wilcoxon/")
}

# source the script with core functions
source(paste0(script_dir, "code_screen_core_fun.R"))

#read the interaction file
for_wilc = read.table(paste0(eps_dir, "interactions.tsv"),
                      sep = "\t", header = T)

#remove self interactions
rows = apply(for_wilc[, c(5,12)], 1, function(i) length(unique(i)) > 1)
for_wilc = for_wilc[rows,]

#test = subset(for_wilc, combid_broad == "CTX_AMX")
for_wilc$recdon = for_wilc$combid_broad
for_wilc$donrec = paste0(for_wilc$donor, "_", for_wilc$Drug)

#df = subset(for_wilc, combid_broad == "CTX_AMX" | combid_broad == "AMX_CTX")

for_wilc$Drug <- as.character(for_wilc$Drug)
for_wilc$donor <- as.character(for_wilc$donor)
for_wilc$unique <- map2_chr(for_wilc$Drug, for_wilc$donor, ~ paste(sort(c(.x, .y)), collapse = "_"))

tmp = dplyr::group_by(for_wilc, unique) %>%
  dplyr::summarise(epsilons = list(eps))

eps_list = tmp$epsilons
names (eps_list) = tmp$unique

#do the same for expected fitness
tmp = dplyr::group_by(for_wilc, unique) %>%
  dplyr::summarise(exp_fit = list(multipl))

exp_fit_list = tmp$exp_fit
names (exp_fit_list) = tmp$unique


#======= Make permutations test for all bugs - Segregated distributions =========
#here removes null interactions
eps_list = eps_list[lapply(eps_list, length) > 0]
  
# *******= 10000 permutations per hit *************
#this is the df from which to make the sampling later
all_bug_data = as.numeric(unlist(eps_list))  

  
sign_cutoff = 0.1
pval_wilc = rep(NaN, length(eps_list))

# first identify strong drug-drug interactions based on the complete epsilon distribution
for(c in 1:length(eps_list)) {
    #see late how to do this
    comb = names(eps_list)[c]
    
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
    
    # compute the permutation test P-value (as in the Methods, see the paper)
    pval_wilc[c] = (sum(p_vals > sign_cutoff) + 1) / (length(p_vals) + 1 )
    
}

  
# determine synergistic / antagonistic interactions
# this assumes a symmetric distribution centered around 0
eps_th = 0.1

syn = lapply(eps_list, quantile, probs=0.25) < - eps_th
ant = lapply(eps_list, quantile, probs=0.75) > eps_th

# check if there are any interactions that are both
# synergetic and antagonistic at the same time
both = (syn & ant)
any(both)
which(both)  

#====================OBTAIN THE RESTRICTED EPSILON DISTRIBUTION=================
syn_th = 0.1
# antagonism cutoff
ant_th = 0.9

# restricted epsilon distribution for synergies
eps_syn = eps_list
eps_ant = eps_list
for(c in 1:length(eps_list))
{
  comb_data = eps_list[[c]]
  if(any(exp_fit_list[[c]] > syn_th)){
    eps_syn[[c]] = comb_data[exp_fit_list[[c]] > syn_th]
  } else eps_syn[[c]] = NA
  
  if(any(exp_fit_list[[c]] < ant_th)) {
    eps_ant[[c]] = comb_data[exp_fit_list[[c]] < ant_th]
  } else eps_ant[[c]] = NA
}


syn_rest = lapply(eps_syn, quantile, 
                  probs=0.25, na.rm=T) < - eps_th
syn_rest[is.na(syn_rest)] = FALSE
sum(syn_rest, na.rm = T)

ant_rest = lapply(eps_ant, quantile, 
                  probs=0.75, na.rm=T) > eps_th
ant_rest[is.na(ant_rest)] = FALSE
sum(ant_rest, na.rm = T)



#======================WRITE P VALUES========================================

pval.df = data.frame(pval = pval_wilc,
                     synergistic = syn,
                     antagonistic = ant,
                     synrest = syn_rest,
                     antrest = ant_rest,
                     comb = names(eps_list))
rownames(pval.df) = NULL

write.table(pval.df, 
            file = paste0(wilc_dir, "pvals-batch", argv[2], ".txt"),
            col.names = T, quote = F,
            row.names = F)


# detect strong drug-drug interactions based on the restricted epsilon distributions
# with only relevant drug concentrations

bliss.vs.expfit = dplyr::rename(for_wilc, comb = unique, expfit = multipl) %>%
  dplyr::inner_join(pval.df) %>%
  dplyr::mutate(type = ifelse(synergistic & !antagonistic, "synergy",
                              ifelse(antagonistic & !synergistic,
                                     "antagonism", "no interaction")))


# figure similar to Fig. Extendeed Data 3d
ext_data_3d = ggplot(bliss.vs.expfit, aes(x = expfit, y = eps, colour = factor(type))) +
  geom_point() + theme_classic() + 
  scale_color_manual(values = c("#C48B3F",
                                "#544E61",
                                "#85BAA1")) + 
  labs(x = "Expected fitness",
       y = "Bliss score",
       colour = "Interaction type")

ggsave(ext_data_3d, filename = paste0(figdir, "SF-3d.pdf"))


