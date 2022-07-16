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

write.table(pval.df, file = paste0(wilc_dir, "htd-pvals-", argv[3], ".txt"),
            col.names = T, quote = F,
            row.names = F)

# commented this out for now
# # determine synergistic / antagonistic interactions
# # this assumes a symmetric distribution centered around 0
# eps_th = 0.1
# 
# syn = lapply(eps_list, quantile, probs=0.25) < - eps_th
# ant = lapply(eps_list, quantile, probs=0.75) > eps_th
# 
# # check if there are any interactions that are both
# # synergetic and antagonistic at the same time
# both = (syn & ant)
# any(both)
# which(both)  
# 
# #====================OBTAIN THE RESTRICTED EPSILON DISTRIBUTION=================
# syn_th = 0.1
# # antagonism cutoff
# ant_th = 0.9
# 
# # restricted epsilon distribution for synergies
# eps_syn = eps_list
# eps_ant = eps_list
# for(c in 1:length(eps_list))
# {
#   comb_data = eps_list[[c]]
#   if(any(exp_fit_list[[c]] > syn_th)){
#     eps_syn[[c]] = comb_data[exp_fit_list[[c]] > syn_th]
#   } else eps_syn[[c]] = NA
#   
#   if(any(exp_fit_list[[c]] < ant_th)) {
#     eps_ant[[c]] = comb_data[exp_fit_list[[c]] < ant_th]
#   } else eps_ant[[c]] = NA
# }
# 
# 
# syn_rest = lapply(eps_syn, quantile, 
#                   probs=0.25, na.rm=T) < - eps_th
# syn_rest[is.na(syn_rest)] = FALSE
# sum(syn_rest, na.rm = T)
# 
# ant_rest = lapply(eps_ant, quantile, 
#                   probs=0.75, na.rm=T) > eps_th
# ant_rest[is.na(ant_rest)] = FALSE
# sum(ant_rest, na.rm = T)
# 
# 
# 
# #======================WRITE P VALUES========================================
# 
# pval.df = data.frame(pval = pval_wilc,
#                      synergistic = syn,
#                      antagonistic = ant,
#                      synrest = syn_rest,
#                      antrest = ant_rest,
#                      comb = names(eps_list))
# rownames(pval.df) = NULL
# 
# write.table(pval.df, 
#             file = paste0(wilc_dir, "pvals-all.txt"),
#             col.names = T, quote = F,
#             row.names = F)
# 
# 
# # detect strong drug-drug interactions based on the restricted epsilon distributions
# # with only relevant drug concentrations
# 
# eps.quantiles <- function(for_wilc) {
#   for_wilc = dplyr::group_by(for_wilc, unique) %>%
#     dplyr::summarise(eps.med = median(eps),
#                      eps.mean = mean(eps),
#                      eps.q1 = quantile(eps, 0.25),
#                      eps.q3 = quantile(eps, 0.75)) %>%
#     dplyr::rename(comb = unique)
#   
#   for_wilc
# }
# 
# eps_all = plyr::ldply(eps, data.frame, .id=NULL)
# eps_all_stats = eps.quantiles(eps_all)
# 
# write.table(eps_all_stats,
#             file = paste0(wilc_dir, "eps-stats-all.txt"),
#             col.names = T, quote = F,
#             row.names = F)
# 
# 
# bliss.vs.expfit = dplyr::rename(eps_all, comb = unique, expfit = multipl) %>%
#   dplyr::group_by(recconc, donid) %>% dplyr::mutate(eps = mean(eps),
#                                                     expfit = mean(expfit)) %>%
#   dplyr::inner_join(pval.df) %>%
#   dplyr::mutate(type = ifelse(synergistic & !antagonistic, "synergy",
#                               ifelse(antagonistic & !synergistic,
#                                      "antagonism", "no interaction")))
# 
# bliss.vs.expfit = dplyr::distinct(bliss.vs.expfit, recconc, donid, .keep_all=T  )
# 
# 
# # figure similar to Fig. Extendeed Data 3d
# batch_num = sort(gsub("batch", "", batches))
# 
# ext_data_3d = ggplot(bliss.vs.expfit, aes(x = expfit, y = eps, colour = factor(type))) +
#   geom_point() + theme_classic() + 
#   scale_color_manual(values = c("#C48B3F",
#                                 "#544E61",
#                                 "#85BAA1")) + 
#   labs(x = "Expected fitness",
#        y = "Bliss score",
#        colour = "Interaction type",
#        title = paste("Interaction scores vs expected fitness, batches",
#                      batch_num[1], "-", batch_num[length(batch_num)])) + 
#   coord_equal() +
#   theme(plot.title = element_text(hjust = 0.5))
# 
# ggsave(ext_data_3d, filename = paste0(figdir, "SF-3d.pdf"))
