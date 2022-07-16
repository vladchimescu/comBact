# Script for multiple testing correction
# significant hits are annotated using the drug classification

argv = commandArgs(trailingOnly = TRUE)
# argv[1] = data/ directory
# argv[2] = strain

library(IHW)
library(dplyr)
library(ggplot2)
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

figdir = "~/Documents/embl/gitlab/drug_comb_screen/figures/wilcoxon/"
if(!dir.exists(figdir)) {
  figdir = paste0("figures/", strain, "/signif/")
}

home_dir = "~/Documents/embl/screen_K/screenK_annotated_reads/"
if(!dir.exists(home_dir)) {
  home_dir = file.path(argv[1], 
                       paste0(strain, "_annotated_reads/"))
}

wilc_dir = "~/Documents/embl/screen_K/wilcoxon/"
if(!dir.exists(wilc_dir)) {
  wilc_dir = file.path(argv[1], "wilcoxon/")
}

# directory for writing epsilon's (interaction scores)
eps_dir = "~/Documents/embl/screen_K/interactions/"
if(!dir.exists(eps_dir)) {
  eps_dir = file.path(argv[1], "interactions/")
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
  for_wilc$comb <- map2_chr(for_wilc$Drug, for_wilc$donor, ~ paste(sort(c(.x, .y)), collapse = "_"))
  
  for_wilc = subset(for_wilc, !Drug %in% c("ONE", "BUG", "NOT"))
  for_wilc = subset(for_wilc, !donor %in% c("ONE", "BUG", "NOT"))
  for_wilc
}

get.bliss.data <- function(eps_dir) {
  batches = grep("batch", list.dirs(eps_dir, recursive = F, 
                                    full.names = F), value = T)
  
  
  list_all_plates <- list.files(paste0(home_dir, batches))
  all_plates = tools::file_path_sans_ext(list_all_plates)
  
  plate_annot = data.frame(plate = all_plates, 
                           donor = unlist(lapply(strsplit(all_plates, "_"), function(x) x[1])),
                           donconc_num = unlist(lapply(strsplit(all_plates, "_"), function(x) x[2])),
                           biol_rep = unlist(lapply(strsplit(all_plates, "_"), function(x) x[4])),
                           batch = unlist(lapply(strsplit(all_plates, "_"), function(x) x[7])),
                           stringsAsFactors = F)
  
  plate_annot = dplyr::mutate(plate_annot,
                              batch = sub("^0+", "", batch))
  
  
  eps = lapply(batches, function(b) {
    eps = read.table(paste0(eps_dir, b, "/interactions.tsv"),
                     sep = "\t", header = T)
    eps = prep.interactions(eps)
    eps
  })
  names(eps) = batches
  eps_df = plyr::ldply(eps, data.frame, .id = "batch")
  # only batch number
  eps_df = dplyr::mutate(eps_df, batch = gsub("batch", "", batch))
  
  eps_df = mutate(eps_df, 
                        recconc_num = gsub("([A-Z]+_)(.+)", "\\2", recid),
                        donconc_num = gsub("([A-Z]+_)(.+)", "\\2", donconc))
  
  eps_df = eps_df %>% inner_join(plate_annot) %>%
    mutate(unique = paste(batch, biol_rep, Replicate, donor, sep = "_"))
  
  eps_df
}


compute.covariates <- function(eps_df) {
  combs = unique(eps_df$comb)
  medpbliss = numeric()
  cov.fit = numeric()
  
  for(i in 1:length(combs)) {
    d1 = sub("([A-Z]+)_(.+)", "\\1", combs[i])
    # estimate variance of individual combination fitness f(AB)
    # for each drugA-drugB-concA-concB combination and
    # compute  mean variance
    
    comb_df = filter(eps_df, comb == combs[i])
    
    comb_df = comb_df %>% mutate(d1_conc = ifelse(donor == d1, 
                                                  as.character(donconc),
                                                  as.character(recid)),
                                 d2_conc = ifelse(donor == d1, as.character(recid),
                                                  as.character(donconc)))
    
    comb_mat = comb_df %>%
      reshape2::dcast(d1_conc + d2_conc  ~ unique, mean, value.var = "fit_comb")
    
    heatmat = as.matrix(comb_mat[,3:ncol(comb_mat)])
    rownames(heatmat) = paste(comb_mat$d1_conc, comb_mat$d2_conc, sep="_")
    
    noise = mean(diag(var(t(heatmat))))
    cov.fit = c(cov.fit, noise)
    
    # estimate effect size as overall effect of median polished
    # matrix with (9 combs x n replicates)
    
    comb_eps_mat = comb_df %>%
      reshape2::dcast(d1_conc + d2_conc  ~ unique, mean, value.var = "eps")
    
    eps_heatmat = as.matrix(comb_eps_mat[,3:ncol(comb_eps_mat)])
    rownames(eps_heatmat) = paste(comb_eps_mat$d1_conc, comb_eps_mat$d2_conc, sep = "_")
    
    nrm =  ifelse(strain == "screenM", T, F)
    mdp = medpolish(eps_heatmat, maxiter = 50, na.rm = nrm)
    
    eps.medp = mdp$overall
    medpbliss = c(medpbliss, eps.medp)
  
  }
  
  data.frame(comb = combs,
             noise = cov.fit,
             medpbliss = medpbliss)
  
  
}

pvals =read.table(paste0(wilc_dir, "pvals-all.txt"),
                  header = T, stringsAsFactors = F)

# unfiltered p-value histogram
p_all = ggplot(pvals, aes(x = pval)) + 
  geom_histogram(bins = 50, fill = "#7EBEC5", colour = "black") + 
  theme_classic() + 
  labs(x = "P value (unadjusted)",
       y = "Count",
       title = "P value histogram of all interactions") + 
  theme(plot.title = element_text(hjust = 0.5))

ggsave(p_all, filename = paste0(figdir, "phist_all.pdf"))


eps_df = get.bliss.data(eps_dir = eps_dir)
eps_stats = compute.covariates(eps_df = eps_df)


# try stratifying by median fitness f_AB
# or by min(med(f_A), med(f_B))
# or by the dispersion of Bliss scores for a given drug combination\

# eps_stats= eps_df %>% group_by(comb) %>%
#   summarise(var.bliss = var(eps))

# eps_stats= eps_df %>% group_by(comb) %>%
#   summarise(med.multipl = median(multipl),
#             med.fit.ab = median(fit_comb),
#             min.fit.single = min(min(fit_rec_av, fit_don_av)))

# explore pval - covariate relation
# take medpolish for example for example
pval_covariate = inner_join(pvals, eps_stats)

# plot the "null" distribution of median polish Bliss scores
h = hist(pval_covariate$medpbliss, breaks = 50)

null.dist = ggplot(pval_covariate, aes(x = medpbliss)) +
  geom_histogram(bins=50) +
  theme_classic() + 
  geom_vline(xintercept = median(pval_covariate$medpbliss),
             colour = "red") + 
  geom_text(aes(label = format(median(pval_covariate$medpbliss), digits = 2),
                y=max(h$counts) + 50, x= median(pval_covariate$medpbliss))) + 
  labs(x = paste("Median polish Bliss scores in strain",
                 strain),
       y = "Count",
       title = "Distribution of median polish Bliss scores") +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(null.dist, filename = paste0(figdir, "null-distribution.pdf"))


eps_centered = pval_covariate$medpbliss - median(pval_covariate$medpbliss)
hist(eps_centered / sd(pval_covariate$medpbliss), 50)
hist(eps_centered, 50)

# 
# # only antagonistic candidates
# eps_ant = filter(pval_covariate, medpbliss > 0)$medpbliss
# hist(eps_ant, 50)



# for K
# ant_th = 0.03
ant_th = 0.04
p_ant = ggplot(filter(pval_covariate, medpbliss > ant_th), aes(x = pval)) + 
  geom_histogram(bins = 50, fill = "#7EBEC5", colour = "black") + 
  theme_classic() + 
  labs(x = "P value (unadjusted)",
       y = "Count",
       title = "P value histogram of antagonisms") + 
  theme(plot.title = element_text(hjust = 0.5))

ggsave(p_ant, filename = paste0(figdir, "phist_antagonisms.pdf"))


# for K
# syn_th = -0.005
syn_th = -0.01
p_syn = ggplot(filter(pval_covariate, medpbliss < -0.005), aes(x = pval)) + 
  geom_histogram(bins = 50, fill = "#7EBEC5", colour = "black") + 
  theme_classic() + 
  labs(x = "P value (unadjusted)",
       y = "Count",
       title = "P value histogram of synergies") + 
  theme(plot.title = element_text(hjust = 0.5))

ggsave(p_syn, filename = paste0(figdir, "phist_synergies.pdf"))



# ihwRes <- ihw(pval ~ min.fit.single,
#               data = pval_covariate,
#               nbins = 6,
#               alpha = 0.05,
#               seed = 3)
# 
# rejections(ihwRes)
# 
# plot(ihwRes, what = "decisionboundary") 
# 
# ihw.pval = adj_pvalues(ihwRes)
# 
# View(pval_covariate[which(ihw.pval <= 0.05),])


# the reason is quite trivial. Most of the interactions
# in this tail are "neutral" (neither synergistic nor antagonistic)

#--------------------STRONG ASSOCIATIONS BASED ON THE COMPLETE DISTRIBUTION------------
# get rid of them first

p_nonneut = ggplot(filter(pval_covariate, medpbliss < syn_th | medpbliss > ant_th), aes(x = pval)) + 
  geom_histogram(bins = 50, fill = "#F2CCC3", colour = "black") + 
  theme_classic() + 
  labs(x = "P value (unadjusted)",
       y = "Count",
       title = "P value histogram of non-neutral interactions") + 
  theme(plot.title = element_text(hjust = 0.5))

ggsave(p_nonneut, filename = paste0(figdir, "phist_nonneutral.pdf"))

pval.without.neutral = filter(pval_covariate,
                              medpbliss < syn_th | medpbliss > ant_th) %>%
  mutate(padj = p.adjust(pval, method = "BH"))


#=================SIGNIFICANT INTERACTIONS==============================
multpTest_cutoff = 0.05
signif = filter(pval.without.neutral, padj <= multpTest_cutoff)


# here 'don' and 'conc' are simply 'Drug 1' and 'Drug 2'
# i.e. there is no inherent order implied by (donor_conc)
signif1 = mutate(signif, rec = gsub("([A-Z]+)(_.*)", "\\1", comb),
                 don = gsub("(.*_)([A-Z]+)", "\\2" , comb))

signif2 = mutate(signif, don = gsub("([A-Z]+)(_.*)", "\\1", comb),
                 rec = gsub("(.*_)([A-Z]+)", "\\2" , comb))

signif = bind_rows(signif1, signif2)

# remove BUG donors and recipients
signif = filter(signif, rec != "BUG" & don != "BUG" & rec != "NOT" & don != "NOT")


signif = mutate(signif, 
                type = ifelse(synergistic & !antagonistic,
                              "synergy",
                              ifelse(antagonistic & !synergistic,
                                     "antagonism", "none")))

## remove "bimodal" interactions
signif = filter(signif, type != "none")

filler = expand.grid(don = signif$don, rec = signif$rec)
signif_hmap = right_join(signif, filler)

signif_hmap = distinct(signif_hmap, don, rec, type, .keep_all = T)

x.labels = sort(unique(signif_hmap$rec))
y.labels = sort(unique(signif_hmap$don))

inter_heatmap = ggplot(signif_hmap, 
                       aes(x = factor(rec, levels = x.labels),
                           y = factor(don, levels = y.labels),
                           fill = factor(type))) + 
  geom_tile(colour = "black") + theme_classic() + 
  scale_fill_manual(values = c("#C9788A", "#78A1A7"),
                    labels = c("Antagonism", "Synergy", "None"),
                    na.value = "#737373") + 
  labs(x = "Drug 1",
       y = "Drug 2",
       fill = "Interaction type",
       title = "Significant drug interactions") + 
  coord_equal() + 
  theme(axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.text.x = element_text(angle = 90),
        plot.title = element_text(hjust = 0.5))

ggsave(inter_heatmap,
       filename = paste0(figdir, "drugheatmap-all.pdf"),
       width = 10, height = 9)


# write strong significant associations
write.table(distinct(signif, comb, .keep_all = T),
            file = paste0(wilc_dir, "strongassoc.tsv"), quote = F,
            sep = "\t", row.names = F)

#----------------RESTRICTED ASSOCIATIONS------------------------
p_nonneut = ggplot(filter(pvals, synrest | antrest), aes(x = pval)) + 
  geom_histogram(bins = 50, fill = "#F2CCC3", colour = "black") + 
  theme_classic() + 
  labs(x = "P value (unadjusted)",
       y = "Count",
       title = "P value histogram of non-neutral interactions") + 
  theme(plot.title = element_text(hjust = 0.5))

ggsave(p_nonneut, filename = paste0(figdir, "phist_rest_nonneutral.pdf"))

pval.without.neutral = filter(pvals, synrest | antrest) %>%
  mutate(padj = p.adjust(pval, method = "BH"))


multpTest_cutoff = 0.05
signif.rest = filter(pval.without.neutral, padj <= multpTest_cutoff)

# here 'don' and 'conc' are simply 'Drug 1' and 'Drug 2'
# i.e. there is no inherent order implied by (donor_conc)
signif1 = mutate(signif.rest, rec = gsub("([A-Z]+)(_.*)", "\\1", comb),
                 don = gsub("(.*_)([A-Z]+)", "\\2" , comb))

signif2 = mutate(signif.rest, don = gsub("([A-Z]+)(_.*)", "\\1", comb),
                 rec = gsub("(.*_)([A-Z]+)", "\\2" , comb))

signif.rest = bind_rows(signif1, signif2)

# remove BUG donors and recipients
signif.rest = filter(signif.rest, rec != "BUG" & don != "BUG" & rec != "NOT" & don != "NOT")


signif.rest = mutate(signif.rest, 
                     type = ifelse(synrest & !antrest,
                                   "synergy",
                                   ifelse(antrest & !synrest,
                                          "antagonism", "none")))



rest_signif_inter = filter(signif.rest, comb %in% setdiff(signif.rest$comb, signif$comb))
rest_signif_inter = distinct(rest_signif_inter, comb, .keep_all = T)
stopifnot(nrow(rest_signif_inter) > 0)
write.table(rest_signif_inter, 
            file = paste0(wilc_dir, "restrictedsignif.tsv"), quote = F,
            sep = "\t", row.names = F)
