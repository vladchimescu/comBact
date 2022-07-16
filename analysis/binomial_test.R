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
library(dplyr)
library(ggExtra)

argv = commandArgs(trailingOnly = TRUE)
# argv[1] = data/ directory
# argv[2] = strain


params = rjson::fromJSON(file = "params.json")

# trim AUC curves at this time point
trimtp = params$trimtp
# maximum normalized fitness (BUG control: fitness = 1)
maxfit = params$maxfit
# effect size thresholds
strongsyn = params$strongsyn
strongant = params$strongant
weaksyn = params$weaksyn
weakant = params$weakant

strain = paste0("screen", argv[2])


script_dir = "~/Documents/embl/gitlab/drug_comb_screen/dataprep/"
if(!dir.exists(script_dir)) script_dir = "dataprep/"

figdir = "~/Documents/embl/gitlab/drug_comb_screen/figures/wilcoxon/"
if(!dir.exists(figdir)) {
  figdir = paste0("figures/", strain, "/signif/")
}

home_dir = "~/Documents/embl/screen_M/screenM_annotated_reads/"
if(!dir.exists(home_dir)) {
  home_dir = file.path(argv[1], 
                       paste0(strain, "_annotated_reads/"))
}

wilc_dir = "~/Documents/embl/screen_M/wilcoxon/"
if(!dir.exists(wilc_dir)) {
  wilc_dir = file.path(argv[1], "wilcoxon/")
}

# directory for writing epsilon's (interaction scores)
eps_dir = "~/Documents/embl/screen_M/interactions/"
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


get.bliss.data <- function(eps_dir, home_dir) {
  batches = grep("batch", list.dirs(eps_dir, recursive = F, 
                                    full.names = F), value = T)
  
  list_all_plates <- list.files(paste0(home_dir, batches))
  all_plates = tools::file_path_sans_ext(list_all_plates)
  
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
                  recconc_num = gsub("([A-Z]+_)(.+)", "\\2", recconc),
                  donconc_num = gsub("([A-Z]+_)(.+)", "\\2", donconc))
  
  eps_df = mutate(eps_df,
                  unique = paste(batch, biol_rep, Replicate, donor, sep = "_"))
  
  eps_df = distinct(eps_df, unique, combid, .keep_all = T)
  
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
                                                  as.character(recconc)),
                                 d2_conc = ifelse(donor == d1, as.character(recconc),
                                                  as.character(donconc)))
    
    comb_mat = comb_df %>%
      reshape2::dcast(d1_conc + d2_conc  ~ unique, mean, value.var = "fit_comb")
    
    heatmat = as.matrix(comb_mat[,3:ncol(comb_mat)])
    rownames(heatmat) = paste(comb_mat$d1_conc, comb_mat$d2_conc, sep="_")
    
    # estimate effect size as overall effect of median polished
    # matrix with (9 combs x n replicates)
    
    comb_eps_mat = comb_df %>%
      reshape2::dcast(d1_conc + d2_conc  ~ unique, mean, value.var = "eps")
    
    eps_heatmat = as.matrix(comb_eps_mat[,3:ncol(comb_eps_mat)])
    rownames(eps_heatmat) = paste(comb_eps_mat$d1_conc, comb_eps_mat$d2_conc, sep = "_")
    
    
    noise = sd(as.numeric(eps_heatmat), na.rm = T)
    cov.fit = c(cov.fit, noise)
    mdp = medpolish(eps_heatmat, maxiter = 50, na.rm = T)
    eps.medp = mdp$overall
    medpbliss = c(medpbliss, eps.medp)
    
  }
  
  data.frame(comb = combs,
             noise = cov.fit,
             medpbliss = medpbliss)
  
  
}

estimate_mode <- function(x) {
  d <- density(x)
  d$x[which.max(d$y)]
}
# batches of the strain
batches = grep("batch", list.dirs(eps_dir, recursive = F, 
                                  full.names = F), value = T)


eps_df = get.bliss.data(eps_dir = eps_dir,
                        home_dir = home_dir)

# nbins = 100
# h = hist(eps_df$eps, breaks = nbins)

null.dist = ggplot(eps_df, aes(x = eps)) +
  geom_density(fill = "blue", colour=NA, alpha=0.3) + 
  theme_classic() + 
  geom_vline(xintercept = median(eps_df$eps),
             colour = "red") + 
  geom_vline(xintercept = estimate_mode(eps_df$eps),
             colour = "black",
             linetype = "dotted") + 
  labs(x = paste("Bliss scores in strain",
                 strain),
       y = "Count",
       title = paste("Strain", argv[2],
                     "distribution of Bliss scores. Median:",
                     format(median(eps_df$eps), digits = 2),
                     "Mode:", format(estimate_mode(eps_df$eps),digits = 2))) +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(null.dist, filename = paste0(figdir, "null-eps-distribution.pdf"))

# median Bliss score
medbliss = estimate_mode(eps_df$eps)

library(dplyr)

eps_stats = eps_df %>% group_by(comb) %>%
  mutate(bernoulli = eps > medbliss) %>%
  dplyr::summarise(sumbern = sum(bernoulli),
            n=n()) %>%
  group_by(comb) %>%
  dplyr::mutate(pval = binom.test(x=sumbern,
                           n=n)$p.value)


eps_stats = ungroup(eps_stats) %>%
            mutate(padj = p.adjust(pval, method = "BH"))

# unfiltered p-value histogram
p_all = ggplot(eps_stats, aes(x = pval)) + 
  geom_histogram(bins = 30, fill = "#7EBEC5", colour = "black") + 
  theme_classic() + 
  labs(x = "P value (unadjusted)",
       y = "Count",
       title = paste("P value histogram of all interactions",
                     "strain", argv[2])) + 
  theme(plot.title = element_text(hjust = 0.5))

ggsave(p_all, filename = paste0(figdir, "phist_binom_all.pdf"))
  
#=================SIGNIFICANT INTERACTIONS==============================
multpTest_cutoff = 0.05
signif = filter(eps_stats, padj <= multpTest_cutoff)


# here 'don' and 'conc' are simply 'Drug 1' and 'Drug 2'
# i.e. there is no inherent order implied by (donor_conc)
signif1 = mutate(signif, rec = gsub("(.+)(_.*)", "\\1", comb),
                 don = gsub("(.*_)(.+)", "\\2" , comb))

signif2 = mutate(signif, don = gsub("(.+)(_.*)", "\\1", comb),
                 rec = gsub("(.*_)(.+)", "\\2" , comb))

signif = bind_rows(signif1, signif2)

# remove BUG donors and recipients
signif = filter(signif, rec != "BUG" & don != "BUG" & rec != "NOT" & don != "NOT")

ef_size = compute.covariates(eps_df = eps_df)

eps_out = inner_join(dplyr::select(eps_stats, comb, pval, padj),
                     dplyr::select(ef_size, comb, medpbliss), by = "comb")

#saveRDS(eps_out, file = paste0(argv[1], "/forshiny/bliss.rds"))

# quantiles of row medians
ef_quant = group_by(eps_df, recconc, donconc) %>%
  dplyr::mutate(rowmed = median(eps, na.rm = T)) %>%
  dplyr::group_by(comb) %>%
  dplyr::summarise(q1rowmed = quantile(rowmed, 0.25),
            q3rowmed = quantile(rowmed, 0.75))

signif = inner_join(signif, ef_size)
signif = inner_join(signif, ef_quant)

signif = mutate(signif, 
                type = ifelse(medpbliss <= -strongsyn & medpbliss < strongant,
                              "synergy",
                              ifelse(medpbliss >= strongant & medpbliss > -strongsyn,
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
        axis.text.x = element_text(vjust = 0.5, angle = 90),
        plot.title = element_text(hjust = 0.5))

ggsave(inter_heatmap,
       filename = paste0(figdir, "binomial-drugheatmap-all.pdf"),
       width = 10, height = 9)


# write strong significant associations
write.table(distinct(signif, comb, .keep_all = T),
            file = paste0(wilc_dir, "binomial-strongassoc.tsv"), quote = F,
            sep = "\t", row.names = F)


bliss.vs.expfit = dplyr::rename(eps_df, expfit = multipl) %>%
  dplyr::group_by(recconc, donconc) %>% dplyr::mutate(eps = mean(eps),
                                                    expfit = mean(expfit)) %>%
  dplyr::inner_join(eps_stats) %>%
  dplyr::inner_join(ef_size) %>%
  dplyr::mutate(type = ifelse(medpbliss < 0 & padj <= multpTest_cutoff, "synergy",
                              ifelse(medpbliss > 0 & padj <= multpTest_cutoff,
                                     "antagonism", "none")))

bliss.vs.expfit = dplyr::distinct(bliss.vs.expfit, recconc, donconc, .keep_all=T  )


# figure similar to Fig. Extendeed Data 3d
batch_num = sort(gsub("batch", "", batches))

ext_data_3d = ggplot(bliss.vs.expfit, aes(x = expfit, y = eps, colour = factor(type))) +
  geom_point() + theme_classic() +
  scale_color_manual(values = c("#C48B3F",
                                "#544E61",
                                "#85BAA1")) +
  labs(x = "Expected fitness",
       y = "Bliss score",
       colour = "Interaction type",
       title = paste("Interaction scores vs expected fitness, batches",
                     batch_num[1], "-", batch_num[length(batch_num)])) +
  coord_equal() +
  theme(plot.title = element_text(hjust = 0.5))


ext_data_3d_marg = ggMarginal(ext_data_3d, type = "histogram", bins=100)

ggsave(ext_data_3d_marg, filename = paste0(figdir, "SF-3d-binomial.pdf"),
       width = 11, height = 11)


# # add interactions based on restricting the distribution of
# # expected fitness
# syn_th = 0.2
# # antagonism cutoff
# ant_th = 0.8
# 
# synergies = filter(distinct(signif, comb, .keep_all = T), type == "synergy")$comb
# antagonisms = filter(distinct(signif, comb, .keep_all = T), type == "antagonism")$comb
# 
# syn_rest = compute.covariates(filter(eps_df, multipl >= syn_th))
# syn_rest_pval = filter(eps_df, multipl >= syn_th) %>%
#   group_by(comb) %>%
#   mutate(bernoulli = eps > medbliss) %>%
#   dplyr::summarise(sumbern = sum(bernoulli),
#                    n=n()) %>%
#   group_by(comb) %>%
#   dplyr::mutate(pval = binom.test(x=sumbern,
#                                   n=n)$p.value) %>%
#   ungroup() %>%  mutate(padj = p.adjust(pval, method = "BH"))
# syn_rest = inner_join(syn_rest, syn_rest_pval)
# syn_rest = mutate(syn_rest, 
#                   type = ifelse(padj <= multpTest_cutoff & medpbliss <= -strongsyn,
#                                 "synergy", NA))
# #syn_rest = filter(syn_rest, type == "synergy") %>% filter(!comb %in% synergies)
# 
# if(nrow(syn_rest)) {
#   write.table(syn_rest,
#               file = paste0(wilc_dir, "binomial-restricted-synergies.tsv"), quote = F,
#               sep = "\t", row.names = F)
# }
# 
# ant_rest = compute.covariates(filter(eps_df, multipl <= ant_th))
# ant_rest_pval = filter(eps_df, multipl <= ant_th) %>%
#   group_by(comb) %>%
#   mutate(bernoulli = eps > medbliss) %>%
#   dplyr::summarise(sumbern = sum(bernoulli),
#                    n=n()) %>%
#   group_by(comb) %>%
#   dplyr::mutate(pval = binom.test(x=sumbern,
#                                   n=n)$p.value) %>%
#   ungroup() %>%  mutate(padj = p.adjust(pval, method = "BH"))
# 
# ant_rest = inner_join(ant_rest, ant_rest_pval)
# ant_rest = mutate(ant_rest, 
#                   type = ifelse(padj <= multpTest_cutoff & medpbliss >= strongant,
#                                 "antagonism", NA))
# #ant_rest = filter(ant_rest, type == "antagonism") %>% filter(!comb %in% antagonisms)
# 
# if(nrow(ant_rest)) {
#   write.table(ant_rest,
#               file = paste0(wilc_dir, "binomial-restricted-antagonisms.tsv"), quote = F,
#               sep = "\t", row.names = F)
# }
# 
