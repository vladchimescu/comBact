# run Wilcoxon test on all batches
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
# effect size thresholds
strongsyn = params$strongsyn
strongant = params$strongant

strain = paste0("screen", argv[2])

script_dir = "~/Documents/embl/gitlab/drug_comb_screen/dataprep/"
if(!dir.exists(script_dir)) script_dir = "dataprep/"

figdir = "~/Documents/embl/gitlab/drug_comb_screen/figures/wilcoxon/"
if(!dir.exists(figdir)) {
  figdir = paste0("figures/", strain, "/signif/")
  dir.create(figdir)
}

# directory for writing epsilon's (interaction scores)
eps_dir = "~/Documents/embl/screen_K/interactions/"
if(!dir.exists(eps_dir)) {
  eps_dir = file.path(argv[1], "interactions/")
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
  medpconc = numeric()
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
    mdp = medpolish(eps_heatmat, maxiter = 1e3, na.rm = T)
    eps.medp = mdp$overall
    medp.conc = mdp$overall + ifelse(mdp$overall > 0, quantile(mdp$row, 0.75), quantile(mdp$row, 0.25))
    
    medpbliss = c(medpbliss, eps.medp)
    medpconc = c(medpconc, medp.conc)
    
  }
  
  data.frame(comb = combs,
             noise = cov.fit,
             medpbliss = medpbliss,
             medpconc = medpconc)
  
  
}

# read in the p-value files of drug combinations
# (that were computed in chunks)
pfiles = grep("^pvals-[0-9]+", list.files(wilc_dir), value = T)
pvals = lapply(pfiles, function(x) read.table(paste0(wilc_dir, x),
                                              header = T, stringsAsFactors = F))
pvals = plyr::ldply(pvals, data.frame)


#======================WRITE P VALUES========================================

# batches of the strain
batches = grep("batch", list.dirs(eps_dir, recursive = F, 
                                  full.names = F), value = T)

eps_df = get.bliss.data(eps_dir = eps_dir,
                        home_dir = home_dir)
ef_size = compute.covariates(eps_df = eps_df)

ef_quant = dplyr::group_by(eps_df, comb) %>%
  dplyr::summarise(eps.med = median(eps),
                   eps.mean = mean(eps),
                   eps.q1 = quantile(eps, 0.25),
                   eps.q3 = quantile(eps, 0.75))

eps.df = inner_join(ef_size, ef_quant)
pval.df = merge(pvals, eps.df, by = "comb")

pval.df = mutate(pval.df, efsize = ifelse(abs(eps.q1) > abs(eps.q3), eps.q1, eps.q3),
                 padj = p.adjust(pval, method = "BH"))

# add p-values of HTD screen in DSM20231
if(argv[2] == 'M') {
  pfiles = grep("^htd-pvals-[0-9]+", list.files(wilc_dir), value = T)
  pvals = lapply(pfiles, function(x) read.table(paste0(wilc_dir, x),
                                                header = T, stringsAsFactors = F))
  pvals = plyr::ldply(pvals, data.frame)
  
  pval.df.human = merge(pvals, eps.df, by = "comb")
  
  pval.df.human = mutate(pval.df.human, efsize = ifelse(abs(eps.q1) > abs(eps.q3), eps.q1, eps.q3),
                   padj = p.adjust(pval, method = "BH"))
  
  pval.df = bind_rows(pval.df, pval.df.human)
}

pval.df = dplyr::mutate(pval.df, d1 = gsub("(.+)_(.+)", "\\1", comb),
                  d2 = gsub("(.+)_(.+)", "\\2", comb))
druglegend = read.table(file.path(argv[1], "legend_gramnegpos.txt"),
                        sep = "\t", header = T,
                        stringsAsFactors = F)
pval.df = dplyr::mutate(pval.df,  d1 = plyr::mapvalues(d1, from = druglegend$code,
                                                          to = druglegend$Drug),
                        d2 = plyr::mapvalues(d2, from = druglegend$code,
                                               to = druglegend$Drug))
pval.df$comb = purrr::map2_chr(as.character(pval.df$d1),
                                as.character(pval.df$d2),
                                ~ paste(sort(c(.x, .y)), collapse = "_"))
pval.df = dplyr::select(pval.df, -c(d1, d2))
write.table(pval.df,
            file = paste0(wilc_dir, "pvals-all.txt"),
            col.names = T, quote = F,
            sep = "\t",
            row.names = F)


eps_df = dplyr::mutate(eps_df, d1 = gsub("(.+)_(.+)", "\\1", comb),
                        d2 = gsub("(.+)_(.+)", "\\2", comb))
eps_df = dplyr::mutate(eps_df,  d1 = plyr::mapvalues(d1, from = druglegend$code,
                                                       to = druglegend$Drug),
                        d2 = plyr::mapvalues(d2, from = druglegend$code,
                                             to = druglegend$Drug))
eps_df$comb = purrr::map2_chr(as.character(eps_df$d1),
                               as.character(eps_df$d2),
                               ~ paste(sort(c(.x, .y)), collapse = "_"))
eps_df = dplyr::select(eps_df, -c(d1, d2))


bliss.vs.expfit = dplyr::rename(eps_df, expfit = multipl) %>%
  dplyr::group_by(recconc, donconc) %>% dplyr::mutate(eps = mean(eps),
                                                    expfit = mean(expfit)) %>%
  dplyr::inner_join(pval.df) %>%
  dplyr::mutate(type = ifelse(padj <= 0.05,ifelse(efsize <= -strongsyn, "synergy",
                              ifelse(efsize >= strongant,
                                     "antagonism", "no interaction")), "no interaction"))

bliss.vs.expfit = dplyr::distinct(bliss.vs.expfit, recconc, donconc, .keep_all=T  )
saveRDS(bliss.vs.expfit, file = paste0(wilc_dir, strain, '-bliss-expfit.rds'))

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

ggsave(ext_data_3d_marg, filename = paste0(figdir, "SF-3d-wilcoxon.pdf"),
       width = 11, height = 11)

