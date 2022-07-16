## Check suspicious drug-drug interactions
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

figdir = "~/Documents/embl/gitlab/drug_comb_screen/figures/wilcoxon/"
if(!dir.exists(figdir)) {
  figdir = paste0("figures/", strain, "/signif/")
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
  for_wilc$unique <- map2_chr(for_wilc$Drug, for_wilc$donor, ~ paste(sort(c(.x, .y)), collapse = "_"))
  
  for_wilc = subset(for_wilc, !Drug %in% c("ONE", "BUG", "NOT"))
  for_wilc = subset(for_wilc, !donor %in% c("ONE", "BUG", "NOT"))
  for_wilc
}


batches = grep("batch", list.dirs(eps_dir, recursive = F, 
                                  full.names = F), value = T)
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


eps_distrib = ggplot(eps_df,
                     aes(x = eps, fill = factor(batch))) + 
  geom_density(alpha=0.6) + 
  theme_classic() + 
  labs(x = "Bliss score",
       y = "Density",
       fill = "Batch",
       title = paste("Bliss score distribution by batches", argv[2])) +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(filename = paste0(figdir, "eps-density-batches.pdf"), eps_distrib)


# load the strongest associations
strong.inter = read.table(file = paste0(wilc_dir, "strongassoc.tsv"),
                          stringsAsFactors = F, header = T)

eps_stats = read.table(file = paste0(wilc_dir, "eps-stats-all.txt"),
                 stringsAsFactors = F, header = T)

strong.inter = dplyr::inner_join(strong.inter, eps_stats)

strong.inter = dplyr::mutate(strong.inter, 
                      efsize = ifelse(type == "synergy",
                                      abs(eps.q1),
                                      eps.q3))

strong.inter = dplyr::arrange(strong.inter, padj, desc(efsize))


pdf(paste0(figdir, "comb-vs-multiplicative.pdf"), onefile = T)
for (i in 1:nrow(strong.inter)) {
 
  d_pair = sort(c(strong.inter$don[i],
                  strong.inter$rec[i]))
  
  check.comb = paste(d_pair[1], d_pair[2], sep="_")
  
  comb.df =dplyr::filter(eps_df, unique == check.comb)
  
  p = ggplot(comb.df,
         aes(x = multipl,
             y= fit_comb,
             colour = factor(batch)) ) + 
    geom_point(aes(shape = factor(Replicate))) + 
    geom_smooth(aes(x = multipl, y = fit_comb), colour = "blue",
                method = 'lm') + 
    geom_abline(slope = 1, colour = "black",
                linetype = "dashed") + 
    stat_smooth_func(geom="text",method='lm',
                     color = 'blue', hjust = 0, parse = T) +
    labs(x = "f(A)*f(B)",
         y = "f(AB)",
         colour = "Batch",
         shape = "Replicate",
         title = paste("Drug combination:", check.comb, strong.inter$type[i])) + 
    theme_bw() + 
    theme(plot.title = element_text(hjust = 0.5))
  print(p)
  
}
dev.off()


## empirical epsilon distributions for each significant combination
pdf(paste0(figdir, "comb-eps-distrib-with-repl.pdf"), 
    onefile = T,
    width = 10,
    height = 5)
for (i in 1:nrow(strong.inter)) {
  
  d_pair = sort(c(strong.inter$don[i],
                  strong.inter$rec[i]))
  
  check.comb = paste(d_pair[1], d_pair[2], sep="_")
  
  comb.df = dplyr::filter(eps_df, unique == check.comb)
  
  p1 = ggplot(comb.df,
              aes(x = eps, fill = factor(batch))) + 
    stat_density(geom = "line", aes(y = 0.02 * ..count.., colour = factor(batch))) +
    geom_dotplot(binwidth = 0.02, alpha=0.75) + 
    theme_classic() + 
    labs(x = "Bliss score",
         y = "Density / Dotplot count",
         fill = "Batch",
         colour = "Batch",
         title = paste("Density / Dotplot:", check.comb, strong.inter$type[i])) + 
    theme(plot.title = element_text(hjust = 0.5),
          aspect.ratio = 1)
  
  p2 = ggplot(comb.df,
             aes(x = eps, fill = factor(batch))) + 
    stat_ecdf(aes(colour = factor(batch))) + 
    theme_classic() + 
    labs(x = "Bliss score",
         y = "CDF(eps)",
         fill = "Batch",
         colour = "Batch",
         title = "ECDF") + 
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = "none",
          aspect.ratio = 1)
  

  cp = cowplot::plot_grid(p1, p2, scale = c(1, 0.9))
  print(cp)

}
dev.off()

# plot eps(AB) vs eps(BA)
pdf(paste0(figdir, "epsABepsBA.pdf"), 
    onefile = T)
for (i in 1:nrow(strong.inter)) {
  d_pair = sort(c(strong.inter$don[i],
                  strong.inter$rec[i]))
  
  check.comb = paste(d_pair[1], d_pair[2], sep="_")
  
  comb.df = dplyr::filter(eps_df, unique == check.comb)
  
  comb.repl = dplyr::inner_join(dplyr::filter(comb.df, Drug == d_pair[1]),
                                dplyr::filter(comb.df, Drug == d_pair[2]) %>%
                                  dplyr::rename(recconc = donid, donid = recconc),
                                by =  c("recconc", "donid", "batch"))
  if(nrow(comb.repl) == 0) {
   comb.repl = dplyr::inner_join(dplyr::filter(comb.df, Drug == d_pair[1]),
                                 dplyr::filter(comb.df, Drug == d_pair[2]) %>%
                                   dplyr::rename(recconc = donid, donid = recconc),
                                 by =  c("recconc", "donid")) %>%
     dplyr::mutate(batch = paste(batch.x, batch.y, sep="_"))
  }
  
  lim.max = max(c(abs(max(c(comb.repl$eps.x, comb.repl$eps.y))),
            abs(min(c(comb.repl$eps.x, comb.repl$eps.y)))))
  
  p = ggplot(comb.repl, aes(x = eps.x,
                        y = eps.y,
                        colour = factor(batch))) + 
    geom_point() + 
    labs(x = paste0("eps(", comb.repl$donrec.x[1], ")"),
         y = paste0("eps(", comb.repl$donrec.y[1], ")"),
         colour = "Batch",
         title = paste("Drug order effect on Bliss score of", check.comb, strong.inter$type[i])) + 
    geom_abline(slope = 1, colour = "red") +
    xlim(c(-lim.max, lim.max)) + ylim(-lim.max, lim.max) + 
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5),
          aspect.ratio = 1)
  print(p)
  
  
}

dev.off()


## now do the analysis for the restricted associations
strong.inter = read.table(file = paste0(wilc_dir, "restrictedsignif.tsv"),
                          stringsAsFactors = F, header = T)


strong.inter = dplyr::inner_join(strong.inter, eps_stats)

strong.inter = dplyr::mutate(strong.inter, 
                      efsize = ifelse(type == "synergy",
                                      abs(eps.q1),
                                      eps.q3))

strong.inter = dplyr::arrange(strong.inter, padj, desc(efsize))


pdf(paste0(figdir, "rest-comb-vs-multiplicative.pdf"), onefile = T)
for (i in 1:nrow(strong.inter)) {
  
  d_pair = sort(c(strong.inter$don[i],
                  strong.inter$rec[i]))
  
  check.comb = paste(d_pair[1], d_pair[2], sep="_")
  
  comb.df =dplyr::filter(eps_df, unique == check.comb)
  
  p = ggplot(comb.df,
             aes(x = multipl,
                 y= fit_comb,
                 colour = factor(batch)) ) + 
    geom_point(aes(shape = factor(Replicate))) + 
    geom_smooth(aes(x = multipl, y = fit_comb), colour = "blue",
                method = 'lm') + 
    geom_abline(slope = 1, colour = "black",
                linetype = "dashed") + 
    stat_smooth_func(geom="text",method='lm',
                     color = 'blue', hjust = 0, parse = T) +
    labs(x = "f(A)*f(B)",
         y = "f(AB)",
         colour = "Batch",
         shape = "Replicate",
         title = paste("Drug combination:", check.comb, strong.inter$type[i])) + 
    theme_bw() + 
    theme(plot.title = element_text(hjust = 0.5))
  print(p)
  
}
dev.off()


## empirical epsilon distributions for each significant combination
pdf(paste0(figdir, "rest-comb-eps-distrib-with-repl.pdf"), 
    onefile = T,
    width = 10,
    height = 5)
for (i in 1:nrow(strong.inter)) {
  
  d_pair = sort(c(strong.inter$don[i],
                  strong.inter$rec[i]))
  
  check.comb = paste(d_pair[1], d_pair[2], sep="_")
  
  comb.df = dplyr::filter(eps_df, unique == check.comb)
  
  p1 = ggplot(comb.df,
              aes(x = eps, fill = factor(batch))) + 
    stat_density(geom = "line", aes(y = 0.02 * ..count.., colour = factor(batch))) +
    geom_dotplot(binwidth = 0.02, alpha=0.75) + 
    theme_classic() + 
    labs(x = "Bliss score",
         y = "Density / Dotplot count",
         fill = "Batch",
         colour = "Batch",
         title = paste("Density / Dotplot:", check.comb, strong.inter$type[i])) + 
    theme(plot.title = element_text(hjust = 0.5),
          aspect.ratio = 1)
  
  p2 = ggplot(comb.df,
              aes(x = eps, fill = factor(batch))) + 
    stat_ecdf(aes(colour = factor(batch))) + 
    theme_classic() + 
    labs(x = "Bliss score",
         y = "CDF(eps)",
         fill = "Batch",
         colour = "Batch",
         title = "ECDF") + 
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = "none",
          aspect.ratio = 1)
  
  
  cp = cowplot::plot_grid(p1, p2, scale = c(1, 0.9))
  print(cp)
  
}
dev.off()

# plot eps(AB) vs eps(BA)
pdf(paste0(figdir, "rest-epsABepsBA.pdf"), 
    onefile = T)
for (i in 1:nrow(strong.inter)) {
  d_pair = sort(c(strong.inter$don[i],
                  strong.inter$rec[i]))
  
  check.comb = paste(d_pair[1], d_pair[2], sep="_")
  
  comb.df = dplyr::filter(eps_df, unique == check.comb)
  
  comb.repl = dplyr::inner_join(dplyr::filter(comb.df, Drug == d_pair[1]),
                                dplyr::filter(comb.df, Drug == d_pair[2]) %>%
                                  dplyr::rename(recconc = donid, donid = recconc),
                                by =  c("recconc", "donid", "batch"))
  if(nrow(comb.repl) == 0) {
    comb.repl = dplyr::inner_join(dplyr::filter(comb.df, Drug == d_pair[1]),
                                  dplyr::filter(comb.df, Drug == d_pair[2]) %>%
                                    dplyr::rename(recconc = donid, donid = recconc),
                                  by =  c("recconc", "donid")) %>%
      dplyr::mutate(batch = paste(batch.x, batch.y, sep="_"))
  }
  
  lim.max = max(c(abs(max(c(comb.repl$eps.x, comb.repl$eps.y))),
                  abs(min(c(comb.repl$eps.x, comb.repl$eps.y)))))
  
  p = ggplot(comb.repl, aes(x = eps.x,
                            y = eps.y,
                            colour = factor(batch))) + 
    geom_point() + 
    labs(x = paste0("eps(", comb.repl$donrec.x[1], ")"),
         y = paste0("eps(", comb.repl$donrec.y[1], ")"),
         colour = "Batch",
         title = paste("Drug order effect on Bliss score of", check.comb, strong.inter$type[i])) + 
    geom_abline(slope = 1, colour = "red") +
    xlim(c(-lim.max, lim.max)) + ylim(-lim.max, lim.max) + 
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5),
          aspect.ratio = 1)
  print(p)
  
  
}

dev.off()




