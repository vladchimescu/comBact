## integrated across-batch analysis of
## Bliss interaction scores and p-values
argv = commandArgs(trailingOnly = TRUE)
# argv[1] = data/ directory
# argv[2] = strain

strain = paste0("screen", argv[2])

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

data_dir = "~/Documents/embl/screen_K/"
if(!dir.exists(data_dir)) data_dir = paste0(argv[1], "/")

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

eps_df = rename(eps_df, comb = unique)

# load the strongest associations
strong.inter = read.table(file = paste0(wilc_dir, "strongassoc.tsv"),
                          stringsAsFactors = F, header = T)

eps_df_sign = inner_join(eps_df, strong.inter)

eps_df_sign = group_by(eps_df_sign, comb, batch) %>%
  mutate(eps.med = median(eps),
         eps.q1 = quantile(eps, 0.25),
         eps.q3 = quantile(eps, 0.75),
         iqr = IQR(eps),
         fit.comb.med = median(fit_comb)) %>%
  distinct(comb, batch, .keep_all = T)



# consistency
eps_cons = mutate(eps_df_sign, syn = eps.q1 < -0.1,
       ant = eps.q3 > 0.1) %>%
  group_by(comb, type) %>%
  summarise(synsum = sum(syn),
            antsum = sum(ant))

filter(eps_cons, type == 'synergy') %>%
  arrange(desc(synsum))


  eps.batch = inner_join(filter(eps_df_sign, batch == 1),
           filter(eps_df_sign, batch == 2),
           by = c("comb", "type"))

eps.batch = mutate(eps.batch,
                   concordant = ifelse(type == "synergy" ,
                                       ifelse(eps.q1.x < -0.1 & eps.q1.y < -0.1, "consistent", "inconsistent"),
                                       ifelse(type == "antagonism" & eps.q3.x > 0.1 & eps.q3.y > 0.1,
                                              "consistent", "inconsistent")))


eps_cor = ggplot(filter(eps.batch, type=='synergy'),
       aes(x = eps.q1.x,
           y = eps.q1.y)) + 
  geom_point(aes(color = factor(concordant))) + 
  xlim(c(-1,1)) + ylim(c(-1,1)) + 
  labs(x = "Bliss score (batch 1)",
       y = "Bliss score (batch 2)",
       title = "Bliss scores of replicated combinations") + 
  geom_abline(slope = 1, colour = 'red') + 
  theme_bw() + 
  coord_equal() + 
  theme(plot.title = element_text(hjust = 0.5))
ggsave(eps_cor, 
       filename = paste0(figdir, "eps-correlation12.pdf"))



