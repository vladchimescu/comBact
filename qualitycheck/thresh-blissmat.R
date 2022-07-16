require(graphics)
require(stats)
require(Hmisc)
require(zoo)
require(gplots)
require(LSD)
require(corrplot)
library(igraph)
require(plotrix)
require(purrr)

argv = commandArgs(trailingOnly = TRUE)
# argv[1] = data/ directory
# argv[2] = strain
# argv = c("~/Documents/embl/screen_M", "M")

params = rjson::fromJSON(file = "params.json")

# trim AUC curves at this time point
trimtp = params$trimtp
# maximum normalized fitness (BUG control: fitness = 1)
maxfit = params$maxfit

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

batches = grep("batch", list.dirs(eps_dir, recursive = F, 
                                  full.names = F), value = T)


list_all_plates <- list.files(paste0(home_dir, batches))
all_plates = tools::file_path_sans_ext(list_all_plates)

plate_annot = data.frame(plate = all_plates, 
                         donor = unlist(lapply(strsplit(all_plates, "_"), function(x) x[1])),
                         biol_rep = unlist(lapply(strsplit(all_plates, "_"), function(x) x[4])),
                         batch = unlist(lapply(strsplit(all_plates, "_"), function(x) x[7])),
                         stringsAsFactors = F)

plate_annot = dplyr::mutate(plate_annot,
                            batch = sub("^0+", "", batch),
                            biol_rep = as.numeric(biol_rep))


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


# load the strongest associations
strong.inter = read.table(file = paste0(wilc_dir, "pvals-all.txt"),
                          stringsAsFactors = F, header = T)
strong.inter = dplyr::mutate(strong.inter, padj = p.adjust(pval, method = "BH"))
# subset only to interactions passing the strong threshold
strong.inter = dplyr::filter(strong.inter, abs(medpbliss) >= 0.05 & padj <= 0.05)
# only interactions close to the threshold
strong.inter = dplyr::filter(strong.inter, abs(medpbliss) < 0.07)
strong.inter = dplyr::mutate(strong.inter,  rec = gsub("(.+)(_.*)", "\\1", comb),
                      don = gsub("(.*_)(.+)", "\\2" , comb))
strong.inter = dplyr::mutate(strong.inter, 
                      type = ifelse(medpbliss < 0,
                                    "synergy",
                                    ifelse(medpbliss > 0,
                                           "antagonism", "none")))
strong.inter = dplyr::inner_join(strong.inter, eps_df)

library(dplyr)

strong.inter = mutate(strong.inter, 
                      recconc_num = gsub("([A-Z]+_)(.+)", "\\2", recconc),
                      donconc_num = gsub("([A-Z]+_)(.+)", "\\2", donconc))
strong.inter = distinct(inner_join(strong.inter, plate_annot))
strong.inter = arrange(strong.inter, padj, desc(abs(medpbliss)))

library(pheatmap)
combs = unique(strong.inter$comb)

pdf(paste0(figdir, "thresh", argv[2], "-blissmat.pdf"),
    width = 11, height = 10,
    onefile = T)
for(i in 1:length(combs)) {
  d1 = sub("([A-Z]+)_(.+)", "\\1", combs[i])
  
  
  # add biol_rep variable
  comb_df = filter(strong.inter, comb == combs[i]) %>%
    mutate(unique = paste(batch, biol_rep, Replicate, donor, sep = "_"))
  
  
  comb_df = comb_df %>% mutate(d1_conc = ifelse(donor == d1, 
                                                as.character(donconc),
                                                as.character(recconc)),
                               d2_conc = ifelse(donor == d1, as.character(recconc),
                                                as.character(donconc)))
  
  comb_eps_mat = comb_df %>%
    reshape2::dcast(d1_conc + d2_conc  ~ unique, mean, value.var = "eps")
  
  eps_heatmat = as.matrix(comb_eps_mat[,3:ncol(comb_eps_mat)])
  rownames(eps_heatmat) = paste(comb_eps_mat$d1_conc, comb_eps_mat$d2_conc, sep = "_")
  
  type = comb_df$type[1]
  eps.medp = comb_df$medpbliss[1]
  noise = comb_df$noise[1]
  
  paletteLength <- 200
  heat_scale <- colorRampPalette(c("#003836", "white", "#a87b00"))(paletteLength)
  
  myBreaks <- c(seq(-1, 0, length.out=ceiling(paletteLength/2) + 1), 
                seq(1/paletteLength, 1,
                    length.out=floor(paletteLength/2)))
  
  p = pheatmap(eps_heatmat,
               cluster_rows = F,
               cluster_cols = F,
               color = heat_scale,
               breaks = myBreaks,
               main = paste(combs[i], 
                            paste0(type, ";"), "medp(eps) =", 
                            format(eps.medp, digits = 2),
                            "cov(f_AB) =", format(noise, digits = 2)))
  
  grid::grid.draw(p$gtable)
  
}

dev.off()

